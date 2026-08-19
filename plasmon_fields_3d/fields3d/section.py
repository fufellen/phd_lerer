"""Поле TM-моды в плоскости распространения и проверка его уравнениями Максвелла.

Постановка. Структура слоиста вдоль z и однородна вдоль y, распространение идёт
вдоль x. У TM-поляризации отличны от нуля H_y, E_x и E_z, поэтому вся картина
задаётся в плоскости (x, z), а третье измерение получается вытягиванием вдоль y.
Так и строятся трёхмерные сцены: хранится сечение, а объём собирается только на
этапе отрисовки. Это экономит память на два порядка и не даёт выдать за
трёхмерный расчёт то, что им не является.

Нормировка наследуется от `lrspp_coupling/slabmodes/tmm.py`: длины в мкм,

    E~ = (omega eps0 / 1e6) E,   H~ = H_y,

и в этих переменных уравнения Максвелла для монохроматического поля вида
exp(-i omega t) принимают вид

    rot E~ = i k0^2 H~,     rot H~ = -i eps E~,

где k0 в мкм^-1, а производные берутся по мкм. Обе невязки считает
`maxwell_residual`, и это главная проверка сцены: она подтверждает, что
нарисовано решение уравнений Максвелла, а не правдоподобная картинка.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

_LRSPP = Path(__file__).resolve().parents[2] / "lrspp_coupling"
if str(_LRSPP) not in sys.path:
    sys.path.insert(0, str(_LRSPP))

from slabmodes import Stack, mode_fields  # noqa: E402


@dataclass
class FieldSection:
    """Комплексные фазоры TM-поля на сетке (x, z), мкм.

    Массивы полей имеют форму (nx, nz). Поле в момент времени t получается как
    Re(F * exp(-i omega t)); за это отвечает `instantaneous`.
    """

    x: np.ndarray
    z: np.ndarray
    ex: np.ndarray
    ez: np.ndarray
    hy: np.ndarray
    eps: np.ndarray
    lambda_um: float
    title: str = ""
    notes: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        shape = (self.x.size, self.z.size)
        for name in ("ex", "ez", "hy", "eps"):
            arr = getattr(self, name)
            if arr.shape != shape:
                raise ValueError(f"{name}: ожидалась форма {shape}, получена {arr.shape}")

    @property
    def k0(self) -> float:
        return 2.0 * np.pi / self.lambda_um

    # ------------------------------------------------------------ построение

    @classmethod
    def from_mode(
        cls,
        stack: Stack,
        neff: complex,
        lambda_um: float,
        x: np.ndarray,
        z: np.ndarray,
        amplitude: complex = 1.0,
        z_offset_um: float = 0.0,
        title: str = "",
        notes: dict | None = None,
    ) -> "FieldSection":
        """Одна планарная мода, вытянутая вдоль x множителем exp(i beta x).

        `z_offset_um` сдвигает стопку относительно сетки: собственная система
        координат `Stack` начинается на нижней границе первого внутреннего слоя,
        а сцене обычно удобно центрировать структуру около нуля.
        """
        k0 = 2.0 * np.pi / lambda_um
        z_local = np.asarray(z, dtype=float) - z_offset_um
        hy_z, ex_z, ez_z = mode_fields(neff, stack, k0, z_local)
        eps_z = stack.eps_at(z_local)

        phase = np.exp(1j * k0 * neff * np.asarray(x, dtype=float))
        outer = amplitude * phase[:, None]
        return cls(
            x=np.asarray(x, dtype=float),
            z=np.asarray(z, dtype=float),
            ex=outer * ex_z[None, :],
            ez=outer * ez_z[None, :],
            hy=outer * hy_z[None, :],
            eps=np.broadcast_to(eps_z[None, :], (np.size(x), np.size(z))).copy(),
            lambda_um=lambda_um,
            title=title,
            notes=dict(notes or {}),
        )

    @classmethod
    def superpose(cls, parts: list["FieldSection"], title: str = "") -> "FieldSection":
        """Сумма нескольких мод одной и той же структуры.

        Линейность уравнений Максвелла делает сумму решений решением, поэтому
        суперпозиция супермод - точное поле связанной пары, а не приближение.
        Требуется совпадение сеток и распределения eps.
        """
        if not parts:
            raise ValueError("нечего складывать")
        base = parts[0]
        for p in parts[1:]:
            if not np.array_equal(p.x, base.x) or not np.array_equal(p.z, base.z):
                raise ValueError("сетки слагаемых не совпадают")
            if not np.allclose(p.eps, base.eps):
                raise ValueError("слагаемые описывают разные структуры")
        return cls(
            x=base.x,
            z=base.z,
            ex=sum(p.ex for p in parts),
            ez=sum(p.ez for p in parts),
            hy=sum(p.hy for p in parts),
            eps=base.eps,
            lambda_um=base.lambda_um,
            title=title or base.title,
            notes=dict(base.notes),
        )

    # -------------------------------------------------------------- величины

    def instantaneous(self, phase_rad: float = 0.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Мгновенные вещественные поля Re(F exp(-i omega t)) при omega t = phase."""
        rot = np.exp(-1j * phase_rad)
        return (
            np.real(self.ex * rot),
            np.real(self.ez * rot),
            np.real(self.hy * rot),
        )

    def intensity(self) -> np.ndarray:
        """|E|^2 = |E_x|^2 + |E_z|^2 в нормировке модуля."""
        return np.abs(self.ex) ** 2 + np.abs(self.ez) ** 2

    def magnitude_e(self) -> np.ndarray:
        return np.sqrt(self.intensity())

    def poynting_x(self) -> np.ndarray:
        """Продольная плотность потока мощности, S_x = Re(-E_z H_y*) / 2.

        Именно эта компонента переносит мощность вдоль волновода. Внутри металла
        Re(eps) < 0 и вклад отрицателен - это известный обратный поток, а не
        ошибка расчёта.
        """
        return 0.5 * np.real(-self.ez * np.conj(self.hy))

    def power_along_x(self) -> np.ndarray:
        """Полная мощность в сечении как функция x."""
        return np.trapezoid(self.poynting_x(), self.z, axis=1)

    def surface_charge(self, z_interface_um: float) -> np.ndarray:
        """Поверхностная плотность заряда на границе, sigma ~ [eps E_z].

        Скачок нормальной компоненты D на границе двух сред равен поверхностной
        плотности заряда. Для ППП это и есть тот колеблющийся заряд, которым
        плазмон отличается от обычной направляемой волны.
        """
        j = int(np.argmin(np.abs(self.z - z_interface_um)))
        below = max(j - 1, 0)
        above = min(j + 1, self.z.size - 1)
        return np.real(self.eps[:, above] * self.ez[:, above] - self.eps[:, below] * self.ez[:, below])

    # ------------------------------------------------------------- проверка

    def _uniform_mask(self) -> np.ndarray:
        """Точки, чей трёхточечный шаблон целиком лежит внутри одного слоя.

        На границе сред E_z разрывна, а производная E_x по z терпит излом,
        поэтому конечная разность там неприменима и такие точки исключаются.
        """
        same_z = np.ones_like(self.eps, dtype=bool)
        same_z[:, 1:-1] = (self.eps[:, :-2] == self.eps[:, 1:-1]) & (
            self.eps[:, 2:] == self.eps[:, 1:-1]
        )
        same_z[:, 0] = False
        same_z[:, -1] = False

        same_x = np.ones_like(self.eps, dtype=bool)
        same_x[1:-1, :] = (self.eps[:-2, :] == self.eps[1:-1, :]) & (
            self.eps[2:, :] == self.eps[1:-1, :]
        )
        same_x[0, :] = False
        same_x[-1, :] = False
        return same_z & same_x

    def maxwell_residual(self) -> dict[str, float]:
        """Относительная невязка уравнений Максвелла на сетке сцены.

        Считаются три числа: закон Фарадея (rot E)_y = i k0^2 H_y и две
        компоненты закона Ампера (rot H)_x = -i eps E_x, (rot H)_z = -i eps E_z.
        Все производные - центральные разности, поэтому ожидаемый уровень
        невязки задан шагом сетки, а не физикой.
        """
        mask = self._uniform_mask()
        if not np.any(mask):
            return {"faraday": float("nan"), "ampere_x": float("nan"), "ampere_z": float("nan")}

        dex_dz = np.gradient(self.ex, self.z, axis=1)
        dez_dx = np.gradient(self.ez, self.x, axis=0)
        dhy_dz = np.gradient(self.hy, self.z, axis=1)
        dhy_dx = np.gradient(self.hy, self.x, axis=0)

        k0 = self.k0
        floor = k0 * np.abs(self.hy)

        rot_e_y = dex_dz - dez_dx
        scale_f = np.maximum.reduce([np.abs(dex_dz), np.abs(dez_dx), k0 * floor])
        res_f = np.abs(rot_e_y - 1j * k0 * k0 * self.hy) / np.maximum(scale_f, 1e-300)

        need_x = -1j * self.eps * self.ex
        scale_x = np.maximum.reduce([np.abs(dhy_dz), np.abs(need_x), floor])
        res_x = np.abs(-dhy_dz - need_x) / np.maximum(scale_x, 1e-300)

        need_z = -1j * self.eps * self.ez
        scale_z = np.maximum.reduce([np.abs(dhy_dx), np.abs(need_z), floor])
        res_z = np.abs(dhy_dx - need_z) / np.maximum(scale_z, 1e-300)

        return {
            "faraday": float(res_f[mask].max()),
            "ampere_x": float(res_x[mask].max()),
            "ampere_z": float(res_z[mask].max()),
            "points_checked": int(mask.sum()),
        }

    def subsample(self, every_x: int = 1, every_z: int = 1) -> "FieldSection":
        """Прореженная копия для отрисовки.

        Расчётная сетка и сетка картинки - разные вещи. Для проверки физики
        нужен мелкий шаг, а объёмная отрисовка миллионов узлов не только
        медленна, но и просто не отображается: у сцены с трассой в сотни
        микрометров сетка доходит до миллионов точек, и объёмный рендер молча
        выдаёт пустой кадр. Прореживание делается после всех проверок.
        """
        sx = max(int(every_x), 1)
        sz = max(int(every_z), 1)
        return FieldSection(
            x=self.x[::sx],
            z=self.z[::sz],
            ex=self.ex[::sx, ::sz],
            ez=self.ez[::sx, ::sz],
            hy=self.hy[::sx, ::sz],
            eps=self.eps[::sx, ::sz],
            lambda_um=self.lambda_um,
            title=self.title,
            notes=dict(self.notes),
        )

    def grid_quality(self) -> dict[str, float]:
        """Насколько мелко сетка разрешает самый быстрый масштаб поля.

        Ограничивает точность любой конечной разности на этой сетке: невязка
        центральной разности растёт как (gamma*dz)^2 / 6, где gamma - наибольшая
        поперечная постоянная затухания среди слоёв. В плазмонных задачах её
        задаёт металл, где поле спадает за десятки нанометров, поэтому именно
        толщина скин-слоя, а не длина волны, диктует шаг сетки.
        """
        dz = float(self.z[1] - self.z[0])
        # sqrt(n_eff^2 - eps) для металла оценивается сверху величиной sqrt(|eps|)
        gamma_max = float(max(self.k0 * np.sqrt(abs(complex(e))) for e in np.unique(self.eps)))
        return {"dz_um": dz, "gamma_max_per_um": gamma_max, "gamma_dz": gamma_max * dz}

    def energy_decay_check(self, neff: complex) -> dict[str, float]:
        """Сверяет спад мощности вдоль x с 2 Im(beta), полученным из n_eff.

        Проверка независима от восстановления поля: поток считается интегралом
        по сечению, а эталон - по собственному значению.
        """
        p = self.power_along_x()
        good = p > 0
        if good.sum() < 10:
            return {"fitted_alpha": float("nan"), "expected_alpha": float("nan"), "rel_error": float("nan")}
        slope = np.polyfit(self.x[good], np.log(p[good]), 1)[0]
        expected = -2.0 * self.k0 * abs(neff.imag)
        rel = abs(slope - expected) / max(abs(expected), 1e-300)
        return {"fitted_alpha": float(slope), "expected_alpha": float(expected), "rel_error": float(rel)}


def residual_convergence(build) -> dict[str, float]:
    """Невязка Максвелла на двух сетках и наблюдаемый порядок сходимости.

    Зачем именно так. Абсолютная величина невязки на сетке сцены сама по себе
    ничего не доказывает: она задана шагом, а шаг у плазмонной задачи диктует
    скин-слой металла в десятки нанометров. Порог пришлось бы подбирать под
    каждую сцену, и он маскировал бы ошибку в самом построении поля.

    Осмысленная проверка - поведение при измельчении. Если поле построено верно,
    остаётся только ошибка центральной разности, и при удвоении числа узлов
    невязка падает вчетверо. Ошибка в формуле поля так себя не ведёт: она даёт
    невязку, не зависящую от шага. Проверяются два числа - наблюдаемый порядок
    около 2 и абсолютный уровень на мелкой сетке.

    Измельчать нужно обе оси сразу. Невязка - максимум по всем уравнениям, и
    производные берутся и по x, и по z, поэтому при измельчении только одной оси
    вклад второй остаётся постоянным и наблюдаемый порядок падает до нуля. Это
    не признак ошибки в поле, а признак неверно поставленной проверки.

    `build(refine)` должна вернуть FieldSection на сетке, у которой обе оси
    измельчены в `refine` раз относительно базовой.
    """
    coarse = build(1)
    fine = build(2)
    r_c = coarse.maxwell_residual()
    r_f = fine.maxwell_residual()

    worst_c = max(r_c["faraday"], r_c["ampere_x"], r_c["ampere_z"])
    worst_f = max(r_f["faraday"], r_f["ampere_x"], r_f["ampere_z"])
    order = float(np.log2(worst_c / worst_f)) if worst_f > 0 else float("inf")
    return {
        "residual_coarse": float(worst_c),
        "residual_fine": float(worst_f),
        "observed_order": order,
        "gamma_dz_fine": fine.grid_quality()["gamma_dz"],
        "points_checked": float(r_f["points_checked"]),
    }
