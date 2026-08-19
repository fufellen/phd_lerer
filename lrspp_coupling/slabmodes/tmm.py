"""Комплексные TM-моды произвольной слоистой структуры (метод матрицы переноса).

Задача. Плоскослоистая структура вдоль z, распространение вдоль x, поля не
зависят от y. Для TM-поляризации отличны от нуля H_y, E_x, E_z, а искомая
величина - комплексный эффективный показатель преломления n_eff = beta / k0.

Формулировка. Внутри слоя с проницаемостью eps и поперечной постоянной
затухания gamma = k0 * sqrt(n_eff^2 - eps) уравнение для H_y имеет вид
H_y'' = gamma^2 H_y. Непрерывны на границах H_y и тангенциальная E_x, поэтому
удобно вести две величины

    u = H_y,    v = (1 / eps) * dH_y/dz    (v пропорциональна E_x),

для которых система первого порядка принимает вид

    du/dz = eps * v,    dv/dz = (gamma^2 / eps) * u.

Матрица переноса через слой толщиной d:

    M = [[cosh(gamma*d),        sinh(gamma*d) / q],
         [q * sinh(gamma*d),    cosh(gamma*d)   ]],   q = gamma / eps.

Дисперсионное уравнение. В нижнем полупространстве оставляем только спадающую
вниз волну, что даёт стартовый вектор [1, q_bot]; в верхнем полупространстве
спадающая вверх волна требует v + q_top * u = 0. Отсюда невязка

    D(n_eff) = v_top + q_top * u_top,   [u_top, v_top] = M_total * [1, q_bot].

Проверка. Для трёх слоёв эта невязка с точностью до множителя q_1 совпадает с
формулой, уже провалидированной в edp/metal_strip_w800 против COMSOL:

    (q0 + q2) * q1 * cosh(g1 d) + (q0*q2 + q1^2) * sinh(g1 d).

Соответствие проверяется в scripts/selftest.py.

Соглашение о потерях. eps = (n + i*k)^2 с k > 0, поэтому у поглощающей среды
Im(eps) > 0, а у затухающей моды Im(n_eff) > 0. То же соглашение принято в
edp/metal_strip_w800/scripts/repeat_lerer_metal_strip.py.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Sequence

import numpy as np

try:  # scipy ускоряет корнепоиск, но не обязателен
    from scipy.optimize import root as _scipy_root
except ModuleNotFoundError:  # pragma: no cover
    _scipy_root = None

# numpy < 2.0 не знает trapezoid, numpy >= 2.0 объявил trapz устаревшим
_trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz")


def trapz(y, x):
    """Интегрирование по правилу трапеций, совместимое с numpy 1.x и 2.x."""
    return _trapz(y, x)


# ---------------------------------------------------------------- геометрия


@dataclass(frozen=True)
class Stack:
    """Слоистая структура: два полупространства и произвольное число слоёв.

    eps        - проницаемости всех слоёв снизу вверх, длина N >= 2;
    thickness  - толщины внутренних слоёв, длина N - 2, в микрометрах;
    names      - подписи слоёв для отчётов и графиков.
    """

    eps: tuple[complex, ...]
    thickness: tuple[float, ...]
    names: tuple[str, ...] = field(default=())

    def __post_init__(self) -> None:
        if len(self.eps) < 2:
            raise ValueError("нужно минимум два полупространства")
        if len(self.thickness) != len(self.eps) - 2:
            raise ValueError(
                f"толщин должно быть {len(self.eps) - 2}, получено {len(self.thickness)}"
            )
        if any(t <= 0 for t in self.thickness):
            raise ValueError("толщины внутренних слоёв должны быть положительными")

    @property
    def n_layers(self) -> int:
        return len(self.eps)

    @property
    def total_thickness(self) -> float:
        return float(sum(self.thickness))

    def interfaces(self) -> np.ndarray:
        """Координаты границ, z = 0 на нижней границе первого внутреннего слоя."""
        return np.concatenate(([0.0], np.cumsum(self.thickness)))

    def label(self, index: int) -> str:
        if index < len(self.names):
            return self.names[index]
        return f"слой {index}"

    def eps_at(self, z: np.ndarray) -> np.ndarray:
        """Проницаемость как функция z (для заливки областей на графиках)."""
        edges = self.interfaces()
        out = np.full(np.shape(z), self.eps[0], dtype=complex)
        out[z > edges[-1]] = self.eps[-1]
        for i, t in enumerate(self.thickness):
            inside = (z >= edges[i]) & (z <= edges[i + 1])
            out[inside] = self.eps[i + 1]
        return out


def decaying_sqrt(value: complex) -> complex:
    """Ветвь корня с Re >= 0: волна должна затухать от структуры, а не расти."""
    root = np.sqrt(complex(value))
    if root.real < 0:
        root = -root
    if abs(root.real) < 1e-14 and root.imag < 0:
        root = -root
    return complex(root)


def gammas(neff: complex, stack: Stack, k0: float) -> list[complex]:
    return [k0 * decaying_sqrt(neff * neff - e) for e in stack.eps]


# ------------------------------------------------------- дисперсионное уравнение


def tm_residual(neff: complex, stack: Stack, k0: float) -> complex:
    """Невязка дисперсионного уравнения; нули отвечают собственным TM-модам."""
    g = gammas(neff, stack, k0)
    q = [g[i] / stack.eps[i] for i in range(stack.n_layers)]

    u, v = 1.0 + 0j, q[0]  # спадающая вниз волна в нижнем полупространстве
    for i, d in enumerate(stack.thickness, start=1):
        gd = g[i] * d
        ch, sh = np.cosh(gd), np.sinh(gd)
        u, v = ch * u + sh / q[i] * v, q[i] * sh * u + ch * v

    return v + q[-1] * u


def tm_residual_3layer(neff: complex, stack: Stack, k0: float) -> complex:
    """Явная трёхслойная формула из edp/metal_strip_w800 - для перекрёстной проверки."""
    if stack.n_layers != 3:
        raise ValueError("формула определена только для трёх слоёв")
    g = gammas(neff, stack, k0)
    q = [g[i] / stack.eps[i] for i in range(3)]
    gd = g[1] * stack.thickness[0]
    return (q[0] + q[2]) * q[1] * np.cosh(gd) + (q[0] * q[2] + q[1] * q[1]) * np.sinh(gd)


# ---------------------------------------------------------------- корнепоиск


def solve_mode(
    stack: Stack,
    k0: float,
    guess: complex,
    residual: Callable[[complex, Stack, float], complex] = tm_residual,
    tol: float = 1e-12,
) -> complex:
    """Уточняет один комплексный корень из начального приближения."""
    if _scipy_root is not None:
        def wrapped(vec: np.ndarray) -> np.ndarray:
            value = residual(complex(vec[0], vec[1]), stack, k0)
            scale = max(abs(value), 1e-300)
            # нормировка убирает экспоненциальный разброс масштаба невязки
            value = value / scale ** 0.0
            return np.array([value.real, value.imag])

        sol = _scipy_root(wrapped, np.array([guess.real, guess.imag]), method="hybr", tol=tol)
        if sol.success:
            neff = complex(sol.x[0], sol.x[1])
            return neff.conjugate() if neff.imag < 0 else neff
    return _newton(residual, guess, stack, k0, tol)


def _newton(residual, guess: complex, stack: Stack, k0: float, tol: float) -> complex:
    """Ньютон с численным якобианом и дроблением шага; работает без scipy."""
    z = complex(guess)
    for _ in range(120):
        value = residual(z, stack, k0)
        f0 = np.array([value.real, value.imag], dtype=float)
        norm0 = float(np.linalg.norm(f0))
        if norm0 < tol:
            break

        hx = max(1e-9, abs(z.real) * 1e-8)
        hy = max(1e-9, abs(z.imag) * 1e-8)
        fx = residual(z + hx, stack, k0)
        fy = residual(z + 1j * hy, stack, k0)
        jac = np.array(
            [
                [(fx.real - value.real) / hx, (fy.real - value.real) / hy],
                [(fx.imag - value.imag) / hx, (fy.imag - value.imag) / hy],
            ]
        )
        try:
            step = np.linalg.solve(jac, -f0)
        except np.linalg.LinAlgError:
            step = np.linalg.lstsq(jac, -f0, rcond=None)[0]

        scale = 1.0
        for _ in range(20):
            candidate = z + scale * complex(step[0], step[1])
            if abs(residual(candidate, stack, k0)) < norm0:
                z = candidate
                break
            scale *= 0.5
        else:
            z += 0.05 * complex(step[0], step[1])

    return z.conjugate() if z.imag < 0 else z


def find_modes(
    stack: Stack,
    k0: float,
    guesses: Sequence[complex],
    n_min: float | None = None,
    n_max: float | None = None,
    unique_tol: float = 1e-7,
    residual_tol: float = 1e-6,
) -> list[complex]:
    """Прогоняет набор приближений и возвращает различные принятые корни.

    Отбираются только связанные моды: Re(n_eff) выше показателя обкладок и ниже
    заданного потолка, невязка мала, мнимая часть отвечает затуханию.
    """
    accepted: list[complex] = []
    for guess in guesses:
        try:
            neff = solve_mode(stack, k0, complex(guess))
        except Exception:
            continue
        if not np.isfinite(neff.real) or not np.isfinite(neff.imag):
            continue

        # невязку нормируем на её масштаб, иначе абсолютный порог бессмыслен
        scale = abs(tm_residual(neff * (1 + 1e-3) + 1e-3, stack, k0)) + 1e-300
        if abs(tm_residual(neff, stack, k0)) / scale > residual_tol:
            continue
        if n_min is not None and neff.real <= n_min:
            continue
        if n_max is not None and neff.real >= n_max:
            continue
        if any(abs(neff - old) < unique_tol for old in accepted):
            continue
        accepted.append(neff)

    accepted.sort(key=lambda z: -z.real)
    return accepted


# ------------------------------------------------------------------- поля


def mode_fields(
    neff: complex, stack: Stack, k0: float, z: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Восстанавливает H_y, E_x, E_z вдоль z для найденной моды.

    В полупространствах поле задаётся аналитической спадающей экспонентой, а не
    результатом численного продолжения: растущая экспонента принудительно
    обнуляется. Иначе малая невязка корня, умноженная на экспоненту, создаёт
    ложный хвост и портит интегралы перекрытия и оценку размера моды.
    """
    z = np.asarray(z, dtype=float)
    g = gammas(neff, stack, k0)
    q = [g[i] / stack.eps[i] for i in range(stack.n_layers)]
    edges = stack.interfaces()

    hy = np.zeros_like(z, dtype=complex)
    ex = np.zeros_like(z, dtype=complex)

    # значения (u, v) на нижней границе каждого внутреннего слоя
    u, v = 1.0 + 0j, q[0]
    nodes = [(u, v)]
    for i, d in enumerate(stack.thickness, start=1):
        gd = g[i] * d
        ch, sh = np.cosh(gd), np.sinh(gd)
        u, v = ch * u + sh / q[i] * v, q[i] * sh * u + ch * v
        nodes.append((u, v))

    below = z < edges[0]
    hy[below] = nodes[0][0] * np.exp(g[0] * (z[below] - edges[0]))
    ex[below] = q[0] * hy[below]

    above = z > edges[-1]
    hy[above] = nodes[-1][0] * np.exp(-g[-1] * (z[above] - edges[-1]))
    ex[above] = -q[-1] * hy[above]

    for i, d in enumerate(stack.thickness, start=1):
        z0 = edges[i - 1]
        inside = (z >= z0) & (z <= edges[i])
        if not np.any(inside):
            continue
        u0, v0 = nodes[i - 1]
        gz = g[i] * (z[inside] - z0)
        ch, sh = np.cosh(gz), np.sinh(gz)
        hy[inside] = ch * u0 + sh / q[i] * v0
        ex[inside] = q[i] * sh * u0 + ch * v0

    eps_z = stack.eps_at(z)
    ex = ex / 1j  # теперь E_x с точностью до положительного множителя 1/(omega eps0)
    ez = -neff * hy / eps_z  # E_z с тем же множителем
    return hy, ex, ez


def power_flux(hy: np.ndarray, neff: complex, eps_z: np.ndarray, z: np.ndarray) -> float:
    """Продольный поток мощности моды.

    Распространение идёт вдоль x, поэтому переносимая мощность определяется
    компонентой S_x = Re(-E_z H_y*) / 2, а не поперечной S_z. Подстановка
    E_z = -(beta / (omega eps0 eps)) H_y даёт

        P = (1 / (2 omega eps0)) * int Re(beta / eps) |H_y|^2 dz,

    что возвращается здесь с точностью до положительного множителя. Внутри
    металла Re(eps) < 0, поэтому вклад отрицателен: это известный обратный
    поток мощности в металле, а не ошибка расчёта.
    """
    return float(_trapz(np.real(neff / eps_z) * np.abs(hy) ** 2, z))


def overlap_power(
    mode_a: tuple[np.ndarray, complex, np.ndarray],
    mode_b: tuple[np.ndarray, complex, np.ndarray],
    z: np.ndarray,
) -> float:
    """Доля мощности, переходящая из моды a в моду b на стыке.

    Взаимная мощностная проекция для слабо отражающего стыка, записанная через
    H_y и eps (для TM это эквивалентно интегралу от E x H* по сечению):

        eta = |int (n_a/eps + n_b*/eps*) H_a H_b* dz|^2
              / (4 * int Re(n_a/eps)|H_a|^2 dz * int Re(n_b/eps)|H_b|^2 dz).

    Нормировка выбрана так, что перекрытие моды с самой собой равно единице.

    Аргументы: (H_y, n_eff, eps(z)) для каждой моды на общей сетке z.
    """
    ha, na, eps_a = mode_a
    hb, nb, eps_b = mode_b
    cross = _trapz((na / eps_a + np.conj(nb / eps_b)) * ha * np.conj(hb), z)
    pa = power_flux(ha, na, eps_a, z)
    pb = power_flux(hb, nb, eps_b, z)
    if pa <= 0 or pb <= 0:
        return float("nan")
    return float(abs(cross) ** 2 / (4.0 * pa * pb))


def orthogonality_residual(
    h1: np.ndarray, h2: np.ndarray, eps_z: np.ndarray, z: np.ndarray
) -> float:
    """Невязка несопряжённой ортогональности TM-мод одной структуры.

    Для мод с потерями задача неэрмитова, и обращается в нуль не сопряжённое
    перекрытие, а несопряжённое:

        int H_m H_n / eps dz = 0,   m != n.

    Возвращается модуль этого интеграла, нормированный на собственные нормы,
    поэтому для разных мод результат должен быть близок к нулю.
    """
    cross = _trapz(h1 * h2 / eps_z, z)
    n1 = _trapz(h1 * h1 / eps_z, z)
    n2 = _trapz(h2 * h2 / eps_z, z)
    denom = np.sqrt(abs(n1) * abs(n2))
    if denom == 0:
        return float("nan")
    return float(abs(cross) / denom)


# ------------------------------------------------------------- метрики моды


def propagation_loss_db_per_cm(neff: complex, lambda_um: float) -> float:
    """Погонные потери по мощности, дБ/см."""
    k0 = 2.0 * np.pi / lambda_um
    alpha_per_um = 2.0 * k0 * abs(neff.imag)  # по мощности
    return float(alpha_per_um * 1e4 * 10.0 / np.log(10.0))


def propagation_length_um(neff: complex, lambda_um: float) -> float:
    """Длина распространения по мощности, L = 1 / (2 Im beta)."""
    k0 = 2.0 * np.pi / lambda_um
    if neff.imag == 0:
        return float("inf")
    return float(1.0 / (2.0 * k0 * abs(neff.imag)))


def mode_width_1e(hy: np.ndarray, z: np.ndarray) -> float:
    """Полная ширина профиля |H_y| по уровню 1/e от максимума."""
    amp = np.abs(hy)
    peak = amp.max()
    if peak <= 0:
        return float("nan")
    mask = amp >= peak / np.e
    if not np.any(mask):
        return float("nan")
    return float(z[mask][-1] - z[mask][0])
