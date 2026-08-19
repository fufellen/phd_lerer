"""Продольный расчёт плавного перехода по локальным модам.

Задача. Поперечное сечение медленно меняется вдоль оси (например, сужается
ширина полоски). Нужно посчитать, сколько мощности дойдёт до конца в рабочей
моде и куда денется остальное.

Метод. Переход разбивается на короткие участки, внутри каждого сечение
считается постоянным. Внутри участка каждая локальная мода набирает фазу и
затухание по exp(i beta L); на стыке двух участков поля сшиваются через
интегралы перекрытия локальных мод.

Три канала потерь разделяются явно:
  - поглощение: из Im(beta) локальных мод;
  - преобразование: перекачка в другие связанные моды на стыках;
  - излучение: доля мощности, не захваченная ни одной связанной модой.

Важное ограничение. Базис здесь состоит только из связанных мод, поэтому
излучение оценивается как дефицит проекции. Для ступенчатой аппроксимации
гладкого перехода этот дефицит убывает при измельчении шага как 1/N и в пределе
стремится к нулю, то есть расчёт сходится к адиабатическому ответу и НЕ даёт
сходящейся оценки излучения. Чтобы не выдавать артефакт за физику, модуль
считает дефицит только диагностикой сходимости, а неадиабатичность оценивает
отдельно - критерием Лава для локального угла сужения.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np

from .tmm import trapz


@dataclass(frozen=True)
class TaperProfile:
    """Продольный профиль перехода: ширина как функция координаты."""

    length_um: float
    width_start_um: float
    width_end_um: float

    def width_at(self, z: np.ndarray) -> np.ndarray:
        t = np.clip(np.asarray(z, dtype=float) / self.length_um, 0.0, 1.0)
        return self.width_start_um + (self.width_end_um - self.width_start_um) * t

    @property
    def half_angle_rad(self) -> float:
        """Полуугол сужения: угол одной кромки к оси."""
        return float(np.arctan(0.5 * abs(self.width_start_um - self.width_end_um) / self.length_um))

    @property
    def half_angle_deg(self) -> float:
        return float(np.degrees(self.half_angle_rad))


def adiabatic_transmission(
    profile: TaperProfile,
    neff_of_width: Callable[[float], complex],
    lambda_um: float,
    steps: int = 400,
) -> dict[str, float]:
    """Адиабатический предел: мода следует за сечением, теряя только на поглощение.

    Возвращает долю прошедшей мощности и накопленные потери в децибелах.
    Интегрируется точно то, что физически накапливается вдоль перехода:

        A = exp( -2 k0 int Im n_eff(w(z)) dz ).
    """
    k0 = 2.0 * np.pi / lambda_um
    z = np.linspace(0.0, profile.length_um, steps + 1)
    widths = profile.width_at(z)
    im = np.array([abs(complex(neff_of_width(float(w))).imag) for w in widths])
    integral = float(trapz(im, z))
    power = float(np.exp(-2.0 * k0 * integral))
    return {
        "transmission": power,
        "loss_db": float(-10.0 * np.log10(max(power, 1e-300))),
        "mean_alpha_db_per_cm": float(
            np.mean(2.0 * k0 * im) * 1e4 * 10.0 / np.log(10.0)
        ),
    }


def love_adiabaticity(
    profile: TaperProfile,
    neff_of_width: Callable[[float], complex],
    n_competing: float,
    lambda_um: float,
    steps: int = 200,
) -> dict[str, float]:
    """Критерий адиабатичности Лава для локального угла сужения.

    Переход считается адиабатическим, если местный угол кромки меньше
    предельного:

        Omega_max(z) = rho(z) * (beta_1 - beta_2) / (2 pi)
                     = (w/2) * (n_eff - n_competing) / lambda,

    где rho - локальная полуширина, а beta_2 отвечает ближайшему конкурирующему
    решению. Для одномодового по ширине перехода конкурентом служит порог
    излучения, то есть показатель обкладки.

    Возвращает предельный угол в самой узкой точке, фактический угол и их
    отношение. Отношение больше единицы означает нарушение критерия.
    """
    z = np.linspace(0.0, profile.length_um, steps + 1)
    widths = profile.width_at(z)
    limits = []
    for w in widths:
        n = complex(neff_of_width(float(w))).real
        limits.append(0.5 * w * max(n - n_competing, 0.0) / lambda_um)
    limits_arr = np.array(limits)
    worst = float(np.min(limits_arr))
    actual = profile.half_angle_rad
    return {
        "omega_limit_rad": worst,
        "omega_limit_deg": float(np.degrees(worst)),
        "omega_actual_rad": actual,
        "omega_actual_deg": profile.half_angle_deg,
        "violation": float(actual / worst) if worst > 0 else float("inf"),
    }


def step_overlap_deficit(
    profile: TaperProfile,
    mode_profile: Callable[[float, np.ndarray], np.ndarray],
    x: np.ndarray,
    steps: int,
) -> float:
    """Диагностика сходимости: суммарный дефицит проекции на ступеньках.

    Не является физической оценкой излучения (см. предупреждение в шапке
    модуля): величина убывает при измельчении шага. Используется только чтобы
    показать, что ступенчатая аппроксимация сошлась.
    """
    z = np.linspace(0.0, profile.length_um, steps + 1)
    widths = profile.width_at(z)
    deficit = 0.0
    prev = mode_profile(float(widths[0]), x)
    prev = prev / np.sqrt(trapz(np.abs(prev) ** 2, x))
    for w in widths[1:]:
        cur = mode_profile(float(w), x)
        cur = cur / np.sqrt(trapz(np.abs(cur) ** 2, x))
        proj = trapz(prev * np.conj(cur), x)
        deficit += max(0.0, 1.0 - float(abs(proj) ** 2))
        prev = cur
    return deficit


def junction_loss_db(eta: float) -> float:
    """Потери резкого стыка по доле переданной мощности."""
    return float(-10.0 * np.log10(max(eta, 1e-300)))
