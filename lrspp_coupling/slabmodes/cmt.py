"""Связанные волноводы: от супермод к перекачке мощности.

Два подхода к одному эффекту.

1. Модовый. Для пары связанных волноводов решается общая слоистая задача и
   находятся две супермоды - чётная и нечётная. Их разность постоянных
   распространения задаёт длину связи

       L_c = pi / |beta_e - beta_o| = lambda / (2 |Re n_e - Re n_o|).

2. Через уравнения связанных волн. Коэффициент связи kappa и расстройка delta
   выражаются через те же супермоды:

       kappa = k0 * |Re n_e - Re n_o| / 2,   delta = k0 * (Re n_1 - Re n_2) / 2,

   где n_1 и n_2 - показатели изолированных волноводов. Предельная доля
   перекачанной мощности равна kappa^2 / (kappa^2 + delta^2): при расстройке
   перекачка неполна даже на оптимальной длине.

Потери учитываются отдельно: супермоды затухают со своими Im(beta), и в первом
приближении огибающая перекачки спадает по среднему затуханию пары.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def extract_kappa(
    lambda_um: float,
    neff_even: complex,
    neff_odd: complex,
    neff_guide_a: complex,
    neff_guide_b: complex,
) -> float:
    """Извлекает коэффициент связи из точного расщепления пары супермод.

    Связь и расстройка входят в расщепление совместно:

        (beta_e - beta_o) / 2 = sqrt(kappa^2 + delta^2),
        delta = (beta_a - beta_b) / 2,

    поэтому по известным супермодам и известным изолированным волноводам
    коэффициент связи восстанавливается однозначно:

        kappa = sqrt( ((beta_e - beta_o)/2)^2 - delta^2 ).

    Это позволяет посчитать kappa там, где задача решается точно (планарная
    вертикальная структура), и затем применить её к реальному прибору с другой
    расстройкой - например, после учёта конечной ширины полосок. Связь задаётся
    перекрытием затухающих хвостов поперёк зазора и от ширины полосок почти не
    зависит, тогда как расстройка зависит от ширины сильно.
    """
    k0 = 2.0 * np.pi / lambda_um
    half_split = k0 * abs(neff_even.real - neff_odd.real) / 2.0
    delta = k0 * (neff_guide_a.real - neff_guide_b.real) / 2.0
    value = half_split * half_split - delta * delta
    return float(np.sqrt(value)) if value > 0 else 0.0


@dataclass(frozen=True)
class CoupledPair:
    """Пара волноводов, заданная коэффициентом связи и расстройкой.

    В отличие от CouplerModel, который берёт готовые супермоды, здесь связь и
    расстройка заданы раздельно: связь получена из точного планарного расчёта, а
    расстройка - из эффективных показателей волноводов конечной ширины.
    """

    lambda_um: float
    kappa_per_um: float
    neff_a: complex
    neff_b: complex

    @property
    def k0(self) -> float:
        return 2.0 * np.pi / self.lambda_um

    @property
    def detuning_per_um(self) -> float:
        return float(self.k0 * (self.neff_a.real - self.neff_b.real) / 2.0)

    @property
    def half_split(self) -> float:
        return float(np.hypot(self.kappa_per_um, self.detuning_per_um))

    @property
    def coupling_length_um(self) -> float:
        """Длина максимальной перекачки без учёта потерь: pi / (2 sqrt(k^2+d^2))."""
        s = self.half_split
        return float(np.pi / (2.0 * s)) if s > 0 else float("inf")

    @property
    def max_transfer(self) -> float:
        k, d = self.kappa_per_um, self.detuning_per_um
        denom = k * k + d * d
        return float(k * k / denom) if denom > 0 else 0.0

    @property
    def mean_loss_db_per_cm(self) -> float:
        im = 0.5 * (abs(self.neff_a.imag) + abs(self.neff_b.imag))
        return float(2.0 * self.k0 * im * 1e4 * 10.0 / np.log(10.0))

    def transfer_curve(self, length_um: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Перекачка мощности с учётом среднего затухания пары."""
        z = np.asarray(length_um, dtype=float)
        s = self.half_split
        alpha = self.k0 * 0.5 * (abs(self.neff_a.imag) + abs(self.neff_b.imag))
        envelope = np.exp(-2.0 * alpha * z)
        cross = self.max_transfer * np.sin(s * z) ** 2 * envelope
        thru = (1.0 - self.max_transfer * np.sin(s * z) ** 2) * envelope
        return cross, thru

    def optimal_length_um(self, z_max_um: float = 6000.0, points: int = 30001) -> tuple[float, float]:
        z = np.linspace(1.0, z_max_um, points)
        cross, _ = self.transfer_curve(z)
        idx = int(np.argmax(cross))
        return float(z[idx]), float(cross[idx])


@dataclass(frozen=True)
class CouplerModel:
    """Параметры связанной пары, восстановленные из супермод."""

    lambda_um: float
    neff_even: complex
    neff_odd: complex
    neff_guide_a: complex | None = None
    neff_guide_b: complex | None = None

    @property
    def k0(self) -> float:
        return 2.0 * np.pi / self.lambda_um

    @property
    def delta_n(self) -> float:
        return abs(self.neff_even.real - self.neff_odd.real)

    @property
    def coupling_length_um(self) -> float:
        """Длина полной перекачки в отсутствие расстройки и потерь."""
        if self.delta_n == 0:
            return float("inf")
        return float(self.lambda_um / (2.0 * self.delta_n))

    @property
    def kappa_per_um(self) -> float:
        return float(self.k0 * self.delta_n / 2.0)

    @property
    def detuning_per_um(self) -> float:
        """Расстройка по изолированным волноводам; 0, если они не заданы."""
        if self.neff_guide_a is None or self.neff_guide_b is None:
            return 0.0
        return float(self.k0 * (self.neff_guide_a.real - self.neff_guide_b.real) / 2.0)

    @property
    def max_transfer(self) -> float:
        """Предельная доля перекачанной мощности без учёта потерь."""
        kap = self.kappa_per_um
        det = self.detuning_per_um
        if kap == 0 and det == 0:
            return 0.0
        return float(kap * kap / (kap * kap + det * det))

    @property
    def mean_loss_db_per_cm(self) -> float:
        """Среднее затухание пары супермод, дБ/см."""
        im = 0.5 * (abs(self.neff_even.imag) + abs(self.neff_odd.imag))
        alpha_per_um = 2.0 * self.k0 * im
        return float(alpha_per_um * 1e4 * 10.0 / np.log(10.0))

    def transfer_curve(self, length_um: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Доли мощности в перекачанном и исходном волноводах.

        Интерференция двух затухающих супермод с равными начальными весами:
        поле в каждом волноводе есть их сумма и разность, поэтому

            P_cross = |exp(i b_e z) - exp(i b_o z)|^2 / 4,
            P_thru  = |exp(i b_e z) + exp(i b_o z)|^2 / 4,

        где b = k0 * n_eff комплексные. Расстройка учитывается множителем
        max_transfer, поскольку при delta != 0 полная перекачка недостижима.
        """
        z = np.asarray(length_um, dtype=float)
        be = self.k0 * self.neff_even
        bo = self.k0 * self.neff_odd
        ee = np.exp(1j * be.real * z) * np.exp(-self.k0 * abs(self.neff_even.imag) * z)
        eo = np.exp(1j * bo.real * z) * np.exp(-self.k0 * abs(self.neff_odd.imag) * z)
        cross = np.abs(ee - eo) ** 2 / 4.0
        thru = np.abs(ee + eo) ** 2 / 4.0
        eta = self.max_transfer
        return cross * eta, thru + cross * (1.0 - eta)


def optimal_length_um(model: CouplerModel, z_max_um: float = 4000.0, points: int = 20001) -> tuple[float, float]:
    """Длина максимальной перекачки с учётом потерь и соответствующая доля.

    Из-за поглощения оптимум лежит раньше формальной длины связи.
    """
    z = np.linspace(1.0, z_max_um, points)
    cross, _ = model.transfer_curve(z)
    idx = int(np.argmax(cross))
    return float(z[idx]), float(cross[idx])
