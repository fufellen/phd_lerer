"""Ввод излучения: перекрытие с модой волокна и распространение в разрыве.

Здесь собраны три вещи, нужные для воспроизведения торцевого ввода и структур
с разрывом волновода:

1. двумерный профиль моды полоски в приближении метода эффективного показателя -
   произведение точного вертикального профиля на горизонтальный;
2. перекрытие двух двумерных полей (скалярное приближение, стандартное для
   расчёта стыковки волокна с волноводом);
3. распространение поля через однородный участок методом углового спектра -
   для разрыва в волноводе, где мода свободно дифрагирует.

Скалярное приближение оправдано тем, что поля волокна и слабо локализованной
LR-моды почти поперечны и почти сонаправлены; для сильнолокализованных мод
нужна векторная форма из tmm.overlap_power.
"""

from __future__ import annotations

import numpy as np

from .tmm import Stack, mode_fields, trapz


def gaussian_field(x: np.ndarray, z: np.ndarray, mfd_x_um: float,
                   mfd_z_um: float | None = None) -> np.ndarray:
    """Поле одномодового волокна как гауссов пучок с заданным диаметром пятна.

    Диаметр модового пятна MFD задаётся по уровню 1/e^2 по интенсивности, как
    это принято в спецификациях волокна, поэтому радиус поля равен MFD/2.
    """
    mfd_z_um = mfd_x_um if mfd_z_um is None else mfd_z_um
    wx, wz = mfd_x_um / 2.0, mfd_z_um / 2.0
    xx, zz = np.meshgrid(x, z, indexing="ij")
    return np.exp(-((xx / wx) ** 2) - ((zz / wz) ** 2))


def strip_mode_field(
    neff_planar: complex,
    neff_strip: complex,
    stack: Stack,
    k0: float,
    width_um: float,
    eps_side: complex,
    x: np.ndarray,
    z: np.ndarray,
) -> np.ndarray:
    """Двумерный профиль моды полоски: вертикальный точный на горизонтальный ЭДП."""
    centre = 0.5 * stack.interfaces()[-1]
    vertical = mode_fields(neff_planar, stack, k0, z + centre)[0]

    half = width_um / 2.0
    u = k0 * np.sqrt(neff_planar**2 - neff_strip * neff_strip + 0j)
    alpha = k0 * np.sqrt(neff_strip * neff_strip - eps_side + 0j)
    if alpha.real < 0:
        alpha = -alpha
    horizontal = np.zeros_like(x, dtype=complex)
    inside = np.abs(x) <= half
    horizontal[inside] = np.cos(u * x[inside])
    horizontal[~inside] = np.cos(u * half) * np.exp(-alpha * (np.abs(x[~inside]) - half))

    return np.outer(horizontal, vertical)


def overlap_2d(a: np.ndarray, b: np.ndarray, x: np.ndarray, z: np.ndarray) -> float:
    """Нормированное перекрытие двух двумерных полей по мощности."""
    num = abs(trapz(trapz(a * np.conj(b), z), x)) ** 2
    da = float(np.real(trapz(trapz(np.abs(a) ** 2, z), x)))
    db = float(np.real(trapz(trapz(np.abs(b) ** 2, z), x)))
    if da <= 0 or db <= 0:
        return float("nan")
    return float(num / (da * db))


def coupling_loss_db(eta: float) -> float:
    return float(-10.0 * np.log10(max(eta, 1e-300)))


def propagate_angular_spectrum(
    field: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    distance_um: float,
    n_medium: float,
    lambda_um: float,
) -> np.ndarray:
    """Распространение поля на заданное расстояние в однородной среде.

    Метод углового спектра: поле раскладывается по плоским волнам, каждая
    набирает свою фазу, затем поле собирается обратно.

        E(x, z; L) = F^-1 { F{E(x, z; 0)} * exp(i k_par L) },
        k_par = sqrt( (2 pi n / lambda)^2 - k_x^2 - k_z^2 ).

    Затухающие составляющие спектра (подкоренное выражение отрицательно)
    экспоненциально подавляются, что и даёт дифракционное расплывание пучка.
    """
    k = 2.0 * np.pi * n_medium / lambda_um
    dx = float(x[1] - x[0])
    dz = float(z[1] - z[0])
    kx = 2.0 * np.pi * np.fft.fftfreq(len(x), d=dx)
    kz = 2.0 * np.pi * np.fft.fftfreq(len(z), d=dz)
    kxx, kzz = np.meshgrid(kx, kz, indexing="ij")

    kpar2 = k * k - kxx**2 - kzz**2
    kpar = np.sqrt(kpar2.astype(complex))
    kpar = np.where(kpar.imag < 0, np.conj(kpar), kpar)  # затухание, а не рост

    spectrum = np.fft.fft2(field)
    return np.fft.ifft2(spectrum * np.exp(1j * kpar * distance_um))


def gap_transmission(
    field: np.ndarray,
    x: np.ndarray,
    z: np.ndarray,
    gap_um: float,
    n_medium: float,
    lambda_um: float,
    output_field: np.ndarray | None = None,
) -> float:
    """Доля мощности, захваченная модой после свободного участка длиной gap.

    Поле моды на входе в разрыв распространяется методом углового спектра, затем
    проектируется на моду выходного волновода.
    """
    propagated = propagate_angular_spectrum(field, x, z, gap_um, n_medium, lambda_um)
    target = field if output_field is None else output_field
    return overlap_2d(propagated, target, x, z)
