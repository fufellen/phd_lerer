"""
Стадия 4: RCWA/FMM для НАСТОЯЩЕЙ 2D-периодической решётки (crossed grating) -
цилиндры на квадратной решётке, как в статье Лерер АМ_1, а не 1D-полоски.

Формулировка - стандартная полновекторная FMM для би-периодических структур
(Moharam et al. 1995, расширение на crossed gratings; конвенции согласованы
с уже провалидированными tmm.py и rcwa_1d_te.py: временной фактор
exp(-i*omega*t), потери Im(eps)>0, ветвь Im(q)>=0, те же интерфейсные и
пропагационные матрицы, только блоки теперь векторные 2M x 2M, M=(2N+1)^2).

Вывод (нормировка: все волновые числа в единицах k0, h = Z0*H):

    из z-компонент уравнений Максвелла:
        hz = Kx ey - Ky ex
        ez = -E^-1 (Kx hy - Ky hx)
    подстановка в x,y-компоненты даёт систему 1-го порядка
        d e/dz' = i F h,   d h/dz' = i G e,   z' = k0 z
    с блочными матрицами (E = свёрточная (Тёплицева по обоим индексам)
    матрица eps(x,y); Kx, Ky - диагональные kx_mn/k0, ky_mn/k0):

        F = [[ Kx E^-1 Ky,      I - Kx E^-1 Kx ],
             [ Ky E^-1 Ky - I,  -Ky E^-1 Kx    ]]
        G = [[ -Kx Ky,          Kx^2 - E       ],
             [ E - Ky^2,        Ky Kx          ]]

    Собственные моды: (F G) W = W diag(q^2), ветвь Im(q)>=0;
    "магнитная" матрица V = G W diag(1/q). Проверка руками: однородная среда,
    нормальное падение, e=(0, ey): hx = -eps*ey/q = -n*ey - совпадает с
    плоской волной E=y*Ey, H=-x*n*Ey/Z0. Знаки согласованы с tmm.py.

Факторизация eps - прямое правило Лорана (Laurent rule). Для crossed
gratings правила Ли (инверсное правило вдоль нормали к разрыву) сходятся
быстрее, но правило Лорана корректно сходится с ростом N; для нашего
умеренного контраста (композит |eps|~3 против воздуха) сходимость
проверяется явно в validate_rcwa_2d.py. Это осознанный компромисс
"проще = меньше риск тонкой ошибки в правилах факторизации".

Мощности: поток Пойнтинга по гармоникам считается из ПОЛНЫХ поперечных
полей гармоники Sz_mn ~ Re(ex hy* - ey hx*), а не из |c|^2 по модам -
в 2D базисные моды одной гармоники (x-pol и y-pol) дают перекрёстные
члены в потоке на косых порядках (Kx*Ky != 0), и |c|^2 там неверен.
"""

from __future__ import annotations

import cmath
import math
from dataclasses import dataclass

import numpy as np


# ---------------------------------------------------------------------------
# Бессель J1 (scipy на этой машине нет) - рациональные аппроксимации
# Abramowitz & Stegun 9.4.4/9.4.6, |err| < 4e-8 - достаточно для
# Фурье-коэффициентов диска.
# ---------------------------------------------------------------------------

def bessel_j1(x: float) -> float:
    ax = abs(x)
    if ax < 8.0:
        y = x * x
        p1 = x * (72362614232.0 + y * (-7895059235.0 + y * (242396853.1
             + y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606))))))
        p2 = 144725228442.0 + y * (2300535178.0 + y * (18583304.74
             + y * (99447.43394 + y * (376.9991397 + y))))
        return p1 / p2
    z = 8.0 / ax
    y = z * z
    xx = ax - 2.356194491
    p1 = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4
         + y * (0.2457520174e-5 + y * (-0.240337019e-6))))
    p2 = 0.04687499995 + y * (-0.2002690873e-3 + y * (0.8449199096e-5
         + y * (-0.88228987e-6 + y * 0.105787412e-6)))
    ans = math.sqrt(0.636619772 / ax) * (math.cos(xx) * p1 - z * math.sin(xx) * p2)
    return ans if x >= 0.0 else -ans


# ---------------------------------------------------------------------------
# Фурье-коэффициенты 2D-профилей eps(x,y), период Lambda по обеим осям
# ---------------------------------------------------------------------------

def fourier_coeffs_cylinder(eps_incl: complex, eps_bg: complex,
                            diameter_over_period: float, N: int) -> np.ndarray:
    """Коэффициенты a[m+N, n+N] для eps(x,y) = eps_bg + (eps_incl-eps_bg)*disk,
    диск радиуса R = (d/Λ)/2 в центре квадратной ячейки с периодом 1.

    Аналитика: интеграл индикатора диска = f_area * 2*J1(g R)/(g R),
    g = 2*pi*sqrt(m^2+n^2), f_area = pi R^2 (доля площади).
    """
    R = 0.5 * diameter_over_period
    f_area = math.pi * R * R
    size = 2 * N + 1
    a = np.zeros((size, size), dtype=complex)
    delta = eps_incl - eps_bg
    for m in range(-N, N + 1):
        for n in range(-N, N + 1):
            if m == 0 and n == 0:
                a[N, N] = eps_bg + delta * f_area
            else:
                g = 2.0 * math.pi * math.hypot(m, n)
                a[m + N, n + N] = delta * f_area * 2.0 * bessel_j1(g * R) / (g * R)
    return a


def fourier_coeffs_stripe_2d(eps_incl: complex, eps_bg: complex,
                             fill: float, N: int) -> np.ndarray:
    """1D-полоска (вдоль y) в 2D-представлении: a[m,n] != 0 только при n=0.
    Нужна для валидационного теста 'полный 2D-код против rcwa_1d_te'."""
    size = 2 * N + 1
    a = np.zeros((size, size), dtype=complex)
    for m in range(-N, N + 1):
        if m == 0:
            a[N, N] = fill * eps_incl + (1.0 - fill) * eps_bg
        else:
            a[m + N, N] = (eps_incl - eps_bg) * math.sin(math.pi * m * fill) / (math.pi * m)
    return a


def convolution_matrix_2d(a: np.ndarray, N: int) -> np.ndarray:
    """Свёрточная матрица E[(m,n),(m',n')] = a[m-m', n-n'], M=(2N+1)^2.

    Линейный индекс гармоники (m,n): idx = (m+N)*(2N+1) + (n+N).
    """
    size = 2 * N + 1
    M = size * size
    E = np.zeros((M, M), dtype=complex)
    for m in range(-N, N + 1):
        for n in range(-N, N + 1):
            row = (m + N) * size + (n + N)
            for mp in range(-N, N + 1):
                dm = m - mp
                if abs(dm) > N:
                    continue
                for np_ in range(-N, N + 1):
                    dn = n - np_
                    if abs(dn) > N:
                        continue
                    col = (mp + N) * size + (np_ + N)
                    E[row, col] = a[dm + N, dn + N]
    return E


# ---------------------------------------------------------------------------
# Слои: собственные моды
# ---------------------------------------------------------------------------

@dataclass
class LayerModes2D:
    W: np.ndarray  # 2M x 2M, формы мод (e = [ex; ey])
    q: np.ndarray  # 2M, нормированные kz/k0 мод
    V: np.ndarray  # 2M x 2M, h = [hx; hy] мод


def _fg_matrices(E: np.ndarray, Kx: np.ndarray, Ky: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    M = E.shape[0]
    I = np.eye(M, dtype=complex)
    Einv_Kx = np.linalg.solve(E, Kx)
    Einv_Ky = np.linalg.solve(E, Ky)
    F = np.block([
        [Kx @ Einv_Ky, I - Kx @ Einv_Kx],
        [Ky @ Einv_Ky - I, -Ky @ Einv_Kx],
    ])
    G = np.block([
        [-Kx @ Ky, Kx @ Kx - E],
        [E - Ky @ Ky, Ky @ Kx],
    ])
    return F, G


def _branch(q: np.ndarray) -> np.ndarray:
    q = q.copy()
    for i in range(len(q)):
        if q[i].imag < 0:
            q[i] = -q[i]
        elif q[i].imag == 0 and q[i].real < 0:
            q[i] = -q[i]
    return q


def layer_eigenmodes_2d(E: np.ndarray, Kx: np.ndarray, Ky: np.ndarray) -> LayerModes2D:
    """Моды периодического слоя: (F G) W = W diag(q^2), V = G W diag(1/q)."""
    F, G = _fg_matrices(E, Kx, Ky)
    eigvals, W = np.linalg.eig(F @ G)
    q = _branch(np.sqrt(eigvals.astype(complex)))
    V = G @ W @ np.diag(1.0 / q)
    return LayerModes2D(W=W, q=q, V=V)


def homogeneous_modes_2d(eps: complex, Kx: np.ndarray, Ky: np.ndarray) -> LayerModes2D:
    """Аналитический базис однородной среды: W = I, q_mn = sqrt(eps-kx^2-ky^2)
    (дважды вырожден по поляризации), V = G @ diag(1/q). Нужен для
    полубесконечных сред: у них должен быть детерминированный базис
    (столбец = чистая ex- или ey-гармоника), чтобы задавать падающую волну
    и читать амплитуды порядков; numpy.linalg.eig при кратных собственных
    значениях вернул бы произвольную смесь поляризаций."""
    M = Kx.shape[0]
    kx = np.diag(Kx)
    ky = np.diag(Ky)
    q1 = np.sqrt(eps - kx ** 2 - ky ** 2 + 0j)
    q1 = _branch(q1)
    q = np.concatenate([q1, q1])  # (ex-блок, ey-блок) - одинаковые q
    I = np.eye(M, dtype=complex)
    E = eps * I
    _, G = _fg_matrices(E, Kx, Ky)
    W = np.eye(2 * M, dtype=complex)
    V = G @ np.diag(1.0 / q)
    return LayerModes2D(W=W, q=q, V=V)


# ---------------------------------------------------------------------------
# Сборка стопки - та же схема интерфейс/пропагация, что в rcwa_1d_te.py
# ---------------------------------------------------------------------------

def _interface_matrix(W_j, V_j, W_j1, V_j1) -> np.ndarray:
    P = np.linalg.solve(W_j, W_j1)
    Qp = np.linalg.solve(V_j, V_j1)
    top = np.hstack([0.5 * (P + Qp), 0.5 * (P - Qp)])
    bot = np.hstack([0.5 * (P - Qp), 0.5 * (P + Qp)])
    return np.vstack([top, bot])


def _propagation_matrix(q: np.ndarray, k0: float, d: float) -> np.ndarray:
    M2 = len(q)
    phase = np.exp(1j * q * k0 * d)
    top = np.hstack([np.diag(1.0 / phase), np.zeros((M2, M2), dtype=complex)])
    bot = np.hstack([np.zeros((M2, M2), dtype=complex), np.diag(phase)])
    return np.vstack([top, bot])


def solve_rcwa_2d(
    layers: list[tuple[np.ndarray, float]],
    eps_in: complex,
    eps_out: complex,
    N: int,
    period_nm: float,
    lam_nm: float,
    pol: str = "y",
) -> tuple[np.ndarray, np.ndarray]:
    """RCWA для стопки 2D-периодических слоёв, нормальное падение.

    layers: [(E_conv | eps_scalar, thickness_nm), ...] сверху вниз. E_conv -
    свёрточная матрица (M x M, M=(2N+1)^2) из convolution_matrix_2d; если
    вместо матрицы передан комплексный скаляр eps, слой однородный и для
    него берётся аналитический базис (быстрее и без вырожденного eig).
    eps_in/eps_out - полубесконечные однородные среды. pol: 'x' или 'y' - поляризация
    падающей волны (для цилиндров на квадратной решётке при нормальном
    падении обе дают одинаковый результат - это валидационный тест).

    Возвращает (R_mn, T_mn) - мощность по гармоникам, массивы (2N+1, 2N+1),
    индекс [m+N, n+N].
    """
    size = 2 * N + 1
    M = size * size
    k0 = 2.0 * cmath.pi / lam_nm
    Gx = 2.0 * cmath.pi / period_nm
    harmonics = np.arange(-N, N + 1)
    # kx зависит только от m, ky только от n (нормальное падение: kx0=ky0=0)
    kx = np.repeat(harmonics, size) * Gx / k0
    ky = np.tile(harmonics, size) * Gx / k0
    Kx = np.diag(kx.astype(complex))
    Ky = np.diag(ky.astype(complex))

    modes_in = homogeneous_modes_2d(eps_in, Kx, Ky)
    modes_out = homogeneous_modes_2d(eps_out, Kx, Ky)

    M_total = np.eye(4 * M, dtype=complex)
    W_prev, V_prev = modes_in.W, modes_in.V
    for E_conv, d in layers:
        if np.isscalar(E_conv):
            modes = homogeneous_modes_2d(complex(E_conv), Kx, Ky)
        else:
            modes = layer_eigenmodes_2d(E_conv, Kx, Ky)
        M_total = M_total @ _interface_matrix(W_prev, V_prev, modes.W, modes.V)
        M_total = M_total @ _propagation_matrix(modes.q, k0, d)
        W_prev, V_prev = modes.W, modes.V
    M_total = M_total @ _interface_matrix(W_prev, V_prev, modes_out.W, modes_out.V)

    M11 = M_total[:2 * M, :2 * M]
    M21 = M_total[2 * M:, :2 * M]

    center = N * size + N  # линейный индекс гармоники (0,0)
    c_in_plus = np.zeros(2 * M, dtype=complex)
    if pol == "x":
        c_in_plus[center] = 1.0            # ex-блок
    elif pol == "y":
        c_in_plus[M + center] = 1.0        # ey-блок
    else:
        raise ValueError("pol must be 'x' or 'y'")

    c_out_plus = np.linalg.solve(M11, c_in_plus)
    c_in_minus = M21 @ c_out_plus

    # Потоки Пойнтинга по гармоникам из полных поперечных полей.
    # Вперёд (+z): e = W c+, h = +V c+. Назад (-z): e = W c-, h = -V c-.
    def flux_per_harmonic(modes: LayerModes2D, c: np.ndarray, direction: float) -> np.ndarray:
        e = modes.W @ c
        h = direction * (modes.V @ c)
        ex, ey = e[:M], e[M:]
        hx, hy = h[:M], h[M:]
        sz = np.real(ex * np.conj(hy) - ey * np.conj(hx))
        return sz

    sz_inc = flux_per_harmonic(modes_in, c_in_plus, +1.0)
    inc_total = sz_inc.sum()  # вся падающая мощность в (0,0)
    sz_ref = flux_per_harmonic(modes_in, c_in_minus, -1.0)
    sz_trn = flux_per_harmonic(modes_out, c_out_plus, +1.0)

    R_lin = -sz_ref / inc_total  # поток отражённой идёт в -z
    T_lin = sz_trn / inc_total
    return R_lin.reshape(size, size), T_lin.reshape(size, size)
