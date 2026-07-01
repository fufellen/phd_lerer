"""
Стадия 1 солвера RCWA/Fourier Modal Method для 2D-периодической решётки
Лерера (см. Obsidian CODEX-чекпоинт: "Контекст задачи Codex - Лерер
метаповерхности 2D solver.md" рядом со статьёй).

Эта стадия - НЕ ещё решётка. Это transfer-matrix method (TMM/Parratt) для
плоской многослойной среды: расчёт R, T для стопки однородных слоёв под
произвольным углом падения, обе поляризации. Причина начинать отсюда: в
RCWA каждый слой раскладывается в ряды Фурье по решётке, но нулевая
Фурье-гармоника периодической решётки - это ровно однородный слой с
эффективной eps, то есть TMM - предельный случай (число гармоник N=0) и
одновременно обязательный строительный блок для сшивки мод между слоями
(S-матрица) в полной RCWA. Здесь TMM реализован и провалидирован независимо
от статей Лерера (учебная физика, формулы Френеля/Эйри), поэтому его
корректность не зависит от того, что могло потеряться при разборе PDF.

Что дальше (см. README.md и чекпоинт): поверх этого добавляется
Фурье-факторизация eps(x,y) периодической решётки, послойные собственные
моды (решение задачи на собственные значения для гармоник Флоке) и сшивка
через S-матрицу - только тогда получится настоящий расчёт двухпериодической
решётки, а не однородного слоя.
"""

from __future__ import annotations

import cmath
from dataclasses import dataclass


def _kz(n: complex, k0: float, kx: float) -> complex:
    """z-компонента волнового вектора в среде с показателем преломления n.

    Ветвь выбрана так, чтобы Im(kz) >= 0 (волна затухает/уходит в сторону
    +z вглубь поглощающей полубесконечной среды, стандартное условие
    излучения для комплексных n).
    """
    val = cmath.sqrt((n * k0) ** 2 - kx ** 2)
    if val.imag < 0:
        val = -val
    return val


@dataclass
class Layer:
    n: complex  # комплексный показатель преломления слоя
    d: float  # толщина, нм (для двух крайних (подложка/выходная среда) не используется)


def solve_rt(layers: list[Layer], lam_nm: float, theta_i_deg: float, pol: str) -> tuple[complex, complex]:
    """TMM/Parratt recursion для амплитудных r, t многослойной стопки.

    layers[0] - среда падения (полубесконечная), layers[-1] - среда выхода
    (полубесконечная), между ними - конечные слои с толщинами d в нм.
    pol: 's' (TE) или 'p' (TM).
    Возвращает (r, t) - амплитудные коэффициенты отражения/прохождения по
    полю (для мощности использовать power_rt()).
    """
    if pol not in ("s", "p"):
        raise ValueError("pol must be 's' or 'p'")
    k0 = 2.0 * cmath.pi / lam_nm
    theta_i = cmath.pi * theta_i_deg / 180.0
    n0 = layers[0].n
    kx = (n0 * k0 * cmath.sin(theta_i)).real  # kx вещественна и сохраняется по всем слоям (Снеллиус)

    kzs = [_kz(layer.n, k0, kx) for layer in layers]

    def fresnel_r(j: int) -> complex:
        n_j, n_j1 = layers[j].n, layers[j + 1].n
        kz_j, kz_j1 = kzs[j], kzs[j + 1]
        if pol == "s":
            return (kz_j - kz_j1) / (kz_j + kz_j1)
        # TM (p): использует kz/eps = kz/n^2
        a_j = kz_j / (n_j ** 2)
        a_j1 = kz_j1 / (n_j1 ** 2)
        return (a_j - a_j1) / (a_j + a_j1)

    def fresnel_t(j: int) -> complex:
        # t = 1 + r в обоих случаях: для TE r определён через "импеданс"
        # Z=kz, для TM - через Z=kz/n^2 (см. fresnel_r); t=1+r - это условие
        # непрерывности поля ψ (Ey для TE, Hy для TM) на границе, оно не
        # зависит от того, что представляет собой Z, поэтому одна и та же
        # формула верна для обеих поляризаций при условии, что power_rt()
        # ниже считает поток мощности в терминах ТОГО ЖЕ ψ (см. там).
        return 1.0 + fresnel_r(j)

    N = len(layers)  # layers[0..N-1], N-2 конечных слоёв (индексы 1..N-2)

    # Рекурсия Парратта: начинаем с последней границы (N-2 -> N-1) и идём назад.
    R = fresnel_r(N - 2)
    for j in range(N - 3, -1, -1):
        d_next = layers[j + 1].d
        phase = cmath.exp(2j * kzs[j + 1] * d_next)
        r_j = fresnel_r(j)
        R = (r_j + R * phase) / (1.0 + r_j * R * phase)

    # Амплитудный t получаем прямым проходом (T-матрица), последовательно
    # накапливая произведение слоевых t_j и коэффициент фазового набега.
    T = complex(1.0)
    R_forward = R  # R на входе в стопку (перед layers[0]/layers[1])
    # Для получения t выполняем прямой пересчёт с использованием матриц
    # переноса (интерфейс + пропагация), это надёжнее, чем восстанавливать
    # t только из рекурсии для R.
    M = _identity_2x2()
    for j in range(N - 1):
        r_j = fresnel_r(j)
        t_j = fresnel_t(j)
        M = _mat_mul(M, _interface_matrix(r_j, t_j))
        if j + 1 <= N - 2:  # пропагация внутри конечного слоя j+1 (не для последней полубесконечной среды)
            M = _mat_mul(M, _propagation_matrix(kzs[j + 1], layers[j + 1].d))

    # M связывает (E0+, E0-) со (EN+, EN-) как [E0+, E0-]^T = M [EN+, 0]^T (нет волны, идущей навстречу из среды выхода)
    r = M[1][0] / M[0][0]
    t = 1.0 / M[0][0]
    return r, t


def _identity_2x2():
    return [[1.0, 0.0], [0.0, 1.0]]


def _mat_mul(a, b):
    return [
        [a[0][0] * b[0][0] + a[0][1] * b[1][0], a[0][0] * b[0][1] + a[0][1] * b[1][1]],
        [a[1][0] * b[0][0] + a[1][1] * b[1][0], a[1][0] * b[0][1] + a[1][1] * b[1][1]],
    ]


def _interface_matrix(r: complex, t: complex):
    # Стандартная матрица переноса границы (Transfer matrix), связывающая
    # амплитуды по обе стороны интерфейса через r, t этого интерфейса.
    return [[1.0 / t, r / t], [r / t, 1.0 / t]]


def _propagation_matrix(kz: complex, d: float):
    phase = cmath.exp(1j * kz * d)
    return [[1.0 / phase, 0.0], [0.0, phase]]


def power_rt(layers: list[Layer], lam_nm: float, theta_i_deg: float, pol: str) -> tuple[float, float]:
    """R, T по мощности (энергии), с учётом разных сред на входе/выходе."""
    r, t = solve_rt(layers, lam_nm, theta_i_deg, pol)
    k0 = 2.0 * cmath.pi / lam_nm
    theta_i = cmath.pi * theta_i_deg / 180.0
    n0, nN = layers[0].n, layers[-1].n
    kx = (n0 * k0 * cmath.sin(theta_i)).real
    kz0 = _kz(n0, k0, kx)
    kzN = _kz(nN, k0, kx)
    R = abs(r) ** 2
    if pol == "s":
        T = (abs(t) ** 2) * (kzN.real / kz0.real)
    else:
        T = (abs(t) ** 2) * ((kzN / nN ** 2).real / (kz0 / n0 ** 2).real)
    return R, T
