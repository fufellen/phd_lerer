"""
Стадия 2: RCWA (Rigorous Coupled-Wave Analysis / Fourier Modal Method) для
1D-периодической решётки (полоски вдоль y, период Λ по x), TE-поляризация
(E вдоль y). Не привязано к статьям Лерера - независимая, отдельно
валидируемая реализация (см. README.md про то, почему не порт его кода).

Матричный формализм - прямое обобщение уже провалидированного скалярного
tmm.py (interface_matrix/propagation_matrix), а не "компактная" gap-medium
формула из учебных заметок по RCWA (первая попытка с ней содержала ошибку,
не совпадала с tmm.py даже в тривиальном однослойном пределе - см. историю
коммитов). Вывод: на границе слоёв j/j+1 из непрерывности Ey и Hx
(TE-поляризация; W - матрица форм мод, V=W@Q - "H-подобная" матрица):

    W_j(c_j+ + c_j-) = W_{j+1}(c_{j+1}+ + c_{j+1}-)
    V_j(c_j+ - c_j-) = V_{j+1}(c_{j+1}+ - c_{j+1}-)

Отсюда матрица интерфейса (P=W_j^-1 W_{j+1}, Q'=V_j^-1 V_{j+1}):

    [c_j+; c_j-] = [[0.5(P+Q'), 0.5(P-Q')], [0.5(P-Q'), 0.5(P+Q')]] @ [c_{j+1}+; c_{j+1}-]

Проверено вручную: при W_j=W_{j+1}=1 (скалярный случай, TE) это в точности
совпадает с уже провалидированной формулой tmm.py (t=1+r, интерфейсная
матрица [[1/t,r/t],[r/t,1/t]]).

Упрощение геометрии: реальная решётка Лерера - двумерно периодическая
(цилиндры по x и y). Здесь - 1D-периодическая версия (полоски, период
только по x). Подробности и обоснование - в README.md.

Только TE (Ey, Hx, Hz) - не нужна факторизация 1/eps (правила Ли для TM),
меньше риск тонких ошибок. Для цели задачи (локализовать резонанс,
сравнить материалы) обеих поляризаций не требуется.
"""

from __future__ import annotations

import cmath
from dataclasses import dataclass

import numpy as np


def fourier_coeffs_strip(eps_incl: complex, eps_bg: complex, fill: float, N: int) -> np.ndarray:
    """Фурье-коэффициенты eps(x) для полоски шириной fill*Lambda, центрированной
    в x=0, период Lambda=1. Возвращает массив длины 2N+1 (индекс i = гармоника i-N).

    c_0 = f*eps_incl + (1-f)*eps_bg; c_i = (eps_incl-eps_bg)*sin(pi*i*f)/(pi*i), i!=0.
    """
    coeffs = np.zeros(2 * N + 1, dtype=complex)
    for idx, i in enumerate(range(-N, N + 1)):
        if i == 0:
            coeffs[idx] = fill * eps_incl + (1.0 - fill) * eps_bg
        else:
            coeffs[idx] = (eps_incl - eps_bg) * cmath.sin(cmath.pi * i * fill) / (cmath.pi * i)
    return coeffs


def toeplitz_from_coeffs(coeffs: np.ndarray) -> np.ndarray:
    """(2N+1)x(2N+1) матрица свёртки (Тёплица) из коэффициентов (индекс i = гармоника i-N)."""
    M = len(coeffs)
    N = (M - 1) // 2
    T = np.zeros((M, M), dtype=complex)
    for row in range(M):
        for col in range(M):
            k = row - col
            idx = k + N
            if 0 <= idx < M:
                T[row, col] = coeffs[idx]
    return T


@dataclass
class LayerModes:
    W: np.ndarray
    Q: np.ndarray
    V: np.ndarray


def layer_eigenmodes_te(eps_toeplitz: np.ndarray, Kx: np.ndarray) -> LayerModes:
    """Собственные моды TE: A = Kx@Kx - eps_toeplitz = W diag(-q^2) W^-1, ветвь Im(q)>=0."""
    A = Kx @ Kx - eps_toeplitz
    eigvals, W = np.linalg.eig(A)
    q = np.sqrt(-eigvals.astype(complex))
    for i in range(len(q)):
        if q[i].imag < 0:
            q[i] = -q[i]
        elif q[i].imag == 0 and q[i].real < 0:
            q[i] = -q[i]
    Q = np.diag(q)
    V = W @ Q
    return LayerModes(W=W, Q=Q, V=V)


def homogeneous_kz(eps: complex, Kx: np.ndarray) -> np.ndarray:
    """Прямой расчёт нормированных kz по гармоникам для однородного слоя (валидация)."""
    kx = np.diag(Kx)
    kz = np.sqrt(eps - kx ** 2 + 0j)
    kz = np.array([v if v.imag >= 0 else -v for v in kz])
    return kz


def _interface_matrix(W_j: np.ndarray, V_j: np.ndarray, W_j1: np.ndarray, V_j1: np.ndarray) -> np.ndarray:
    """[c_j+;c_j-] = IM @ [c_{j+1}+;c_{j+1}-] (см. докстринг модуля)."""
    M = W_j.shape[0]
    P = np.linalg.solve(W_j, W_j1)
    Qp = np.linalg.solve(V_j, V_j1)
    top = np.hstack([0.5 * (P + Qp), 0.5 * (P - Qp)])
    bot = np.hstack([0.5 * (P - Qp), 0.5 * (P + Qp)])
    return np.vstack([top, bot])


def _propagation_matrix(Q: np.ndarray, k0: float, d: float) -> np.ndarray:
    """[c(entry)+;c(entry)-] = PM @ [c(exit)+;c(exit)-] внутри одного слоя толщиной d."""
    M = Q.shape[0]
    phase = np.exp(1j * np.diag(Q) * k0 * d)
    top = np.hstack([np.diag(1.0 / phase), np.zeros((M, M), dtype=complex)])
    bot = np.hstack([np.zeros((M, M), dtype=complex), np.diag(phase)])
    return np.vstack([top, bot])


def solve_rcwa_te(
    layers: list[tuple[np.ndarray, float]],
    eps_in: complex,
    eps_out: complex,
    N: int,
    period_nm: float,
    lam_nm: float,
    theta_i_deg: float,
) -> tuple[np.ndarray, np.ndarray]:
    """RCWA TE для стопки periodic-слоёв (прямое обобщение tmm.solve_rt на матричный случай).

    layers: [(eps_toeplitz, thickness_nm), ...] сверху вниз. eps_in/eps_out -
    полубесконечные однородные среды по краям. Возвращает (R_harm, T_harm) -
    мощность по каждой гармонике (сумма по пропагирующим гармоникам физична:
    sum(R)+sum(T)=1 в беспотерьном случае).
    """
    Mh = 2 * N + 1
    k0 = 2.0 * cmath.pi / lam_nm
    theta_i = cmath.pi * theta_i_deg / 180.0
    kx0 = (cmath.sqrt(eps_in) * k0 * cmath.sin(theta_i)).real
    Gx = 2.0 * cmath.pi / period_nm
    harmonics = np.arange(-N, N + 1)
    kx_i = kx0 + harmonics * Gx
    Kx = np.diag(kx_i / k0)

    kz_in = homogeneous_kz(eps_in, Kx)
    kz_out = homogeneous_kz(eps_out, Kx)
    W_in, V_in = np.eye(Mh, dtype=complex), np.diag(kz_in)
    W_out, V_out = np.eye(Mh, dtype=complex), np.diag(kz_out)

    modes_list = [layer_eigenmodes_te(eps_t, Kx) for eps_t, _ in layers]

    # Полная цепочка матриц: [c_in+;c_in-] = M_total @ [c_out+; 0]
    M_total = np.eye(2 * Mh, dtype=complex)
    W_prev, V_prev = W_in, V_in
    for (eps_t, d), modes in zip(layers, modes_list):
        M_total = M_total @ _interface_matrix(W_prev, V_prev, modes.W, modes.V)
        M_total = M_total @ _propagation_matrix(modes.Q, k0, d)
        W_prev, V_prev = modes.W, modes.V
    M_total = M_total @ _interface_matrix(W_prev, V_prev, W_out, V_out)

    M11 = M_total[:Mh, :Mh]
    M21 = M_total[Mh:, :Mh]

    c_in_plus = np.zeros(Mh, dtype=complex)
    c_in_plus[N] = 1.0  # источник: только (0,0) гармоника, TE, амплитуда 1

    c_out_plus = np.linalg.solve(M11, c_in_plus)
    c_in_minus = M21 @ c_out_plus

    R_harm = (np.abs(c_in_minus) ** 2) * (kz_in.real / kz_in[N].real)
    T_harm = (np.abs(c_out_plus) ** 2) * (kz_out.real / kz_in[N].real)
    return R_harm, T_harm
