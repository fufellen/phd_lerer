"""Второй шаг метода эффективного показателя: переход от плёнки к полоске.

Первый шаг - вертикальная планарная задача (slabmodes.tmm), она даёт комплексный
n_eff бесконечно широкой структуры. Второй шаг решает горизонтальную задачу для
полоски конечной ширины: сердцевина с eps = n_eff^2, обкладки с показателем
области вне полоски.

Ограничение метода. Это редуцированная модель: она не описывает векторную
гибридизацию на углах полоски и вблизи отсечки завышает локализацию. Границы
применимости для плазмонных волноводов разобраны в заметке ваута
"Итоговая проверка применимости ЭДП Лерера к плазмонным нановолноводам".
"""

from __future__ import annotations

import numpy as np

from .tmm import decaying_sqrt, solve_mode

try:
    from scipy.optimize import root as _scipy_root
except ModuleNotFoundError:  # pragma: no cover
    _scipy_root = None


def horizontal_residual_tm(
    neff: complex, k0: float, eps_core: complex, eps_side: complex, width_um: float
) -> complex:
    """Симметричная горизонтальная задача, чётное TM-подобное решение.

    Для LR-моды доминирующая компонента электрического поля лежит вдоль
    горизонтальной стратификации, поэтому берётся форма с весом 1/eps, как в
    ЭДП Лерера (см. edp/metal_strip_w800).
    """
    half = width_um / 2.0
    u = k0 * np.sqrt(eps_core - neff * neff + 0j)
    alpha = k0 * decaying_sqrt(neff * neff - eps_side)
    return (u / eps_core) * np.tan(u * half) - alpha / eps_side


def horizontal_residual_te(
    neff: complex, k0: float, eps_core: complex, eps_side: complex, width_um: float
) -> complex:
    """То же, но без весов 1/eps - чётное TE-подобное решение."""
    half = width_um / 2.0
    u = k0 * np.sqrt(eps_core - neff * neff + 0j)
    alpha = k0 * decaying_sqrt(neff * neff - eps_side)
    return u * np.tan(u * half) - alpha


def solve_strip(
    k0: float,
    eps_core: complex,
    eps_side: complex,
    width_um: float,
    polarization: str = "tm",
    extra_guesses: tuple[complex, ...] = (),
) -> complex:
    """Находит основную моду полоски конечной ширины.

    Возвращает комплексный n_eff. Корень ищется из набора приближений и
    отбирается по двум условиям: он лежит между показателем боковой области и
    планарным n_eff, и он максимален по действительной части (основная мода).
    """
    residual = horizontal_residual_tm if polarization == "tm" else horizontal_residual_te
    n_side = float(np.sqrt(eps_side).real)
    n_planar = complex(eps_core) ** 0.5

    guesses: list[complex] = list(extra_guesses)
    for delta in (1e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 0.1):
        guesses.append(n_planar - delta)
        guesses.append(complex(n_planar.real - delta, n_planar.imag * 1.2))
    guesses.append(complex((n_planar.real + n_side) / 2.0, n_planar.imag))

    roots: list[complex] = []
    for guess in guesses:
        try:
            if _scipy_root is not None:
                def wrapped(vec: np.ndarray) -> np.ndarray:
                    value = residual(complex(vec[0], vec[1]), k0, eps_core, eps_side, width_um)
                    return np.array([value.real, value.imag])

                sol = _scipy_root(wrapped, np.array([guess.real, guess.imag]), method="hybr", tol=1e-12)
                if not sol.success:
                    continue
                candidate = complex(sol.x[0], sol.x[1])
            else:
                candidate = _newton_scalar(residual, guess, k0, eps_core, eps_side, width_um)
        except Exception:
            continue

        if candidate.imag < 0:
            candidate = candidate.conjugate()
        if not np.isfinite(candidate.real) or not np.isfinite(candidate.imag):
            continue
        if candidate.real <= n_side or candidate.real > n_planar.real + 1e-9:
            continue
        if any(abs(candidate - old) < 1e-9 for old in roots):
            continue
        roots.append(candidate)

    if not roots:
        raise RuntimeError(
            f"связанная мода полоски шириной {width_um} мкм не найдена "
            f"(планарный n_eff = {n_planar:.6f}, боковой показатель = {n_side:.4f})"
        )
    return max(roots, key=lambda z: z.real)


def _newton_scalar(residual, guess, k0, eps_core, eps_side, width_um) -> complex:
    z = complex(guess)
    for _ in range(100):
        value = residual(z, k0, eps_core, eps_side, width_um)
        if abs(value) < 1e-13:
            break
        h = max(1e-9, abs(z) * 1e-8)
        deriv = (residual(z + h, k0, eps_core, eps_side, width_um) - value) / h
        if deriv == 0:
            break
        z = z - value / deriv
    return z


def strip_cutoff_width(
    k0: float, eps_core: complex, eps_side: complex, w_lo: float = 0.05, w_hi: float = 20.0
) -> float:
    """Оценивает минимальную ширину, при которой полоска ещё ведёт моду."""
    lo, hi = w_lo, w_hi
    try:
        solve_strip(k0, eps_core, eps_side, hi)
    except RuntimeError:
        return float("nan")
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        try:
            solve_strip(k0, eps_core, eps_side, mid)
            hi = mid
        except RuntimeError:
            lo = mid
    return hi
