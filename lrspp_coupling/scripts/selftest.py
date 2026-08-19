"""Самопроверка ядра slabmodes на задачах с известным ответом.

Запуск:
    python lrspp_coupling/scripts/selftest.py

Проверяется пять независимых вещей:
1. общая N-слойная невязка совпадает с трёхслойной формулой из edp/metal_strip_w800,
   уже сверенной с COMSOL;
2. толстая металлическая плёнка даёт ППП одиночной границы с аналитическим n_eff;
3. симметричная плёнка Si|Au|Si воспроизводит корни уравнений Майера (2.29a,b),
   посчитанные независимо вещественным методом бисекции;
4. симметричный диэлектрический слой совпадает со стандартным трансцендентным
   уравнением для чётной TM-моды;
5. восстановленное поле непрерывно на границах, спадает в полупространствах,
   а перекрытие моды с самой собой равно единице.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from slabmodes import (  # noqa: E402
    Stack,
    materials,
    mode_fields,
    orthogonality_residual,
    overlap_power,
    power_flux,
    solve_mode,
    tm_residual,
    tm_residual_3layer,
)
from slabmodes.tmm import decaying_sqrt  # noqa: E402

FAILURES: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    mark = "OK  " if ok else "СБОЙ"
    print(f"  [{mark}] {name}" + (f"  {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


# ---------------------------------------------------------------- проверка 1


def test_three_layer_equivalence() -> None:
    print("1. Общая невязка против трёхслойной формулы edp/metal_strip_w800")
    lam = 1.55
    k0 = materials.k0_from_lambda(lam)
    stack = Stack(
        eps=(materials.EPS_ZPU450, materials.EPS_AU_1550, materials.EPS_ZPU450),
        thickness=(0.020,),
        names=("полимер", "Au", "полимер"),
    )
    worst = 0.0
    for neff in (1.4505 + 1e-5j, 1.46 + 1e-4j, 1.5 + 1e-3j, 2.0 + 0.01j):
        general = tm_residual(neff, stack, k0)
        explicit = tm_residual_3layer(neff, stack, k0)
        g = k0 * decaying_sqrt(neff * neff - stack.eps[1])
        q1 = g / stack.eps[1]
        # общая невязка отличается ровно на множитель q1
        rel = abs(general * q1 - explicit) / max(abs(explicit), 1e-30)
        worst = max(worst, rel)
    check("совпадение с точностью до множителя q1", worst < 1e-10, f"макс. отн. отличие {worst:.2e}")


# ---------------------------------------------------------------- проверка 2


def test_single_interface_limit() -> None:
    print("2. Толстая плёнка -> ППП одиночной границы")
    lam = 1.55
    k0 = materials.k0_from_lambda(lam)
    eps_m = materials.EPS_AU_1550
    eps_d = materials.EPS_ZPU450
    analytic = complex(np.sqrt(eps_m * eps_d / (eps_m + eps_d)))
    if analytic.imag < 0:
        analytic = analytic.conjugate()

    stack = Stack(eps=(eps_d, eps_m, eps_d), thickness=(1.0,))  # 1 мкм - практически полубесконечный
    neff = solve_mode(stack, k0, analytic + 1e-3)
    err = abs(neff - analytic) / abs(analytic)
    check(
        "n_eff сходится к sqrt(eps_m eps_d/(eps_m+eps_d))",
        err < 1e-9,
        f"расчёт {neff.real:.9f}{neff.imag:+.3e}i, аналитика {analytic.real:.9f}{analytic.imag:+.3e}i",
    )


# ---------------------------------------------------------------- проверка 3


def test_maier_symmetric_film() -> None:
    print("3. Симметричная плёнка Si|Au|Si против уравнений Майера (2.29a,b)")
    lam0 = 1.55
    k0 = materials.k0_from_lambda(lam0)
    eps_d = 3.47**2
    eps_m = -118.89142236  # без потерь, как в заметке ваута про точные профили ДМД
    t = 0.020
    a = t / 2.0

    def kap(neff: float, eps: float) -> float:
        return k0 * np.sqrt(neff * neff - eps)

    def maier_a(neff: float) -> float:  # H_y чётная
        return np.tanh(kap(neff, eps_m) * a) + kap(neff, eps_d) * eps_m / (kap(neff, eps_m) * eps_d)

    def maier_b(neff: float) -> float:  # H_y нечётная
        return np.tanh(kap(neff, eps_m) * a) + kap(neff, eps_m) * eps_d / (kap(neff, eps_d) * eps_m)

    def bisect(f, lo: float, hi: float) -> float:
        flo = f(lo)
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            fm = f(mid)
            if flo * fm <= 0:
                hi = mid
            else:
                lo, flo = mid, fm
        return 0.5 * (lo + hi)

    def scan(f) -> float:
        grid = np.linspace(3.4701, 30.0, 40000)
        prev = f(grid[0])
        for x in grid[1:]:
            cur = f(x)
            if np.isfinite(prev) and np.isfinite(cur) and prev * cur < 0:
                return bisect(f, x - (grid[1] - grid[0]), x)
            prev = cur
        raise RuntimeError("корень не найден")

    stack = Stack(eps=(complex(eps_d), complex(eps_m), complex(eps_d)), thickness=(t,))
    for label, func in (("2.29a (H_y чётная)", maier_a), ("2.29b (H_y нечётная)", maier_b)):
        reference = scan(func)
        neff = solve_mode(stack, k0, complex(reference * 1.02, 0.0))
        err = abs(neff.real - reference) / reference
        check(
            f"корень {label}",
            err < 1e-8 and abs(neff.imag) < 1e-9,
            f"матрица переноса {neff.real:.9f}, уравнение Майера {reference:.9f}",
        )


# ---------------------------------------------------------------- проверка 4


def test_dielectric_slab() -> None:
    print("4. Симметричный диэлектрический слой против трансцендентного уравнения")
    lam = 1.55
    k0 = materials.k0_from_lambda(lam)
    n_core, n_clad, d = 1.44, 1.43, 2.8
    eps_core, eps_clad = complex(n_core**2), complex(n_clad**2)

    # чётная TM-мода: (u/eps_core) tan(u d/2) = alpha/eps_clad
    def residual(neff: float) -> float:
        u = k0 * np.sqrt(eps_core.real - neff * neff)
        alpha = k0 * np.sqrt(neff * neff - eps_clad.real)
        return (u / eps_core.real) * np.tan(u * d / 2.0) - alpha / eps_clad.real

    lo, hi = n_clad + 1e-9, n_core - 1e-9
    grid = np.linspace(hi, lo, 20000)
    reference = None
    prev = residual(grid[0])
    for x in grid[1:]:
        cur = residual(x)
        if np.isfinite(prev) and np.isfinite(cur) and prev * cur < 0:
            a, b = x, x + (grid[0] - grid[1])
            for _ in range(200):
                m = 0.5 * (a + b)
                if residual(a) * residual(m) <= 0:
                    b = m
                else:
                    a = m
            reference = 0.5 * (a + b)
            break
        prev = cur
    if reference is None:
        check("найден опорный корень", False)
        return

    stack = Stack(eps=(eps_clad, eps_core, eps_clad), thickness=(d,))
    neff = solve_mode(stack, k0, complex(reference + 1e-4, 0.0))
    err = abs(neff.real - reference)
    check(
        "основная чётная TM-мода",
        err < 1e-10,
        f"матрица переноса {neff.real:.10f}, трансцендентное {reference:.10f}",
    )


# ---------------------------------------------------------------- проверка 5


def test_field_reconstruction() -> None:
    print("5. Восстановление поля: непрерывность, затухание, нормировка, ортогональность")
    lam = 1.55
    k0 = materials.k0_from_lambda(lam)
    stack = Stack(
        eps=(materials.EPS_ZPU450, materials.EPS_AU_1550, materials.EPS_ZPU450),
        thickness=(0.014,),
    )
    neff = solve_mode(stack, k0, 1.4512 + 1e-5j)

    # непрерывность проверяем аналитически: значения по обе стороны границы,
    # а не конечной разностью, иначе измеряется быстрое изменение поля в металле
    edges = stack.interfaces()
    eps_side = 1e-9
    worst_h, worst_e = 0.0, 0.0
    for edge in edges:
        zz = np.array([edge - eps_side, edge + eps_side])
        h, e, _ = mode_fields(neff, stack, k0, zz)
        worst_h = max(worst_h, abs(h[1] - h[0]) / max(abs(h[0]), 1e-30))
        worst_e = max(worst_e, abs(e[1] - e[0]) / max(abs(e[0]), 1e-30))
    check("H_y непрерывна на границах", worst_h < 1e-6, f"макс. разрыв {worst_h:.2e}")
    check("E_x непрерывна на границах", worst_e < 1e-6, f"макс. разрыв {worst_e:.2e}")

    # окно интегрирования привязываем к длине спадания, а не к фиксированным мкм
    gamma_clad = k0 * np.sqrt(neff * neff - stack.eps[0])
    decay_um = 1.0 / abs(gamma_clad.real)
    half = 14.0 * decay_um
    z = np.linspace(-half, half + stack.total_thickness, 200001)
    hy, ex, ez = mode_fields(neff, stack, k0, z)
    eps_z = stack.eps_at(z)

    tail = abs(hy[0]) / abs(hy).max()
    check(
        "поле спадает на краю окна",
        tail < 1e-5,
        f"|H| на краю / макс = {tail:.2e} при окне {half:.1f} мкм",
    )

    p = power_flux(hy, neff, eps_z, z)
    check("продольный поток мощности положителен", p > 0, f"P = {p:.6e}")

    eta = overlap_power((hy, neff, eps_z), (hy, neff, eps_z), z)
    check("перекрытие моды с собой равно 1", abs(eta - 1.0) < 1e-9, f"eta = {eta:.12f}")

    # ортогональность разных мод одной структуры: несопряжённая форма
    stack_sym = Stack(
        eps=(complex(3.47**2), complex(-118.89142236), complex(3.47**2)), thickness=(0.020,)
    )
    n_lr = solve_mode(stack_sym, k0, 3.55 + 0j)
    n_sr = solve_mode(stack_sym, k0, 4.40 + 0j)
    check(
        "найдены две различные моды плёнки",
        abs(n_lr - n_sr) > 1e-3,
        f"n_LR = {n_lr.real:.6f}, n_SR = {n_sr.real:.6f}",
    )
    if abs(n_lr - n_sr) > 1e-3:
        zz = np.linspace(-2.0, 2.0 + stack_sym.total_thickness, 400001)
        h_lr = mode_fields(n_lr, stack_sym, k0, zz)[0]
        h_sr = mode_fields(n_sr, stack_sym, k0, zz)[0]
        eps_zz = stack_sym.eps_at(zz)
        resid = orthogonality_residual(h_lr, h_sr, eps_zz, zz)
        check(
            "две моды несопряжённо ортогональны",
            resid < 1e-6,
            f"|int H1 H2 / eps| / нормы = {resid:.2e}",
        )


def main() -> int:
    print("Самопроверка ядра slabmodes\n")
    test_three_layer_equivalence()
    test_single_interface_limit()
    test_maier_symmetric_film()
    test_dielectric_slab()
    test_field_reconstruction()
    print()
    if FAILURES:
        print(f"ПРОВАЛЕНО проверок: {len(FAILURES)}")
        for name in FAILURES:
            print(f"  - {name}")
        return 1
    print("Все проверки пройдены.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
