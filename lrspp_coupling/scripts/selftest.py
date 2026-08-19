"""Самопроверка ядра slabmodes на задачах с известным ответом.

Запуск:
    python lrspp_coupling/scripts/selftest.py

Проверяется шесть независимых вещей:
1. общая N-слойная невязка совпадает с трёхслойной формулой из edp/metal_strip_w800,
   уже сверенной с COMSOL;
2. толстая металлическая плёнка даёт ППП одиночной границы с аналитическим n_eff;
3. симметричная плёнка Si|Au|Si воспроизводит корни уравнений Майера (2.29a,b),
   посчитанные независимо вещественным методом бисекции;
4. симметричный диэлектрический слой совпадает со стандартным трансцендентным
   уравнением для чётной TM-моды;
5. восстановленное поле непрерывно на границах, спадает в полупространствах,
   а перекрытие моды с самой собой равно единице;
6. восстановленные E и H подставляются в уравнения Максвелла: невязка законов
   Фарадея и Ампера должна быть на уровне ошибки конечной разности.
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


# ---------------------------------------------------------------- проверка 6


def test_maxwell_residual() -> None:
    """Прямая подстановка восстановленного поля в уравнения Максвелла.

    В нормировке mode_fields (длины в мкм, E~ = omega eps0 E / 1e6, H~ = H_y)
    уравнения имеют вид rot E~ = i k0^2 H~ и rot H~ = -i eps E~. Для TM-моды,
    зависящей от x как exp(i k0 n_eff x), отличны от нуля только

        (rot E)_y = dE_x/dz - dE_z/dx,   (rot H)_x = -dH_y/dz,
                                         (rot H)_z =  dH_y/dx,

    поэтому проверка сводится к трём числам. Производная по x берётся
    аналитически (i k0 n_eff), по z - центральной разностью строго внутри
    одного слоя: на границе dE_x/dz терпит разрыв, а E_z разрывна сама.

    Именно эта проверка ловит потерю множителя k0 у E_z: без него невязка
    Фарадея составляет около 94 %, а не 1e-7.
    """
    print("6. Невязка уравнений Максвелла для восстановленного поля")
    lam = 1.55
    k0 = materials.k0_from_lambda(lam)
    stack = Stack(
        eps=(materials.EPS_ZPU450, materials.EPS_AU_1550, materials.EPS_ZPU450),
        thickness=(0.014,),
        names=("полимер", "Au", "полимер"),
    )
    neff = solve_mode(stack, k0, 1.4512 + 1e-5j)
    edges = stack.interfaces()

    # точки строго внутри слоёв, с запасом на шаг разности. Середину плёнки
    # брать нельзя: у чётной моды там dH_y/dz = 0 тождественно, и относительная
    # невязка делится на нуль. Поэтому металл зондируется со смещением.
    h = 2e-5
    probes = [
        ("обкладка снизу", -0.60),
        ("металл, четверть толщины", edges[0] + 0.25 * stack.thickness[0]),
        ("обкладка сверху", edges[-1] + 0.60),
    ]

    worst_far, worst_amp_x, worst_amp_z = 0.0, 0.0, 0.0
    for _, z0 in probes:
        zz = np.array([z0 - h, z0, z0 + h])
        hy, ex, ez = mode_fields(neff, stack, k0, zz)
        eps_here = stack.eps_at(np.array([z0]))[0]

        dex_dz = (ex[2] - ex[0]) / (2.0 * h)
        dhy_dz = (hy[2] - hy[0]) / (2.0 * h)
        dez_dx = 1j * k0 * neff * ez[1]
        dhy_dx = 1j * k0 * neff * hy[1]

        # каждое уравнение нормируется на масштаб входящих в него слагаемых,
        # с полом k0|H|: иначе в точке, где слагаемое случайно обращается в
        # нуль, относительная невязка теряет смысл
        floor = k0 * abs(hy[1])

        # Фарадей: (rot E)_y = i k0^2 H_y
        rot_e_y = dex_dz - dez_dx
        scale_e = max(abs(dex_dz), abs(dez_dx), k0 * floor)
        worst_far = max(worst_far, abs(rot_e_y - 1j * k0 * k0 * hy[1]) / scale_e)

        # Ампер: (rot H)_x = -i eps E_x, (rot H)_z = -i eps E_z
        scale_x = max(abs(dhy_dz), abs(eps_here * ex[1]), floor)
        scale_z = max(abs(dhy_dx), abs(eps_here * ez[1]), floor)
        worst_amp_x = max(worst_amp_x, abs(-dhy_dz + 1j * eps_here * ex[1]) / scale_x)
        worst_amp_z = max(worst_amp_z, abs(dhy_dx + 1j * eps_here * ez[1]) / scale_z)

    # порог 1e-6 отвечает ошибке центральной разности: обе эти невязки содержат
    # численную производную по z. Компонента z берёт производную аналитически,
    # поэтому её порог машинный
    check("закон Фарадея выполнен", worst_far < 1e-6, f"макс. невязка {worst_far:.2e}")
    check("закон Ампера, компонента x", worst_amp_x < 1e-6, f"макс. невязка {worst_amp_x:.2e}")
    check("закон Ампера, компонента z", worst_amp_z < 1e-12, f"макс. невязка {worst_amp_z:.2e}")

    # независимая проверка масштаба E_z: в обкладке |E_z / E_x| = beta / kappa_d
    z_probe = np.array([edges[-1] + 0.60])
    hy, ex, ez = mode_fields(neff, stack, k0, z_probe)
    kappa_d = k0 * decaying_sqrt(neff * neff - stack.eps[-1])
    expected = abs(k0 * neff / kappa_d)
    got = abs(ez[0] / ex[0])
    check(
        "отношение |E_z / E_x| в обкладке равно beta / kappa",
        abs(got - expected) / expected < 1e-9,
        f"расчёт {got:.6f}, аналитика {expected:.6f}",
    )


def main() -> int:
    print("Самопроверка ядра slabmodes\n")
    test_three_layer_equivalence()
    test_single_interface_limit()
    test_maier_symmetric_film()
    test_dielectric_slab()
    test_field_reconstruction()
    test_maxwell_residual()
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
