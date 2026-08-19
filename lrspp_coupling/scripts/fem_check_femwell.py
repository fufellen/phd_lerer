"""Полновекторная проверка метода эффективного показателя открытым FEM-решателем.

COMSOL на этом ПК недоступен, поэтому контрольный расчёт выполняется открытым
пакетом femwell (FEM-решатель мод фотонных волноводов на scikit-fem и gmsh,
https://github.com/HelgeGehring/femwell).

Проверяется самое слабое звено аналитической модели - второй шаг метода
эффективного показателя, то есть переход от бесконечной плёнки к полоске
конечной ширины. Именно этим шагом получены числа, на которых держатся
воспроизведения обеих статей Park.

Сравниваются три расчёта для полоски золота 14 нм в полимере ZPU450:
  1) точная планарная задача (матрица переноса) - предел бесконечной ширины;
  2) метод эффективного показателя для конечной ширины;
  3) полновекторный двумерный FEM (femwell).

Запуск:
    python lrspp_coupling/scripts/fem_check_femwell.py
"""

from __future__ import annotations

import csv
import sys
from collections import OrderedDict
from pathlib import Path

import numpy as np
from shapely.geometry import box

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from slabmodes import materials, propagation_loss_db_per_cm, solve_mode, solve_strip  # noqa: E402
from slabmodes.tmm import Stack  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "results"
OUT.mkdir(parents=True, exist_ok=True)

LAMBDA_UM = 1.55
K0 = materials.k0_from_lambda(LAMBDA_UM)
EPS_CLAD = materials.EPS_ZPU450
EPS_AU = materials.EPS_AU_1550
N_CLAD = float(np.sqrt(EPS_CLAD).real)
T_AU_UM = 0.014

WIDTHS = (3.0, 6.0)
DOMAIN_W = 70.0
DOMAIN_H = 70.0


def analytic(width_um: float) -> tuple[complex, complex]:
    """Планарный предел и результат метода эффективного показателя."""
    stack = Stack(eps=(EPS_CLAD, EPS_AU, EPS_CLAD), thickness=(T_AU_UM,))
    n_planar = solve_mode(stack, K0, complex(1.4512, 2e-5))
    n_strip = solve_strip(K0, n_planar**2, EPS_CLAD, width_um)
    return n_planar, n_strip


def fem(width_um: float, n_modes: int = 4, refine: float = 1.0) -> list[complex]:
    """Полновекторный расчёт двумерного сечения в femwell."""
    from femwell.maxwell.waveguide import compute_modes
    from femwell.mesh import mesh_from_OrderedDict
    from skfem import Basis, ElementTriP0
    from skfem.io.meshio import from_meshio

    core = box(-width_um / 2, -T_AU_UM / 2, width_um / 2, T_AU_UM / 2)
    clad = box(-DOMAIN_W / 2, -DOMAIN_H / 2, DOMAIN_W / 2, DOMAIN_H / 2)

    polygons = OrderedDict(core=core, clad=clad)
    resolutions = {
        "core": {"resolution": 0.004 * refine, "distance": 0.6},
        "clad": {"resolution": 1.0 * refine, "distance": 8.0},
    }
    mesh = from_meshio(
        mesh_from_OrderedDict(polygons, resolutions, default_resolution_max=4.0 * refine)
    )

    basis0 = Basis(mesh, ElementTriP0())
    epsilon = basis0.zeros(dtype=complex)
    epsilon[basis0.get_dofs(elements="clad")] = EPS_CLAD
    epsilon[basis0.get_dofs(elements="core")] = EPS_AU

    modes = compute_modes(
        basis0, epsilon, wavelength=LAMBDA_UM, num_modes=n_modes,
        order=1, n_guess=1.4515,
    )
    out = []
    for m in modes:
        n = complex(m.n_eff)
        if n.imag < 0:
            n = n.conjugate()
        out.append(n)
    out.sort(key=lambda z: -z.real)
    return out, mesh.p.shape[1]


def main() -> int:
    lines: list[str] = []
    add = lines.append
    add("Полновекторная проверка ЭДП открытым решателем femwell")
    add(f"Полоска золота {T_AU_UM * 1000:.0f} нм в полимере n = {N_CLAD:.3f}, "
        f"lambda = {LAMBDA_UM} мкм")
    add(f"Расчётная область {DOMAIN_W:g} x {DOMAIN_H:g} мкм")
    add("")

    rows = []
    for width in WIDTHS:
        n_planar, n_strip = analytic(width)
        try:
            fem_modes, n_nodes = fem(width)
        except Exception as exc:  # pragma: no cover
            add(f"FEM для ширины {width} мкм не сошёлся: {type(exc).__name__}: {exc}")
            continue

        # выбираем моду, ближайшую к аналитической: у полоски есть также
        # краевые и высшие решения, номер моды сам по себе ничего не значит
        n_fem = min(fem_modes, key=lambda z: abs(z - n_strip))

        add(f"Ширина {width:g} мкм   (узлов сетки {n_nodes})")
        add(f"   планарный предел      n_eff = {n_planar.real:.6f}{n_planar.imag:+.3e}i"
            f"   {propagation_loss_db_per_cm(n_planar, LAMBDA_UM):8.3f} дБ/см")
        add(f"   ЭДП, конечная ширина  n_eff = {n_strip.real:.6f}{n_strip.imag:+.3e}i"
            f"   {propagation_loss_db_per_cm(n_strip, LAMBDA_UM):8.3f} дБ/см")
        add(f"   FEM, полновекторный   n_eff = {n_fem.real:.6f}{n_fem.imag:+.3e}i"
            f"   {propagation_loss_db_per_cm(n_fem, LAMBDA_UM):8.3f} дБ/см")
        d_re = n_strip.real - n_fem.real
        d_loss = (propagation_loss_db_per_cm(n_strip, LAMBDA_UM)
                  / max(propagation_loss_db_per_cm(n_fem, LAMBDA_UM), 1e-30))
        add(f"   отличие ЭДП от FEM: Re n_eff {d_re:+.6f}, потери в {d_loss:.2f} раза")
        add(f"   все найденные FEM-моды: " + ", ".join(f"{m.real:.6f}" for m in fem_modes))
        if n_fem.real <= N_CLAD:
            add("   ВНИМАНИЕ: FEM-мода лежит ниже показателя обкладки - решение "
                "ограничено размером области, а не физикой; нужна большая область")
        else:
            decay = 1.0 / (K0 * float(np.sqrt(max(n_fem.real**2 - N_CLAD**2, 1e-12))))
            add(f"   длина спадания поля в обкладке {decay:.1f} мкм при полуширине "
                f"области {DOMAIN_W / 2:.0f} мкм")
        add("")
        rows.append({
            "width": width, "planar": n_planar, "eim": n_strip, "fem": n_fem,
            "d_re": d_re, "loss_ratio": d_loss, "nodes": n_nodes,
        })

    if rows:
        with (OUT / "fem_check_femwell.csv").open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["width_um", "planar_re", "planar_im", "eim_re", "eim_im",
                        "fem_re", "fem_im", "delta_re", "loss_ratio_eim_over_fem", "mesh_nodes"])
            for r in rows:
                w.writerow([f"{r['width']:.2f}",
                            f"{r['planar'].real:.9f}", f"{r['planar'].imag:.6e}",
                            f"{r['eim'].real:.9f}", f"{r['eim'].imag:.6e}",
                            f"{r['fem'].real:.9f}", f"{r['fem'].imag:.6e}",
                            f"{r['d_re']:.6e}", f"{r['loss_ratio']:.4f}", r["nodes"]])

    (OUT / "fem_check_femwell.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    for line in lines:
        print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
