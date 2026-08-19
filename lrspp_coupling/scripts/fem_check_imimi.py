"""Полновекторная проверка потерь ДМДМД против ДМД открытым решателем femwell.

Зачем. При воспроизведении Park 2011 отношение потерь двух волноводов не сошлось
с измерениями: расчёт даёт для ДМДМД шириной 3 мкм 10.05 дБ/см против ДМД 0.79,
то есть в 12.7 раза, тогда как измерено 5.01 против 1.89, то есть в 2.65 раза.
Причём у ДМД расчёт НИЖЕ измерения (обычная картина для реального металла), а у
ДМДМД - ВЫШЕ, то есть ошибки идут в разные стороны и одной технологической
поправкой не объясняются.

Проверено заранее и исключено:
  - пропущенная мода: двумерный скан по комплексной плоскости находит ровно одну
    связанную моду ДМДМД, вторая супермода вытеснена под показатель обкладки;
  - толщина центрального слоя: чтобы получить измеренное отношение, она должна
    быть около 2.1 мкм вместо заявленных 500 нм.

Здесь двумерная задача решается методом конечных элементов без приближения
эффективного показателя, чтобы понять, ошибка это редуцированной модели или нет.

Запуск:
    python lrspp_coupling/scripts/fem_check_imimi.py
"""

from __future__ import annotations

import csv
import sys
from collections import OrderedDict
from pathlib import Path

import numpy as np
from shapely.geometry import box

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from slabmodes import Stack, materials, propagation_loss_db_per_cm, solve_mode, solve_strip  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "results"
OUT.mkdir(parents=True, exist_ok=True)

LAMBDA_UM = 1.55
K0 = materials.k0_from_lambda(LAMBDA_UM)
EPS_CLAD = materials.EPS_ZPU450
EPS_AU = materials.EPS_AU_1550
N_CLAD = float(np.sqrt(EPS_CLAD).real)
T_AU_UM = 0.014
T_GAP_UM = 0.500

WIDTH_UM = 3.0
DOMAIN = 30.0

PAPER = {"imi": 1.89, "imimi": 5.01}  # дБ/см, измерено при ширине 3 мкм


def analytic() -> dict[str, complex]:
    imi = Stack(eps=(EPS_CLAD, EPS_AU, EPS_CLAD), thickness=(T_AU_UM,))
    imimi = Stack(eps=(EPS_CLAD, EPS_AU, EPS_CLAD, EPS_AU, EPS_CLAD),
                  thickness=(T_AU_UM, T_GAP_UM, T_AU_UM))
    n_imi = solve_mode(imi, K0, complex(1.4512, 2e-5))
    n_imimi = solve_mode(imimi, K0, complex(1.4534, 7e-5))
    return {
        "imi_planar": n_imi,
        "imimi_planar": n_imimi,
        "imi_strip": solve_strip(K0, n_imi**2, EPS_CLAD, WIDTH_UM),
        "imimi_strip": solve_strip(K0, n_imimi**2, EPS_CLAD, WIDTH_UM),
    }


def fem(kind: str, n_guess: float):
    from femwell.maxwell.waveguide import compute_modes
    from femwell.mesh import mesh_from_OrderedDict
    from skfem import Basis, ElementTriP0
    from skfem.io.meshio import from_meshio

    half = WIDTH_UM / 2.0
    clad = box(-DOMAIN / 2, -DOMAIN / 2, DOMAIN / 2, DOMAIN / 2)
    if kind == "imi":
        metals = OrderedDict(m1=box(-half, -T_AU_UM / 2, half, T_AU_UM / 2))
    else:
        c = T_GAP_UM / 2 + T_AU_UM / 2
        metals = OrderedDict(
            m1=box(-half, c - T_AU_UM / 2, half, c + T_AU_UM / 2),
            m2=box(-half, -c - T_AU_UM / 2, half, -c + T_AU_UM / 2),
        )

    polygons = OrderedDict(**metals, clad=clad)
    resolutions = {k: {"resolution": 0.006, "distance": 0.12} for k in metals}
    resolutions["clad"] = {"resolution": 0.8, "distance": 4.0}
    mesh = from_meshio(mesh_from_OrderedDict(polygons, resolutions, default_resolution_max=2.5))

    basis0 = Basis(mesh, ElementTriP0())
    eps = basis0.zeros(dtype=complex)
    eps[basis0.get_dofs(elements="clad")] = EPS_CLAD
    for k in metals:
        eps[basis0.get_dofs(elements=k)] = EPS_AU

    modes = compute_modes(basis0, eps, wavelength=LAMBDA_UM, num_modes=4, order=1, n_guess=n_guess)
    out = []
    for m in modes:
        n = complex(m.n_eff)
        out.append(n.conjugate() if n.imag < 0 else n)
    out.sort(key=lambda z: -z.real)
    return out, mesh.p.shape[1]


def main() -> int:
    a = analytic()
    lines: list[str] = []
    add = lines.append
    add("Потери ДМД против ДМДМД: метод эффективного показателя и полновекторный FEM")
    add(f"Ширина полоски {WIDTH_UM:g} мкм, золото {T_AU_UM * 1000:.0f} нм, "
        f"центральный слой ДМДМД {T_GAP_UM * 1000:.0f} нм, область {DOMAIN:g} мкм")
    add("")

    rows = []
    for kind, key in (("imi", "imi"), ("imimi", "imimi")):
        n_eim = a[f"{key}_strip"]
        try:
            fem_modes, nodes = fem(kind, float(n_eim.real))
        except Exception as exc:
            add(f"{kind}: FEM не сошёлся: {type(exc).__name__}: {exc}")
            continue
        n_fem = min(fem_modes, key=lambda z: abs(z - n_eim))
        bound = n_fem.real > N_CLAD
        decay = (1.0 / (K0 * float(np.sqrt(max(n_fem.real**2 - N_CLAD**2, 1e-12))))
                 if bound else float("inf"))
        add(f"{kind.upper()}   (узлов {nodes})")
        add(f"   планарный предел      {a[f'{key}_planar'].real:.6f}   "
            f"{propagation_loss_db_per_cm(a[f'{key}_planar'], LAMBDA_UM):8.3f} дБ/см")
        add(f"   ЭДП, полоска          {n_eim.real:.6f}   "
            f"{propagation_loss_db_per_cm(n_eim, LAMBDA_UM):8.3f} дБ/см")
        add(f"   FEM, полновекторный   {n_fem.real:.6f}   "
            f"{propagation_loss_db_per_cm(n_fem, LAMBDA_UM):8.3f} дБ/см")
        add(f"   измерено в статье                    {PAPER[key]:8.2f} дБ/см")
        add(f"   длина спадания поля {decay:.1f} мкм при полуширине области {DOMAIN / 2:.0f} мкм"
            + ("" if bound else "   ВНИМАНИЕ: решение не связано"))
        add(f"   моды FEM: " + ", ".join(f"{m.real:.6f}" for m in fem_modes))
        add("")
        rows.append({"kind": kind, "eim": n_eim, "fem": n_fem, "nodes": nodes, "decay": decay})

    if len(rows) == 2:
        eim_ratio = (propagation_loss_db_per_cm(rows[1]["eim"], LAMBDA_UM)
                     / propagation_loss_db_per_cm(rows[0]["eim"], LAMBDA_UM))
        fem_ratio = (propagation_loss_db_per_cm(rows[1]["fem"], LAMBDA_UM)
                     / propagation_loss_db_per_cm(rows[0]["fem"], LAMBDA_UM))
        add("Отношение потерь ДМДМД к ДМД при ширине 3 мкм")
        add(f"   ЭДП        {eim_ratio:6.2f}")
        add(f"   FEM        {fem_ratio:6.2f}")
        add(f"   измерено   {PAPER['imimi'] / PAPER['imi']:6.2f}")
        with (OUT / "fem_check_imimi.csv").open("w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["kind", "eim_re", "eim_im", "fem_re", "fem_im",
                        "eim_db_cm", "fem_db_cm", "measured_db_cm", "mesh_nodes"])
            for r in rows:
                w.writerow([r["kind"], f"{r['eim'].real:.9f}", f"{r['eim'].imag:.6e}",
                            f"{r['fem'].real:.9f}", f"{r['fem'].imag:.6e}",
                            f"{propagation_loss_db_per_cm(r['eim'], LAMBDA_UM):.4f}",
                            f"{propagation_loss_db_per_cm(r['fem'], LAMBDA_UM):.4f}",
                            PAPER[r["kind"]], r["nodes"]])

    (OUT / "fem_check_imimi.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    for line in lines:
        print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
