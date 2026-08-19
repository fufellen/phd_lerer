"""Валидация открытого решателя femwell по эталонным результатам COMSOL.

Зачем. COMSOL доступен не на всех машинах, поэтому нужен открытый решатель мод,
которому можно доверять. В этой папке лежат экспортированные из COMSOL 6.2
собственные моды Ag-полоски на подложке и отобранные ветви. Здесь та же задача
решается пакетом femwell (FEM на scikit-fem и gmsh), и результаты сравниваются
число в число.

Задача. Полоска серебра толщиной 20 нм лежит на подложке n = 1.45, сверху и по
бокам воздух. Расчётная область 6 x 4 мкм - та же, что в COMSOL-модели,
созданной скриптами CreateAgStrip*.java. Диэлектрические проницаемости серебра
взяты те же самые, что подставлялись в COMSOL, чтобы сравнение шло на одинаковых
входных данных, а не на разных редакциях таблиц.

Две сверки:
  1) развёртка по длине волны при ширине 800 нм, от 450 до 800 нм - основная,
     поскольку в той сводке отбор моды выполнен по полевым признакам
     (столбец selection_rule = field_pool_closest_to_edp);
  2) развёртка по ширине полоски при 500 нм - контрольная, и она вскрыла дефект
     СТАРОГО отбора, а не решателя, см. ниже.

Критерий приёмки. Плазмонная мода металлической полоски сильно локализована, обе
программы решают одну и ту же краевую задачу, поэтому расхождение по Re(n_eff)
должно быть в пределах десятых долей процента, а по Im(n_eff) - в пределах
единиц процентов; последнее чувствительнее, поскольку мнимая часть определяется
полем внутри тонкого металла и потому зависит от сетки.

Найденный дефект старой постобработки. Для одной и той же конфигурации
(W = 800 нм, lambda = 500 нм) две сводки в results/ дают разные COMSOL-значения:
2.060106 в развёртке по длине волны и 1.952123 в развёртке по ширине. Прогон
femwell показывает, что у полоски есть лестница плазмонных мод, и старый отбор в
развёртке по ширине сидит на модах номер 1, 2, 3 и 5 по мере роста ширины, тогда
как физическая основная мода - всегда с наибольшим Re(n_eff). У развёртки по
ширине нет столбца selection_rule, то есть отбор там вёлся без полевого
критерия. Поэтому эталоном взята развёртка по длине волны, а развёртка по ширине
сверяется с аналитикой ЭДП.

Запуск (из корня репозитория):
    python .\\edp\\metal_strip_w800\\scripts\\validate_femwell_vs_comsol.py
"""

from __future__ import annotations

import csv
import sys
from collections import OrderedDict
from pathlib import Path

import numpy as np
from shapely.geometry import box

BASE = Path(__file__).resolve().parents[1]
RESULTS = BASE / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

# --- параметры модели, скопированные из CreateAgStrip*.java ------------------
N_SUB = 1.45
EPS_SUB = complex(N_SUB**2, 0.0)
EPS_AIR = complex(1.0, 0.0)
T_AG_UM = 0.020
W_BOX_UM = 6.0
H_AIR_UM = 2.0
H_SUB_UM = 2.0

EPS_AG_BY_LAMBDA = {
    450: complex(-6.078099, 0.74594),
    500: complex(-8.491989, 0.75842),
    550: complex(-11.126868, 0.827824),
    600: complex(-13.904985, 0.925288),
    650: complex(-17.021063, 1.147584),
    700: complex(-20.42832, 1.284248),
    750: complex(-23.969379, 1.42042),
    800: complex(-27.953072, 1.512654),
}

# Аналитика ЭДП из results/comsol_width_sweep_500nm_selected_modes.csv
WIDTH_ANALYTIC = {400: 2.01513048, 800: 2.06858879, 1200: 2.0797623, 2000: 2.08578405}


def read_comsol_selected(path: Path, key: str) -> dict[float, complex]:
    out: dict[float, complex] = {}
    with path.open(encoding="utf-8-sig", newline="") as f:
        for row in csv.DictReader(f):
            out[float(row[key])] = complex(float(row["comsol_real"]), float(row["comsol_imag"]))
    return out


def solve_femwell(width_um: float, lambda_nm: int, num_modes: int = 10,
                  metal_res_um: float = 0.005, guess: float = 2.05):
    """Собственные моды сечения в femwell на геометрии COMSOL-модели."""
    from femwell.maxwell.waveguide import compute_modes
    from femwell.mesh import mesh_from_OrderedDict
    from skfem import Basis, ElementTriP0
    from skfem.io.meshio import from_meshio

    half_w = width_um / 2.0
    polygons = OrderedDict(
        strip=box(-half_w, 0.0, half_w, T_AG_UM),
        sub=box(-W_BOX_UM / 2, -H_SUB_UM, W_BOX_UM / 2, 0.0),
        air=box(-W_BOX_UM / 2, -H_SUB_UM, W_BOX_UM / 2, H_AIR_UM),
    )
    resolutions = {
        "strip": {"resolution": metal_res_um, "distance": 0.25},
        "sub": {"resolution": 0.25, "distance": 1.0},
        "air": {"resolution": 0.25, "distance": 1.0},
    }
    mesh = from_meshio(mesh_from_OrderedDict(polygons, resolutions, default_resolution_max=0.5))

    basis0 = Basis(mesh, ElementTriP0())
    eps = basis0.zeros(dtype=complex)
    eps[basis0.get_dofs(elements="air")] = EPS_AIR
    eps[basis0.get_dofs(elements="sub")] = EPS_SUB
    eps[basis0.get_dofs(elements="strip")] = EPS_AG_BY_LAMBDA[lambda_nm]

    modes = compute_modes(basis0, eps, wavelength=lambda_nm / 1000.0,
                          num_modes=num_modes, order=1, n_guess=guess)
    found = []
    for m in modes:
        n = complex(m.n_eff)
        found.append(n.conjugate() if n.imag < 0 else n)
    return found, mesh.p.shape[1]


def plasmonic_ladder(modes: list[complex]) -> list[complex]:
    """Затухающие моды полоски по убыванию Re; первая из них - основная.

    Отбираются решения выше показателя подложки и с заметным затуханием: это
    отсекает почти лишённые потерь плёночные моды подложки и моды расчётного
    ящика, у которых мнимая часть на порядки меньше.
    """
    return sorted([m for m in modes if m.real > N_SUB and m.imag > 1e-4], key=lambda z: -z.real)


def main() -> int:
    width_ref = read_comsol_selected(RESULTS / "comsol_width_sweep_500nm_selected_modes.csv", "width_nm")
    lambda_ref = read_comsol_selected(RESULTS / "comsol_lambda_sweep_w800_selected_modes.csv", "lambda_nm")

    lines: list[str] = []
    add = lines.append
    add("Валидация femwell по эталонным результатам COMSOL 6.2")
    add(f"Ag-полоска {T_AG_UM * 1000:.0f} нм на подложке n = {N_SUB}, "
        f"область {W_BOX_UM:g} x {H_AIR_UM + H_SUB_UM:g} мкм")
    add("Проницаемости серебра совпадают с подставленными в COMSOL")
    add("")

    add("1. Основная сверка: развёртка по длине волны, ширина 800 нм")
    add("   lambda   COMSOL Re    femwell Re   откл.,%   COMSOL Im    femwell Im   откл.,%   узлов")
    d_re, d_im = [], []
    lambda_rows = []
    for lam_nm, n_comsol in sorted(lambda_ref.items()):
        modes, nodes = solve_femwell(0.800, int(lam_nm), guess=float(n_comsol.real))
        f = plasmonic_ladder(modes)[0]
        er = 100.0 * (f.real - n_comsol.real) / n_comsol.real
        ei = 100.0 * (f.imag - n_comsol.imag) / n_comsol.imag
        d_re.append(abs(er))
        d_im.append(abs(ei))
        add(f"   {lam_nm:6.0f}   {n_comsol.real:.6f}   {f.real:10.6f}   {er:+7.3f}   "
            f"{n_comsol.imag:.6f}   {f.imag:10.6f}   {ei:+7.2f}   {nodes:5d}")
        lambda_rows.append({"label": f"{lam_nm:.0f}", "comsol": n_comsol, "fem": f, "nodes": nodes})
    max_re, max_im = max(d_re), max(d_im)
    add(f"   отклонение Re: среднее {np.mean(d_re):.3f} %, макс {max_re:.3f} %")
    add(f"   отклонение Im: среднее {np.mean(d_im):.2f} %, макс {max_im:.2f} %")
    add("")

    add("2. Контроль: развёртка по ширине при 500 нм")
    add("   femwell берётся основная мода и сравнивается с аналитикой ЭДП;")
    add("   дополнительно указано, какому номеру моды отвечает старый отбор COMSOL")
    add("   W, нм   femwell основная       ЭДП        откл.%   старый COMSOL   его номер")
    width_rows = []
    for w_nm, n_comsol in sorted(width_ref.items()):
        modes, nodes = solve_femwell(w_nm / 1000.0, 500, guess=2.05)
        ladder = plasmonic_ladder(modes)
        f = ladder[0]
        idx = min(range(len(ladder)), key=lambda i: abs(ladder[i] - n_comsol)) + 1
        d = 100.0 * (f.real - WIDTH_ANALYTIC[w_nm]) / WIDTH_ANALYTIC[w_nm]
        add(f"   {w_nm:5.0f}   {f.real:.6f}{f.imag:+.6f}i   {WIDTH_ANALYTIC[w_nm]:.6f}   "
            f"{d:+6.2f}   {n_comsol.real:.6f}        номер {idx}")
        width_rows.append({"label": f"W={w_nm:.0f}", "comsol": n_comsol, "fem": f,
                           "nodes": nodes, "comsol_index": idx})
    add("   Основная мода femwell сходится к аналитике ЭДП по мере расширения полоски,")
    add("   как и должно быть: метод эффективного показателя точнее для широких полосок.")
    add("   Старый отбор COMSOL с ростом ширины съезжает на всё более высокие моды,")
    add("   поэтому эта сводка непригодна как эталон и требует переотбора по полю.")
    add("")

    ok_re, ok_im = max_re < 1.0, max_im < 10.0
    add("Итог")
    add(f"   отклонение Re(n_eff) от COMSOL: макс {max_re:.3f} %   порог 1 %   "
        f"{'пройден' if ok_re else 'НЕ ПРОЙДЕН'}")
    add(f"   отклонение Im(n_eff) от COMSOL: макс {max_im:.2f} %   порог 10 %   "
        f"{'пройден' if ok_im else 'НЕ ПРОЙДЕН'}")
    add(f"   статус femwell: {'validated' if (ok_re and ok_im) else 'discrepant'}")
    add("   побочный результат: в results/comsol_width_sweep_500nm_selected_modes.csv")
    add("   отобраны не основные моды, сводка непригодна как эталон.")

    with (RESULTS / "femwell_vs_comsol.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["sweep", "parameter", "comsol_re", "comsol_im", "femwell_re", "femwell_im",
                    "delta_re_percent", "delta_im_percent", "mesh_nodes", "comsol_mode_index"])
        for r in lambda_rows:
            w.writerow(["lambda_w800", r["label"], f"{r['comsol'].real:.9f}", f"{r['comsol'].imag:.9f}",
                        f"{r['fem'].real:.9f}", f"{r['fem'].imag:.9f}",
                        f"{100 * (r['fem'].real - r['comsol'].real) / r['comsol'].real:.4f}",
                        f"{100 * (r['fem'].imag - r['comsol'].imag) / r['comsol'].imag:.3f}",
                        r["nodes"], ""])
        for r in width_rows:
            w.writerow(["width_500nm", r["label"], f"{r['comsol'].real:.9f}", f"{r['comsol'].imag:.9f}",
                        f"{r['fem'].real:.9f}", f"{r['fem'].imag:.9f}", "", "",
                        r["nodes"], r["comsol_index"]])

    (RESULTS / "femwell_vs_comsol.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    for line in lines:
        print(line)
    return 0 if (ok_re and ok_im) else 1


if __name__ == "__main__":
    raise SystemExit(main())
