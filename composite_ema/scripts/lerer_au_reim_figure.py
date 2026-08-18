"""Рисунок для Au в формате самого Лерера: Re и Im eps_eff, три модели.

Повторяет его рисунок из «Сравнени3 формул расчета эпс.docx»: шесть кривых на
одном поле, конвенция exp(+i w t), поэтому Im < 0. Нумерация и цвета кривых
повторяют его легенду, чтобы кривые можно было наложить напрямую:

    1  - КМ  - Re (чёрный)      1' - Im (красный)
    2  - Б   - Re (зелёный)     2' - Im (синий)
    3  - МГ  - Re (оливковый)   3' - Im (тёмно-синий)

Строится дважды: на его таблицах luxpop и на Джонсоне-Кристи, чтобы был виден
вклад источника оптических постоянных. Параметры - его: n_m = 1,77, C = 10 %,
R_p = 30 нм в MLWA, диапазон 400-750 нм.

Дополнительно пишется CSV со всеми кривыми для построения в Origin.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import composite_ema as ce
from lerer_compare_diagnosis import bruggeman, eps_luxpop, load_luxpop, mg, mlwa

OUT = Path(__file__).resolve().parent.parent / "results"

N_M = 1.77
EPS_M = complex(N_M**2, 0.0)
C = 0.10
R_P = 30.0
LAM = np.arange(400.0, 750.0 + 1e-9, 1.0)

STYLE = [
    ("1 - КМ - Re", "Re", "black"),
    ("1' -      Im", "Im", "red"),
    ("2 - Б - Re", "Re", "#00c000"),
    ("2' -     Im", "Im", "blue"),
    ("3 - МГ - Re", "Re", "#808000"),
    ("3' -     Im", "Im", "#000080"),
]


def curves(eps_fn) -> dict[str, np.ndarray]:
    eps_p = [eps_fn(float(l)) for l in LAM]
    km = np.array([mg(EPS_M, e, C) for e in eps_p])
    br = np.array([bruggeman(EPS_M, e, C) for e in eps_p])
    ml = np.array([mlwa(EPS_M, e, C, float(l), R_P) for e, l in zip(eps_p, LAM)])
    return {"КМ": km, "Б": br, "МГ": ml}


def draw(ax, data: dict[str, np.ndarray], title: str) -> None:
    order = [("КМ", "black", "red"), ("Б", "#00c000", "blue"), ("МГ", "#808000", "#000080")]
    labels = iter(STYLE)
    for key, c_re, c_im in order:
        vals = data[key]
        lab_re = next(labels)[0]
        ax.plot(LAM, vals.real, color=c_re, linewidth=1.8, label=lab_re)
        lab_im = next(labels)[0]
        ax.plot(LAM, vals.imag, color=c_im, linewidth=1.8, label=lab_im)

    ax.set_xlim(350, 800)
    ax.set_xlabel(r"$\lambda$, nm", fontsize=11)
    ax.set_ylabel(r"$\varepsilon$", fontsize=12, rotation=0, labelpad=12)
    ax.grid(True, alpha=0.25)
    ax.axhline(0.0, color="0.6", linewidth=0.8)
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=8, loc="upper left", framealpha=0.95)
    ax.text(0.06, 0.10, "Au", transform=ax.transAxes, fontsize=15, fontweight="bold")


def main() -> None:
    au_lux = load_luxpop("Au_.c")
    sources = [
        ("данные luxpop (как у А. М. Лерера)", lambda l: eps_luxpop(au_lux, l)),
        ("данные Johnson & Christy 1972", ce.eps_au),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.2))
    all_data = {}
    for ax, (title, fn) in zip(axes, sources):
        data = curves(fn)
        all_data[title] = data
        draw(ax, data, title)

    fig.suptitle(
        r"$\varepsilon_{\mathrm{eff}}$ композита с наночастицами Au: "
        r"$n_m$=1,77, $C$=10 %, MLWA при $R_p$=30 нм. "
        r"Конвенция $\exp(+i\omega t)$, поэтому Im $<$ 0",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    png = OUT / "lerer_au_reim_three_models.png"
    fig.savefig(png, dpi=170)
    print("сохранено: %s" % png)

    # CSV для построения в Origin
    csv_path = OUT / "lerer_au_reim_three_models.csv"
    lux, jc = all_data[sources[0][0]], all_data[sources[1][0]]
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter=";")
        w.writerow([
            "lambda_nm",
            "luxpop_KM_Re", "luxpop_KM_Im",
            "luxpop_Brugg_Re", "luxpop_Brugg_Im",
            "luxpop_MLWA30_Re", "luxpop_MLWA30_Im",
            "JC_KM_Re", "JC_KM_Im",
            "JC_Brugg_Re", "JC_Brugg_Im",
            "JC_MLWA30_Re", "JC_MLWA30_Im",
        ])
        for i, l in enumerate(LAM):
            row = [f"{l:.0f}"]
            for src in (lux, jc):
                for key in ("КМ", "Б", "МГ"):
                    row.append(f"{src[key][i].real:.5f}")
                    row.append(f"{src[key][i].imag:.5f}")
            w.writerow(row)
    print("сохранено: %s" % csv_path)

    # опорные числа
    for name, data in (("luxpop", lux), ("Johnson-Christy", jc)):
        print("\n%s, Au:" % name)
        for key, label in (("КМ", "Клаузиус-Моссотти"), ("Б", "Бруггеман"), ("МГ", "MLWA R=30 нм")):
            v = data[key]
            i_im = int(np.argmax(np.abs(v.imag)))
            i_re = int(np.argmax(v.real))
            print("   %-18s Im: %4.0f нм (%+.3f)   макс Re: %4.0f нм (%+.3f)"
                  % (label, LAM[i_im], v.imag[i_im], LAM[i_re], v.real[i_re]))


if __name__ == "__main__":
    main()
