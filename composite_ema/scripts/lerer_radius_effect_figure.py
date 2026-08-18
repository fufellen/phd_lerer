"""Рисунок к ответу Лереру: роль радиуса частицы в модели MLWA.

Показывает, что вывод «Бруггеман и MLWA близки, Клаузиус-Моссотти отличается»
держится на значении R_p = 30 нм, переданном в COMPOSITE_MLWA как радиус.
При R_p = 15 нм (если 30 нм имелись в виду как диаметр) картина меняется на
противоположную: MLWA идёт рядом с Клаузиусом-Моссотти.

Данные - таблицы luxpop самого Лерера (Au_.c, Ag_.c).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from lerer_compare_diagnosis import bruggeman, eps_luxpop, load_luxpop, mg, mlwa

OUT = Path(__file__).resolve().parent.parent / "results"


def main() -> None:
    n_m = 1.77
    eps_m = complex(n_m**2, 0.0)
    C = 0.10
    lam = np.arange(400.0, 751.0, 1.0)

    tables = {"Au": load_luxpop("Au_.c"), "Ag": load_luxpop("Ag_.c")}

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4))

    for ax, (mat, table) in zip(axes, tables.items()):
        eps_p = np.array([eps_luxpop(table, float(l)) for l in lam])

        curves = [
            ("Клаузиус-Моссотти", np.array([mg(eps_m, e, C) for e in eps_p]), "k", "-", 2.0),
            ("Бруггеман", np.array([bruggeman(eps_m, e, C) for e in eps_p]), "tab:red", "-", 2.0),
            ("MLWA, $R_p$=30 нм (как в расчёте)",
             np.array([mlwa(eps_m, e, C, float(l), 30.0) for e, l in zip(eps_p, lam)]),
             "tab:green", "-", 2.0),
            ("MLWA, $R_p$=15 нм (диаметр 30 нм)",
             np.array([mlwa(eps_m, e, C, float(l), 15.0) for e, l in zip(eps_p, lam)]),
             "tab:blue", "--", 2.0),
        ]

        for label, vals, color, ls, lw in curves:
            ax.plot(lam, np.abs(vals.imag), color=color, linestyle=ls, linewidth=lw, label=label)
            i = int(np.argmax(np.abs(vals.imag)))
            ax.plot(lam[i], np.abs(vals.imag)[i], marker="o", color=color, markersize=5)

        ax.set_title("%s, $n_m$=1,77, $C$=10 %%" % mat)
        ax.set_xlabel(r"$\lambda$, нм")
        ax.set_ylabel(r"$|\mathrm{Im}\,\varepsilon_{\mathrm{eff}}|$")
        ax.grid(True, alpha=0.3)
        ax.set_xlim(400, 750)

    axes[0].legend(fontsize=8, loc="upper left")
    fig.suptitle(
        "Положение резонанса в трёх моделях (данные luxpop). "
        "Радиус частицы решает, к какой кривой прижмётся MLWA",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))

    OUT.mkdir(exist_ok=True)
    png = OUT / "lerer_radius_effect_AuAg.png"
    fig.savefig(png, dpi=170)
    print("сохранено: %s" % png)

    # печатаем пики для подписи к рисунку
    for mat, table in tables.items():
        eps_p = [eps_luxpop(table, float(l)) for l in lam]
        print("\n%s:" % mat)
        for label, fn in (
            ("Клаузиус-Моссотти", lambda e, l: mg(eps_m, e, C)),
            ("Бруггеман", lambda e, l: bruggeman(eps_m, e, C)),
            ("MLWA R=30", lambda e, l: mlwa(eps_m, e, C, l, 30.0)),
            ("MLWA R=15", lambda e, l: mlwa(eps_m, e, C, l, 15.0)),
        ):
            vals = np.array([fn(e, float(l)) for e, l in zip(eps_p, lam)])
            i = int(np.argmax(np.abs(vals.imag)))
            print("   %-20s пик %4.0f нм, |Im| = %6.3f" % (label, lam[i], abs(vals.imag)[i]))


if __name__ == "__main__":
    main()
