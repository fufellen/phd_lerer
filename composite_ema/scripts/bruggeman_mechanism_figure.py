"""Рисунок к ответу: почему Бруггеман так далеко от Максвелла-Гарнетта.

Левая панель - обе модели при разных концентрациях: при малых C кривые
совпадают (модели тождественны в первом порядке по C), расхождение растёт
вместе с C.

Правая панель - причина: модуль резонансного знаменателя. У Максвелла-Гарнетта
включение сидит в прозрачной матрице, |eps_p + 2 eps_m| проваливается до 2,0 -
резонанс острый. У Бруггемана фон - сама смесь, комплексная и поглощающая,
|eps_p + 2 eps_eff| нигде не опускается ниже 5,2 - острому резонансу взяться
неоткуда.

Данные - таблицы luxpop А. М. Лерера, Au.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Рисунок идёт в письмо на ширину колонки A4, поэтому шрифты крупнее обычного
matplotlib.rcParams.update({
    "font.size": 11,
    "axes.titlesize": 11.5,
    "axes.labelsize": 11.5,
    "xtick.labelsize": 10.5,
    "ytick.labelsize": 10.5,
})

from lerer_compare_diagnosis import bruggeman, eps_luxpop, load_luxpop, mg

OUT = Path(__file__).resolve().parent.parent / "results"

N_M = 1.77
EPS_M = complex(N_M**2, 0.0)


def main() -> None:
    au = load_luxpop("Au_.c")
    lam = np.arange(400.0, 1000.0 + 1e-9, 1.0)
    eps_p = np.array([eps_luxpop(au, float(l)) for l in lam])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.2, 4.8))

    # --- панель 1: сходимость при малых C -------------------------------
    colors = {0.01: "#1f77b4", 0.05: "#2ca02c", 0.10: "#d62728", 0.20: "#9467bd"}
    for C, col in colors.items():
        km = np.array([mg(EPS_M, e, C) for e in eps_p])
        br = np.array([bruggeman(EPS_M, e, C) for e in eps_p])
        ax1.plot(lam, np.abs(km.imag), color=col, linewidth=2.0,
                 label="Клаузиус-Моссотти, $C$=%g %%" % (C * 100))
        ax1.plot(lam, np.abs(br.imag), color=col, linewidth=2.0, linestyle="--",
                 label="Бруггеман, $C$=%g %%" % (C * 100))

    ax1.set_xlabel(r"$\lambda$, нм")
    ax1.set_ylabel(r"$|\mathrm{Im}\,\varepsilon_{\mathrm{eff}}|$")
    ax1.set_title("Модели тождественны при малых $C$\nи расходятся с ростом $C$", fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=8.5, ncol=1, loc="upper right")
    ax1.set_xlim(400, 1000)
    ax1.set_ylim(0, 12)

    # --- панель 2: резонансный знаменатель ------------------------------
    C = 0.10
    br = np.array([bruggeman(EPS_M, e, C) for e in eps_p])
    den_mg = np.abs(eps_p + 2.0 * EPS_M)
    den_br = np.abs(eps_p + 2.0 * br)

    ax2.plot(lam, den_mg, color="black", linewidth=2.2,
             label=r"МГ: $|\varepsilon_p + 2\varepsilon_m|$, фон прозрачный")
    ax2.plot(lam, den_br, color="#d62728", linewidth=2.2,
             label=r"Б: $|\varepsilon_p + 2\varepsilon_{\mathrm{eff}}|$, фон поглощающий")

    i1, i2 = int(np.argmin(den_mg)), int(np.argmin(den_br))
    ax2.plot(lam[i1], den_mg[i1], "o", color="black", markersize=7)
    ax2.plot(lam[i2], den_br[i2], "o", color="#d62728", markersize=7)
    # Подписи разносим в свободную область справа сверху, иначе при печати
    # на ширину колонки они наезжают друг на друга
    ax2.annotate("минимум %.2f (%.0f нм)" % (den_br[i2], lam[i2]),
                 xy=(lam[i2], den_br[i2]), xytext=(660, 20.5),
                 fontsize=10, color="#d62728", fontweight="bold",
                 arrowprops=dict(arrowstyle="->", color="#d62728", lw=1.3))
    ax2.annotate("минимум %.2f (%.0f нм)" % (den_mg[i1], lam[i1]),
                 xy=(lam[i1], den_mg[i1]), xytext=(660, 14.0),
                 fontsize=10, color="black", fontweight="bold",
                 arrowprops=dict(arrowstyle="->", color="black", lw=1.3))

    ax2.set_xlabel(r"$\lambda$, нм")
    ax2.set_ylabel("модуль резонансного знаменателя")
    ax2.set_title("Почему у Бруггемана нет острого резонанса:\n"
                  "знаменатель нигде не становится малым", fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=10, loc="upper left")
    ax2.set_xlim(400, 1000)
    ax2.set_ylim(0, 32)

    fig.suptitle("Композит с наночастицами Au, $n_m$=1,77, данные luxpop", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))

    png = OUT / "bruggeman_vs_mg_mechanism.png"
    fig.savefig(png, dpi=170)
    print("сохранено: %s" % png)


if __name__ == "__main__":
    main()
