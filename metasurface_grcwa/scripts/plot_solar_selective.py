"""Рисунки к solar_selective.py и ri_sensor_check.py.
Выход: results/solar_selective.png, results/ri_sensor.png."""
from __future__ import annotations

import csv
import io
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import lerer_grcwa as L  # noqa: E402
import solar_selective as S  # noqa: E402

COL = {"original": "#ff7f0e", "improved": "#1f77b4", "improved_no_mirror": "#2ca02c"}
LAB = {"original": "исходная (подложка)", "improved": "предложенная (золотое зеркало)",
       "improved_no_mirror": "предложенная без зеркала (подложка)"}


def read_solar():
    rows = list(csv.DictReader(io.open(L.RESULTS / "solar_selective_spectra.csv", encoding="utf-8")))
    out = {}
    for r in rows:
        out.setdefault((r["design"], r["au_model"]), []).append(
            (float(r["lambda_nm"]), float(r["R"]), float(r["T"]), float(r["A_layers"]), float(r["A_total"])))
    return {k: np.array(sorted(v)) for k, v in out.items()}


def plot_solar(data):
    fig, axs = plt.subplots(1, 3, figsize=(15.5, 7.6))
    wl, g = S.am15g()
    ax = axs[0]
    for d in ("original", "improved_no_mirror", "improved"):
        arr = data[(d, "lerer_fit")]
        sel = arr[:, 0] <= 4000
        a = arr[sel, 4] if d == "improved" else arr[sel, 3]
        ax.plot(arr[sel, 0] / 1000, a, color=COL[d], lw=1.8, label=LAB[d])
    ax2 = ax.twinx()
    ax2.fill_between(wl / 1000, g, color="#f2c94c", alpha=0.35, label="спектр AM1.5G")
    ax2.set_ylabel("спектральная плотность AM1.5G, Вт/(м²·нм)")
    ax2.set_ylim(0, 2.2)
    ax.set_xlim(0.3, 4.0); ax.set_ylim(0, 1.02)
    ax.set_xlabel("длина волны в свободном пространстве, мкм")
    ax.set_ylabel("поглощение (с зеркалом 1 − R, на подложке 1 − R − T)")
    ax.set_title("Солнечная полоса 0,3…4 мкм", fontsize=10)
    ax.grid(alpha=0.3)
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=7.5, loc="upper right")
    ax = axs[1]
    for d in ("original", "improved_no_mirror", "improved"):
        arr = data[(d, "lerer_fit")]
        sel = arr[:, 0] >= 2500
        ax.plot(arr[sel, 0] / 1000, arr[sel, 4], color=COL[d], lw=1.8, label=LAB[d] + ", 1 − R")
        arr2 = data[(d, "lerer_ordal")]
        ax.plot(arr2[:, 0] / 1000, arr2[:, 4], color=COL[d], lw=1.0, ls=":", label=LAB[d] + ", золото по Ordal")
    ax2 = ax.twinx()
    lam = np.geomspace(2500, 50000, 400)
    for T, c in ((373.0, "#d62728"), (573.0, "#8c564b")):
        b = S.planck(lam, T); ax2.plot(lam / 1000, b / b.max(), color=c, lw=1.2, ls="--", label="Планк, %.0f К (норм.)" % T)
    ax2.set_ylabel("спектр Планка, отн. ед.")
    ax2.set_ylim(0, 1.05)
    ax.set_xscale("log"); ax.set_xlim(2.5, 50); ax.set_ylim(0, 1.02)
    ax.set_xlabel("длина волны в свободном пространстве, мкм")
    ax.set_ylabel("излучательная способность ε = 1 − R (закон Кирхгофа)")
    ax.set_title("Тепловой ИК 2,5…50 мкм", fontsize=10)
    ax.grid(alpha=0.3, which="both")
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=1)
    ax = axs[2]
    Ts = np.linspace(320, 800, 25)
    for d in ("original", "improved_no_mirror", "improved"):
        sol = data[(d, "lerer_fit")]
        vis = sol[sol[:, 0] <= 4000]; ir = sol[sol[:, 0] >= 2500]
        a_sol = vis[:, 4] if d == "improved" else vis[:, 3]
        alpha, _ = S.solar_absorptance(vis[:, 0], a_sol)
        lam_all = np.concatenate([vis[:, 0], ir[:, 0]]); a_all = np.concatenate([vis[:, 4], ir[:, 4]])
        for C, ls in ((1, "-"), (10, "--")):
            eta = [alpha - S.emittance(lam_all, a_all, T) * S.SIGMA * T ** 4 / (C * S.Q_SUN) for T in Ts]
            ax.plot(Ts - 273.15, eta, color=COL[d], lw=1.8, ls=ls, label="%s, %d %s" % (LAB[d], C, "солнце" if C == 1 else "солнц"))
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("температура поглотителя, °C")
    ax.set_ylabel("к.п.д. фототермического преобразования η")
    ax.set_ylim(-0.2, 1.0)
    ax.set_title("η = α_s − ε_th σT⁴ / (C·1000 Вт/м²)", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=1)
    fig.suptitle("Столбиковый поглотитель как селективный солнечный поглотитель: нормальное падение, "
                 "матрица и прослойка без потерь в ИК, золото Au.txt + Друде", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(L.RESULTS / "solar_selective.png", dpi=170)
    plt.close(fig)


def plot_ri():
    rows = list(csv.DictReader(io.open(L.RESULTS / "ri_sensor_spectra.csv", encoding="utf-8")))
    data = {}
    for r in rows:
        data.setdefault((r["design"], float(r["n_super"])), []).append(
            (float(r["lambda_nm"]), float(r["R"]), float(r["A_layers"]), float(r["A_total"])))
    data = {k: np.array(v) for k, v in data.items()}
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
    for ax, d in zip(axs, ("original", "improved")):
        for ns, ls in ((1.00, "-"), (1.33, "--"), (1.35, ":")):
            arr = data[(d, ns)]
            ax.plot(arr[:, 0], arr[:, 1], color="#d62728", ls=ls, lw=1.8, label="отражение R, среда над столбиками n = %.2f" % ns)
            a = arr[:, 3] if d == "improved" else arr[:, 2]
            ax.plot(arr[:, 0], a, color="#1f77b4", ls=ls, lw=1.8, label="поглощение, n = %.2f" % ns)
        ax.set_title(LAB[d], fontsize=10)
        ax.set_xlabel("длина волны в свободном пространстве, нм")
        ax.set_ylabel("доля мощности")
        ax.set_ylim(0, 1.02)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, loc="center left")
    fig.suptitle("Поглотитель как рефрактометрический датчик: воздух (сплошные), вода 1,33 (штрих) и раствор 1,35 (пунктир) "
                 "над столбиками и в зазорах", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(L.RESULTS / "ri_sensor.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("solar", "all"):
        plot_solar(read_solar())
    if which in ("ri", "all"):
        plot_ri()
    print("рисунки записаны")
