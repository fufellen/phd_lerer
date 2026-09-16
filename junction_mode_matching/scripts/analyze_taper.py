"""Сведение переходов по ширине золота: трёхмерная модель COMSOL (RunTaper.java) против EME по местным модам
(taper_eme.py) и адиабатического эталона exp(-int alpha dz).

Модовая передача из S-параметров, как в analyze_reverse.py: T = |S|^2 / (exp(-alpha_1 l_1) exp(-alpha_2 l_2)),
только теперь на участке от стыка до порта 2 металл имеет полную ширину не всюду - опорная плоскость порта 2
стоит при z = l_pcm, а переход занимает z от 0 до L. Делить на затухание рабочей моды по всей длине l_pcm
означает отнести потери самого перехода (его лишнее поглощение и преобразование) к "потерям перехода" -
это та же мера, что и в заметке о трёх соединениях: сколько дошло бы, будь переход идеальным и без металла
сверх рабочего, минус то, что дошло.

Запуск: python analyze_taper.py <папка с taper_results.csv> <папка вывода> [taper_eme_results.csv]
"""
import csv
import io
import math
import os
import sys

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

K0 = 2 * math.pi / 1.55
RU = {"linear": "линейный", "quad": "парабола, острый", "sqrt": "парабола, тупой", "ellipse": "эллипс, тупой",
      "ellipsec": "эллипс, острый", "expo": "экспонента", "gauss": "гаусс", "rcos": "приподнятый косинус",
      "klop": "Клопфенштейн", "fastwin": "окно за 5 %", "slowwin": "окно за 50 %", "step": "без перехода (скачок при z = L)",
      "cubic": "куб", "quart": "корень четвёртой степени", "table": "подобранный"}


def main():
    src = sys.argv[1]
    dst = sys.argv[2] if len(sys.argv) > 2 else src
    eme_path = sys.argv[3] if len(sys.argv) > 3 else None
    rows = list(csv.DictReader(io.open(os.path.join(src, "taper_results.csv"), encoding="utf-8")))
    cases = {}
    for r in rows:
        cases[r["case"]] = r
    res = {}
    for name, r in sorted(cases.items()):
        n1 = complex(float(r["neff1_re"]), float(r["neff1_im"]))
        n2 = complex(float(r["neff2_re"]), float(r["neff2_im"]))
        a1, a2 = 2 * K0 * abs(n1.imag), 2 * K0 * abs(n2.imag)
        lf, lp = float(r["l_feed_um"]), float(r["l_pcm_um"])
        decay = math.exp(-a1 * lf) * math.exp(-a2 * lp)
        if r["dir"] == "fwd":
            s_t = complex(float(r["S21_re"]), float(r["S21_im"]))
            s_r = complex(float(r["S11_re"]), float(r["S11_im"]))
        else:
            s_t = complex(float(r["S12_re"]), float(r["S12_im"]))
            s_r = complex(float(r["S22_re"]), float(r["S22_im"]))
        T = abs(s_t) ** 2 / decay
        key = (r["shape"], float(r["taper_um"]), r["dir"], r["mesh"])
        res[key] = dict(T=T, R=abs(s_r) ** 2, a2=a2, p28=float(r["p_pcm_end"]), elements=r["elements"], case=name)
        print("%-26s T = %.4f (%.3f дБ)  |S_refl|^2 = %.2e  P(2.8) = %.4f  el = %s" % (name, T, -10 * math.log10(T), abs(s_r) ** 2, float(r["p_pcm_end"]), r["elements"]))
    eme = {}
    if eme_path and os.path.exists(eme_path):
        for r in csv.DictReader(io.open(eme_path, encoding="utf-8")):
            eme[(r["shape"], float(r["L_um"]))] = r
    out = io.open(os.path.join(dst, "taper_summary.csv"), "w", encoding="utf-8", newline="\n")
    out.write("shape,shape_ru,L_um,mesh,T_fwd_comsol,T_rev_comsol,ratio_comsol,loss_fwd_dB,loss_rev_dB,T_fwd_eme,T_rev_eme,T_ad_eme,loss_fwd_eme_dB,loss_ad_dB\n")
    print("\nсводка по формам (COMSOL вперёд / назад, EME вперёд, адиабатический эталон):")
    table = []
    for (shape, L, d, mesh), s in sorted(res.items()):
        if d != "fwd":
            continue
        rev = res.get((shape, L, "rev", mesh))
        e = eme.get((shape, L))
        Tf, Tr = s["T"], (rev["T"] if rev else float("nan"))
        Te = float(e["T_fwd_eq"]) if e else float("nan")
        Ter = float(e["T_rev_eq"]) if e else float("nan")
        Ta = float(e["T_ad_eq"]) if e else float("nan")
        print("  %-24s L=%.2f: COMSOL %.4f / %.4f (%.3f / %.3f дБ, отн. %.4f)   EME %.4f / %.4f   адиаб. %.4f"
              % (RU.get(shape, shape), L, Tf, Tr, -10 * math.log10(Tf), -10 * math.log10(Tr) if Tr == Tr else float("nan"),
                 Tr / Tf if Tr == Tr else float("nan"), Te, Ter, Ta))
        out.write("%s,%s,%.3f,%s,%.6f,%.6f,%.6f,%.4f,%.4f,%.6f,%.6f,%.6f,%.4f,%.4f\n"
                  % (shape, RU.get(shape, shape), L, mesh, Tf, Tr, Tr / Tf if Tr == Tr else float("nan"),
                     -10 * math.log10(Tf), -10 * math.log10(Tr) if Tr == Tr else float("nan"), Te, Ter, Ta,
                     -10 * math.log10(Te) if Te == Te else float("nan"), -10 * math.log10(Ta) if Ta == Ta else float("nan")))
        table.append((shape, L, Tf, Tr, Te, Ter, Ta))
    out.close()

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    # рисунок 1: профили острия
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from taper_eme import profile
    u = np.linspace(0, 1, 400)
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    for sh, ls in (("linear", "-"), ("quad", "-"), ("sqrt", "--"), ("ellipse", "--"), ("ellipsec", "-"), ("expo", "-."), ("klop", ":")):
        ax.plot(u, 450 * profile(sh, u) / 2, ls, label=RU[sh])
    ax.axhspan(105, 145, color="gray", alpha=0.25, label="окно перестройки 210…290 нм (полуширина 105…145)")
    ax.set_xlabel("z / L вдоль перехода")
    ax.set_ylabel("полуширина золота, нм")
    ax.set_title("Профили острия металлического перехода (вид сверху, половина по симметрии)", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0.0)
    fig.tight_layout()
    fig.savefig(os.path.join(dst, "taper_profiles.png"), dpi=150, bbox_inches="tight")
    # рисунок 2: потери по формам для каждой длины
    Ls = sorted(set(t[1] for t in table))
    fig, axes = plt.subplots(1, len(Ls), figsize=(6.5 * len(Ls), 4.6), squeeze=False)
    for ax, L in zip(axes[0], Ls):
        rows_L = [t for t in table if t[1] == L]
        rows_L.sort(key=lambda t: -t[2])
        names = [RU.get(t[0], t[0]) for t in rows_L]
        x = np.arange(len(rows_L))
        wdt = 0.2
        ax.bar(x - 1.5 * wdt, [-10 * math.log10(t[2]) for t in rows_L], wdt, label="COMSOL, вперёд", color="#1f77b4")
        ax.bar(x - 0.5 * wdt, [-10 * math.log10(t[3]) if t[3] == t[3] else 0 for t in rows_L], wdt, label="COMSOL, назад", color="#d62728")
        ax.bar(x + 0.5 * wdt, [-10 * math.log10(t[4]) if t[4] == t[4] else 0 for t in rows_L], wdt, label="EME по местным модам, вперёд", color="#2ca02c")
        ax.bar(x + 1.5 * wdt, [-10 * math.log10(t[6]) if t[6] == t[6] else 0 for t in rows_L], wdt, label="адиабатический эталон exp(−∫α dz)", color="#9467bd")
        ax.set_xticks(x)
        ax.set_xticklabels(names, rotation=30, ha="right", fontsize=8)
        ax.set_ylabel("потери перехода, дБ")
        ax.set_title("длина перехода L = %.1f мкм" % L, fontsize=10)
        ax.grid(alpha=0.3, axis="y")
        ax.legend(fontsize=7)
    fig.suptitle("Потери перехода по ширине золота против формы острия: сколько не дошло в рабочую моду относительно перехода без металла", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(dst, "taper_losses.png"), dpi=150)
    print("рисунки: taper_profiles.png, taper_losses.png")


if __name__ == "__main__":
    main()
