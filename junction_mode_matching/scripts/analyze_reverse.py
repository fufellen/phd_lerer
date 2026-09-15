"""Обработка прямых и обратных прогонов COMSOL (RunReverse.java): модовая передача стыка из S-параметров.

Порт 1 стоит в подводящем волноводе на z = -l_feed, порт 2 - в фазовращателе на z = +l_pcm. S-параметр
численного порта - амплитуда в моде порта, отнесённая к падающей, с опорными плоскостями на портах, поэтому
    |S21|^2 = T_fwd exp(-alpha_1 l_feed) exp(-alpha_2 l_pcm),   |S12|^2 = T_rev exp(-alpha_2 l_pcm) exp(-alpha_1 l_feed),
где alpha_i = 2 k0 |Im n_i| - затухание по мощности моды порта i, взятое из самой модели (столбцы neff1, neff2).
T_fwd - доля мощности, вошедшая из моды подвода в плазмонную моду; T_rev - из плазмонной в моду подвода.

Запуск: python analyze_reverse.py <папка с reverse_results.csv и reverse_profiles.csv> [выходная папка]
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


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "."
    dst = sys.argv[2] if len(sys.argv) > 2 else src
    rows = list(csv.DictReader(io.open(os.path.join(src, "reverse_results.csv"), encoding="utf-8")))
    # последняя строка на случай выигрывает (повторные прогоны дописываются)
    cases = {}
    for r in rows:
        cases[r["case"]] = r
    out = io.open(os.path.join(dst, "reverse_summary.csv"), "w", encoding="utf-8", newline="\n")
    out.write("case,scheme,dir,state,mesh,l_feed,l_pcm,elements,neff1_re,neff1_im,neff2_re,neff2_im,alpha1_1um,alpha2_1um,"
              "S_refl_abs2,S_trans_abs2,T_joint,loss_dB,R_joint_dB\n")
    summary = {}
    for name, r in sorted(cases.items()):
        n1 = complex(float(r["neff1_re"]), float(r["neff1_im"]))
        n2 = complex(float(r["neff2_re"]), float(r["neff2_im"]))
        a1 = 2 * K0 * abs(n1.imag)
        a2 = 2 * K0 * abs(n2.imag)
        lf, lp = float(r["l_feed_um"]), float(r["l_pcm_um"])
        decay = math.exp(-a1 * lf) * math.exp(-a2 * lp)
        if r["dir"] == "fwd":
            s_t = complex(float(r["S21_re"]), float(r["S21_im"]))
            s_r = complex(float(r["S11_re"]), float(r["S11_im"]))
        else:
            s_t = complex(float(r["S12_re"]), float(r["S12_im"]))
            s_r = complex(float(r["S22_re"]), float(r["S22_im"]))
        T = abs(s_t) ** 2 / decay
        R = abs(s_r) ** 2
        # отражение измерено на порту: у обратного прогона оно прошло путь 2 l_pcm в затухающей моде
        R_joint = R / math.exp(-2 * (a2 if r["dir"] == "rev" else a1) * (lp if r["dir"] == "rev" else lf))
        state = "c" if float(r["n_pcm"]) > 3.0 else "a"
        out.write("%s,%s,%s,%s,%s,%.3f,%.3f,%s,%.7f,%.4e,%.7f,%.4e,%.5f,%.5f,%.4e,%.6f,%.5f,%.4f,%.1f\n"
                  % (name, r["scheme"], r["dir"], state, r["mesh"], lf, lp, r["elements"], n1.real, n1.imag, n2.real, n2.imag,
                     a1, a2, R, abs(s_t) ** 2, T, -10 * math.log10(T), 10 * math.log10(max(R_joint, 1e-30))))
        summary[name] = dict(T=T, R=R_joint, scheme=r["scheme"], dir=r["dir"], state=state, mesh=r["mesh"], lf=lf, a1=a1, a2=a2)
        print("%-24s T = %.5f (%.3f дБ)   R_joint = %.2e   alpha1 = %.5f alpha2 = %.5f  |S_t|^2 = %.5f"
              % (name, T, -10 * math.log10(T), R_joint, a1, a2, abs(s_t) ** 2))
    out.close()
    # пары вперёд/назад одной геометрии
    print("\nпары вперёд / назад (одна геометрия, одна сетка):")
    with io.open(os.path.join(dst, "reverse_pairs.csv"), "w", encoding="utf-8", newline="\n") as f:
        f.write("scheme,state,mesh,l_feed,T_fwd,T_rev,ratio_rev_fwd,loss_fwd_dB,loss_rev_dB\n")
        for name, s in sorted(summary.items()):
            if s["dir"] != "fwd":
                continue
            partner = name.replace("_fwd_", "_rev_")
            if partner in summary:
                p = summary[partner]
                print("  %-10s %s сетка %s l_feed %.1f: T_fwd = %.5f, T_rev = %.5f, T_rev/T_fwd = %.5f"
                      % (s["scheme"], s["state"], s["mesh"], s["lf"], s["T"], p["T"], p["T"] / s["T"]))
                f.write("%s,%s,%s,%.3f,%.6f,%.6f,%.6f,%.4f,%.4f\n" % (s["scheme"], s["state"], s["mesh"], s["lf"], s["T"], p["T"],
                        p["T"] / s["T"], -10 * math.log10(s["T"]), -10 * math.log10(p["T"])))

    # профили
    prof = {}
    for r in csv.DictReader(io.open(os.path.join(src, "reverse_profiles.csv"), encoding="utf-8")):
        prof.setdefault(r["case"], []).append((float(r["z_um"]), float(r["p_z"])))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.6))
    for ax, key, title in ((axes[0], "sync405_fwd_a", "вперёд: мода подводящего волновода с порта 1 (z = −0,9 мкм)"),
                           (axes[1], "sync405_rev_a", "назад: плазмонная мода с порта 2 (z = +3 мкм)"),
                           (axes[2], "sync405_rev_a_lf250", "назад, подводящий волновод удлинён до 2,5 мкм")):
        if key not in prof:
            continue
        z, p = np.array(sorted(prof[key])).T
        p = np.abs(p)
        ax.plot(z, p, "o-", color="#1f77b4", label="поток мощности через плоскость, трёхмерная модель")
        s = summary.get(key)
        if s is not None:
            zz = np.linspace(min(z), max(z), 300)
            lf = s["lf"]
            if s["dir"] == "fwd":
                ref = np.where(zz < 0, np.exp(-s["a1"] * (zz + lf)), s["T"] * np.exp(-s["a1"] * lf) * np.exp(-s["a2"] * zz))
                lab = "мода: exp(−α z), доля в моде T = %.3f по S21" % s["T"]
            else:
                lp = 3.0
                ref = np.where(zz > 0, np.exp(-s["a2"] * (lp - zz)), s["T"] * np.exp(-s["a2"] * lp) * np.exp(-s["a1"] * (-zz)))
                lab = "мода: exp(−α z), доля в моде T = %.3f по S12" % s["T"]
            ax.plot(zz, ref, "--", color="#d62728", label=lab)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_xlabel("z вдоль волновода, мкм (стык при z = 0)")
        ax.set_ylabel("доля запущенной мощности")
        ax.set_title(title, fontsize=10)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc="lower left" if s and s["dir"] == "fwd" else "upper left")
    fig.suptitle("Стык sync405, аморфное состояние: продольные профили мощности в двух направлениях (трёхмерная модель COMSOL, базовая сетка)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(dst, "reverse_profiles.png"), dpi=150)
    print("рисунок:", os.path.join(dst, "reverse_profiles.png"))


if __name__ == "__main__":
    main()
