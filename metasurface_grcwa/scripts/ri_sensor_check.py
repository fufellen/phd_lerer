"""Проверка поглотителя как рефрактометрического датчика: что меняется, когда воздух над столбиками
и в зазорах между ними сменяется водой (n = 1,33) и водным раствором (1,34; 1,35).

Приложение «сенсоры и детекторы» из введения статьи Лерер АМ_1 и первый приоритет контрольной точки
«цели формула КМ поглотители нелинейность Au». Меряется сдвиг минимума отражения и изменение
поглощения на единицу показателя преломления (нм/RIU) для исходной и предложенной конструкций.

Запуск: python ri_sensor_check.py
Выход: results/ri_sensor_spectra.csv, results/ri_sensor_summary.txt
"""
from __future__ import annotations

import io
import sys
import time

import numpy as np

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import lerer_grcwa as L  # noqa: E402

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

AU = L.AuModel("lerer")
LAMS = list(range(400, 751, 5))
SUPER = [1.00, 1.33, 1.34, 1.35]
DESIGNS = {"original": (L.ORIGINAL, 221), "improved": (L.IMPROVED, 121)}


def main():
    t0 = time.time()
    out = io.open(L.RESULTS / "ri_sensor_spectra.csv", "w", encoding="utf-8", newline="\n")
    out.write("design,n_super,lambda_nm,R,T,A_layers,A_total\n")
    lines = []
    for design, (d, nG) in DESIGNS.items():
        res = {}
        for ns in SUPER:
            rows = []
            for lam in LAMS:
                layers, eps_sub = L.stack_pillar_absorber(float(lam), AU, eps_super=complex(ns ** 2, 0.0), **d)
                s = L.solve(float(lam), d["period_nm"] / 1000.0, layers, eps_sub=eps_sub, nG=nG,
                            eps_super=complex(ns ** 2, 0.0))
                rows.append((lam, s.R, s.T, s.A, 1.0 - s.R))
                out.write("%s,%.2f,%d,%.6f,%.6f,%.6f,%.6f\n" % (design, ns, lam, s.R, s.T, s.A, 1.0 - s.R))
            out.flush()
            res[ns] = np.array(rows)
            print("%s n=%.2f: %.0f с" % (design, ns, time.time() - t0))
        lines.append("%s" % design)
        for ns in SUPER:
            r = res[ns]
            a = r[:, 4] if d["mirror"] else r[:, 3]
            i = int(np.argmin(r[:, 1]))
            lines.append("  n_super = %.2f: среднее A(400-750) = %.4f, min A = %.4f при %d нм, "
                         "минимум R = %.4f при %d нм, R(650) = %.4f"
                         % (ns, a.mean(), a.min(), int(r[np.argmin(a), 0]), r[i, 1], int(r[i, 0]), r[LAMS.index(650), 1]))
        # чувствительность по положению минимума отражения и по уровню отражения
        lam_min = {ns: int(res[ns][np.argmin(res[ns][:, 1]), 0]) for ns in SUPER}
        s_pos = (lam_min[1.35] - lam_min[1.33]) / 0.02
        dR = {lam: (res[1.35][LAMS.index(lam), 1] - res[1.33][LAMS.index(lam), 1]) / 0.02 for lam in (450, 550, 650, 750)}
        lines.append("  сдвиг минимума R между n = 1,33 и 1,35: %d -> %d нм, чувствительность %.0f нм/RIU"
                     % (lam_min[1.33], lam_min[1.35], s_pos))
        lines.append("  dR/dn при 450/550/650/750 нм: %s (доля/RIU)" % " ".join("%+.3f" % dR[l] for l in (450, 550, 650, 750)))
        lines.append("  ширина минимума R при n = 1,33 (уровень удвоенного минимума): %s"
                     % width_at(res[1.33][:, 0], res[1.33][:, 1]))
    out.close()
    txt = "\n".join(lines) + "\n"
    print(txt)
    io.open(L.RESULTS / "ri_sensor_summary.txt", "w", encoding="utf-8").write(txt)


def width_at(lam, R):
    i = int(np.argmin(R))
    level = 2.0 * R[i]
    lo = lam[i]
    for j in range(i, -1, -1):
        if R[j] > level:
            break
        lo = lam[j]
    hi = lam[i]
    for j in range(i, len(lam)):
        if R[j] > level:
            break
        hi = lam[j]
    return "%d...%d нм (%d нм)" % (lo, hi, hi - lo)


if __name__ == "__main__":
    main()
