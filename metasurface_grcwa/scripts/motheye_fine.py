"""Сетка 10 нм для конусов «глаз мотылька» предложенной конструкции: есть ли у них провал между точками серии.
Запуск: python motheye_fine.py -> results/motheye_fine.csv, results/motheye_fine_summary.txt"""
import io, sys, time
import numpy as np
sys.path.insert(0, ".")
import lerer_grcwa as L
import pillar_shapes as PS
import pillar_shapes_extra as PX
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
LAMS = list(range(400, 751, 10))
out = io.open(L.RESULTS / "motheye_fine.csv", "w", encoding="utf-8", newline="\n")
out.write("design,variant,lambda_nm,R,T,A_layers,A_total\n")
lines = []
t0 = time.time()
for design in ("improved", "original"):
    d = PS.DESIGNS[design]
    for variant in ("motheye_h300", "motheye_h500"):
        pl = PX.extra_layers(variant, design)
        A = []
        for lam in LAMS:
            layers, eps_sub = L.stack_pillar_absorber(float(lam), PS.AU, pillar_layers=pl, **{k: v for k, v in d.items() if k not in ("pillar_nm", "h_pil_nm")})
            s = L.solve(float(lam), d["period_nm"] / 1000.0, layers, eps_sub=eps_sub, nG=PS.NG[design])
            A.append(s.A)
            out.write("%s,%s,%d,%.6f,%.6f,%.6f,%.6f\n" % (design, variant, lam, s.R, s.T, s.A, 1 - s.R))
        out.flush()
        A = np.array(A); j = int(np.argmin(A))
        lines.append("%s %s: среднее по сетке 10 нм %.4f, наименьшее %.4f при %d нм (%.0f с)" % (design, variant, A.mean(), A[j], LAMS[j], time.time() - t0))
        print(lines[-1])
out.close()
io.open(L.RESULTS / "motheye_fine_summary.txt", "w", encoding="utf-8").write("\n".join(lines) + "\n")
