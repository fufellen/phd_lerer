"""Самосогласованный расчёт нелинейной метаповерхности с золотыми наночастицами (среднеполевой уровень).

Закрывает разрыв, названный в заметке «Самосогласованный расчёт нелинейной метаповерхности»: во
флагманской фигуре нелинейной статьи (grating_2d_rcwa/scripts/rcwa_nonlinear_au_absorption.py)
нелинейная проницаемость композита считалась по ПАДАЮЩЕЙ интенсивности, как будто композит видит
плоскую волну в матрице с той же интенсивностью, а затем один раз решалась линейная дифракция
(расцепленное приближение). Здесь поле, сформированное самой решёткой, возвращается в модель материала:

1. решается дифракция при текущей eps_eff цилиндров и слоя (grcwa, двумерная решётка цилиндров);
2. по поглощению каждого слоя A_j = k0 Im eps_j <|E|^2>_j t_j f_j / |E0|^2 восстанавливается среднее
   |E|^2 в композите цилиндров и слоя (потоковый баланс, без объёмного интеграла);
3. по <|E|^2> пересчитывается проницаемость частицы (ветка A, феноменологическая chi^(3), локальное поле
   квазистатической сферы) и эффективная проницаемость композита по Максвеллу-Гарнетту;
4. обновление с демпфированием, повтор до сходимости eps, R, T, P.

Геометрия статьи Лерер АМ_1: период 500 нм, цилиндры D = 400 нм высотой 100 нм в воздухе, слой 100 нм,
подложка n = 1,45; композит Au 10 % (Джонсон-Кристи) в матрице 1,77; chi^(3) - benchmark Hache 1988
(пс, у резонанса), как во флагманской фигуре. Интенсивности 1e6, 5e6, 2e7 Вт/см^2 - те же.

Запуск: python nonlinear_selfconsistent.py [conv|main|hyst|all]
Выход: results/nonlinear_selfconsistent*.csv|txt
"""
from __future__ import annotations

import io
import math
import sys
import time

import numpy as np

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import lerer_grcwa as L  # noqa: E402

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
import nonlinear_au_sweep as nl  # noqa: E402  (composite_ema/scripts уже в sys.path через lerer_grcwa)

PERIOD, DIAM, H_CYL, H_FILM = 500.0, 400.0, 100.0, 100.0
EPS_H = complex(L.N_HOST ** 2, 0.0)
CHI_OPT = nl.BENCHMARKS[0].chi_opt_si          # Hache 1988, оптическая конвенция
INTENSITIES = [1.0e6, 5.0e6, 2.0e7]            # Вт/см^2
LAMS = list(range(400, 701, 10))
NG = 121
NX = 120
FILL_CYL = None   # доля площади цилиндра на сетке решателя, вычисляется в setup()
ETA = 0.5
TOL = 1e-6
MAX_IT = 80


def e0sq_from_intensity(i_wcm2: float) -> float:
    """|E0|^2 падающей волны в воздухе, В^2/м^2, из интенсивности в Вт/см^2."""
    return 2.0 * i_wcm2 * 1.0e4 / (L.EPS0 * L.C_LIGHT)


def eps_p_nonlinear(lam_nm: float, e2: float, chi_opt: complex = CHI_OPT):
    """Проницаемость частицы Au при среднем квадрате макроскопического поля e2 в композите
    (оптическая конвенция): eps_p = eps_L + (3/4) chi |f|^2 e2, f = 3 eps_h / (eps_p + 2 eps_h)."""
    eps_l = L.eps_au_jc(lam_nm)
    if e2 == 0.0:
        return eps_l, True
    eps_p = eps_l
    for _ in range(400):
        f = 3.0 * EPS_H / (eps_p + 2.0 * EPS_H)
        new = eps_l + 0.75 * chi_opt * abs(f) ** 2 * e2
        if abs(new - eps_p) < 1e-12 * max(1.0, abs(new)):
            return new, True
        eps_p = 0.5 * eps_p + 0.5 * new
    return eps_p, False


def eps_eff_decoupled(lam_nm: float, i_wcm2: float) -> complex:
    """Расцепленное приближение флагманской фигуры: |E|^2 = 2I/(n_h eps0 c) - плоская волна в матрице."""
    e2 = 2.0 * i_wcm2 * 1.0e4 / (L.N_HOST * L.EPS0 * L.C_LIGHT)
    eps_p, _ = eps_p_nonlinear(lam_nm, e2)
    return L.mg(EPS_H, eps_p, L.C_FILL)


def self_test():
    """Согласие с nonlinear_au_sweep (конвенция Лерера) и линейный предел."""
    for lam in (450.0, 532.0, 650.0):
        for i in (1.0e6, 2.0e7):
            e2 = 2.0 * i * 1.0e4 / (L.N_HOST * L.EPS0 * L.C_LIGHT)
            mine, _ = eps_p_nonlinear(lam, e2)
            ref, _ = nl.eps_p_nonlinear(lam, i, nl.BENCHMARKS[0].chi_lerer_si, EPS_H)
            assert abs(mine - ref.conjugate()) < 1e-9 * abs(ref), (lam, i, mine, ref)
    assert eps_p_nonlinear(532.0, 0.0)[0] == L.eps_au_jc(532.0)
    print("[self_test] eps_p(I) совпадает с nonlinear_au_sweep во всех точках; линейный предел точен")


def solve_at(lam, eps_c, eps_f, nG=NG):
    layers, eps_sub = L.stack_cylinder_metasurface(lam, eps_c, eps_f, PERIOD, DIAM, H_CYL, H_FILM)
    s = L.solve(lam, PERIOD / 1000.0, layers, eps_sub=eps_sub, want_layers=True, nG=nG, nx=NX)
    return s


def setup():
    global FILL_CYL
    p = L.Pillar(H_CYL / 1000.0, DIAM / 1000.0, DIAM / 1000.0, 1.0 + 0j, 1.0 + 0j, "circle")
    FILL_CYL = L.fill_fraction(p, PERIOD / 1000.0, NX)
    print("[setup] доля площади цилиндра на сетке %d x %d: %.5f (точно pi D^2/4P^2 = %.5f)"
          % (NX, NX, FILL_CYL, math.pi * DIAM ** 2 / 4 / PERIOD ** 2))


def field_ratios(lam, sol, eps_c, eps_f):
    """<|E|^2>/|E0|^2 в композите цилиндров и слоя."""
    r_c = L.mean_e2_ratio(sol.A_layers[0], lam, eps_c.imag, H_CYL / 1000.0, FILL_CYL)
    r_f = L.mean_e2_ratio(sol.A_layers[1], lam, eps_f.imag, H_FILM / 1000.0, 1.0)
    return r_c, r_f


def solve_selfconsistent(lam, i_wcm2, start=None, eta=ETA, tol=TOL, max_it=MAX_IT, log=None):
    """Возвращает (sol, eps_c, eps_f, r_c, r_f, iterations, converged)."""
    e0sq = e0sq_from_intensity(i_wcm2)
    eps_lin = L.mg(EPS_H, L.eps_au_jc(lam), L.C_FILL)
    eps_c, eps_f = (eps_lin, eps_lin) if start is None else start
    conv = False
    it = 0
    for it in range(1, max_it + 1):
        sol = solve_at(lam, eps_c, eps_f)
        r_c, r_f = field_ratios(lam, sol, eps_c, eps_f)
        ep_c, ok1 = eps_p_nonlinear(lam, r_c * e0sq)
        ep_f, ok2 = eps_p_nonlinear(lam, r_f * e0sq)
        new_c, new_f = L.mg(EPS_H, ep_c, L.C_FILL), L.mg(EPS_H, ep_f, L.C_FILL)
        d = max(abs(new_c - eps_c) / max(1.0, abs(new_c)), abs(new_f - eps_f) / max(1.0, abs(new_f)))
        if log is not None:
            log.append((it, d, sol.A))
        eps_c = (1.0 - eta) * eps_c + eta * new_c
        eps_f = (1.0 - eta) * eps_f + eta * new_f
        if d < tol and ok1 and ok2:
            conv = True
            break
    sol = solve_at(lam, eps_c, eps_f)
    r_c, r_f = field_ratios(lam, sol, eps_c, eps_f)
    return sol, eps_c, eps_f, r_c, r_f, it, conv


def run_conv():
    lines = ["сходимость по гармоникам, линейный композит, A / <|E|^2>n_h/|E0|^2 (цилиндры, слой)"]
    for nG in (61, 121, 221):
        t0 = time.time()
        row = []
        for lam in (450.0, 550.0, 650.0):
            eps = L.mg(EPS_H, L.eps_au_jc(lam), L.C_FILL)
            s = solve_at(lam, eps, eps, nG)
            r_c, r_f = field_ratios(lam, s, eps, eps)
            row.append("%d нм: A=%.4f, %.3f/%.3f" % (lam, s.A, r_c * L.N_HOST, r_f * L.N_HOST))
        lines.append("nG = %3d (%.0f с): %s" % (nG, time.time() - t0, " | ".join(row)))
        print(lines[-1])
    io.open(L.RESULTS / "nonlinear_selfconsistent_convergence.txt", "w", encoding="utf-8").write("\n".join(lines) + "\n")


def run_main():
    out = io.open(L.RESULTS / "nonlinear_selfconsistent.csv", "w", encoding="utf-8", newline="\n")
    out.write("lambda_nm,I_wcm2,P_lin,P_dec,P_sc,dP_dec,dP_sc,ratio_cyl_lin,ratio_film_lin,ratio_cyl_sc,ratio_film_sc,"
              "Re_eps_cyl_sc,Im_eps_cyl_sc,Re_eps_film_sc,Im_eps_film_sc,Re_eps_dec,Im_eps_dec,iterations,converged\n")
    t0 = time.time()
    for lam in LAMS:
        eps_lin = L.mg(EPS_H, L.eps_au_jc(lam), L.C_FILL)
        s_lin = solve_at(lam, eps_lin, eps_lin)
        r_c0, r_f0 = field_ratios(lam, s_lin, eps_lin, eps_lin)
        start = None
        msg = []
        for i_wcm2 in INTENSITIES:
            eps_d = eps_eff_decoupled(lam, i_wcm2)
            s_dec = solve_at(lam, eps_d, eps_d)
            s_sc, eps_c, eps_f, r_c, r_f, it, conv = solve_selfconsistent(lam, i_wcm2, start)
            start = (eps_c, eps_f)
            out.write("%d,%.3e,%.6f,%.6f,%.6f,%+.6f,%+.6f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%d,%d\n"
                      % (lam, i_wcm2, s_lin.A, s_dec.A, s_sc.A, s_dec.A - s_lin.A, s_sc.A - s_lin.A,
                         r_c0 * L.N_HOST, r_f0 * L.N_HOST, r_c * L.N_HOST, r_f * L.N_HOST,
                         eps_c.real, eps_c.imag, eps_f.real, eps_f.imag, eps_d.real, eps_d.imag, it, int(conv)))
            msg.append("I=%.0e: dP dec %+.4f, sc %+.4f (%d ит.%s)" % (i_wcm2, s_dec.A - s_lin.A, s_sc.A - s_lin.A, it,
                                                                   "" if conv else " НЕ СОШЛОСЬ"))
        out.flush()
        print("%d нм: P_lin=%.4f, n_h<|E|^2>/|E0|^2 цил %.3f слой %.3f | %s  [%.0f с]"
              % (lam, s_lin.A, r_c0 * L.N_HOST, r_f0 * L.N_HOST, "; ".join(msg), time.time() - t0))
    out.close()


def run_hysteresis():
    """Прямой и обратный проходы по интенсивности при трёх длинах волн плюс проверка более сильных полей."""
    lines = ["прямой (от малой I) и обратный (от большой I) проходы; P_sc и число итераций"]
    for lam in (530.0, 550.0, 600.0):
        fwd, start = [], None
        for i_wcm2 in INTENSITIES:
            s, eps_c, eps_f, r_c, r_f, it, conv = solve_selfconsistent(lam, i_wcm2, start)
            start = (eps_c, eps_f)
            fwd.append((i_wcm2, s.A, it, conv))
        bwd, start = [], start
        for i_wcm2 in reversed(INTENSITIES):
            s, eps_c, eps_f, r_c, r_f, it, conv = solve_selfconsistent(lam, i_wcm2, start)
            start = (eps_c, eps_f)
            bwd.append((i_wcm2, s.A, it, conv))
        bwd = list(reversed(bwd))
        for (i1, p1, it1, c1), (i2, p2, it2, c2) in zip(fwd, bwd):
            lines.append("%d нм I=%.0e: прямой P=%.6f (%d ит.), обратный P=%.6f (%d ит.), разница %.1e%s"
                         % (lam, i1, p1, it1, p2, it2, abs(p1 - p2), "" if (c1 and c2) else "  НЕ СОШЛОСЬ"))
            print(lines[-1])
    lines.append("проверка более сильных полей при 550 нм (за пределами возмущательной ветки benchmark)")
    start = None
    for i_wcm2 in (1.0e8, 5.0e8):
        log = []
        s, eps_c, eps_f, r_c, r_f, it, conv = solve_selfconsistent(550.0, i_wcm2, start, log=log)
        start = (eps_c, eps_f)
        lines.append("I=%.0e: P_sc=%.4f, %d ит., сошлось=%s, eps_cyl=%.3f%+.3fi, eps_film=%.3f%+.3fi, "
                     "последние невязки %s" % (i_wcm2, s.A, it, conv, eps_c.real, eps_c.imag, eps_f.real, eps_f.imag,
                                                " ".join("%.1e" % d for _, d, _ in log[-4:])))
        print(lines[-1])
    io.open(L.RESULTS / "nonlinear_selfconsistent_hysteresis.txt", "w", encoding="utf-8").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    t0 = time.time()
    self_test()
    setup()
    if which in ("conv", "all"):
        run_conv()
    if which in ("main", "all"):
        run_main()
    if which in ("hyst", "all"):
        run_hysteresis()
    print("готово за %.0f с" % (time.time() - t0))
