"""
Валидация tmm.py: сравнение с независимыми (учебными, не из статей Лерера)
формулами Френеля и Эйри, плюс первый (нестрогий) прогон через
composite_ema, чтобы увидеть, ведёт ли себя резонанс Фрёлиха ожидаемо уже
на уровне однородного (непериодического) слоя - до подключения настоящей
2D-периодической дифракции.
"""

from __future__ import annotations

import cmath
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "composite_ema" / "scripts"))

from tmm import Layer, power_rt  # noqa: E402
import composite_ema as ema  # noqa: E402


def fresnel_normal(n1: complex, n2: complex) -> tuple[float, float]:
    """Формулы Френеля при нормальном падении (учебник), R и T по мощности."""
    r = (n1 - n2) / (n1 + n2)
    R = abs(r) ** 2
    T = 1.0 - R  # при нормальном падении и вещественных n2/n1 (беспотерьные среды)
    return R, T


def fresnel_oblique(n1: float, n2: float, theta_i_deg: float, pol: str) -> tuple[float, float]:
    """Формулы Френеля под углом (учебник, беспотерьные вещественные среды)."""
    theta_i = math.radians(theta_i_deg)
    sin_t = n1 * math.sin(theta_i) / n2
    theta_t = math.asin(sin_t)
    if pol == "s":
        r = (n1 * math.cos(theta_i) - n2 * math.cos(theta_t)) / (n1 * math.cos(theta_i) + n2 * math.cos(theta_t))
    else:
        r = (n2 * math.cos(theta_i) - n1 * math.cos(theta_t)) / (n2 * math.cos(theta_i) + n1 * math.cos(theta_t))
    R = r ** 2
    T = 1.0 - R
    return R, T


def airy_thin_film(n0: complex, n1: complex, n2: complex, d_nm: float, lam_nm: float) -> tuple[float, float]:
    """Классическая формула Эйри для тонкой плёнки n1 толщиной d между n0 и n2, нормальное падение."""
    r01 = (n0 - n1) / (n0 + n1)
    r12 = (n1 - n2) / (n1 + n2)
    t01 = 2 * n0 / (n0 + n1)
    t12 = 2 * n1 / (n1 + n2)
    beta = 2 * cmath.pi * n1 * d_nm / lam_nm
    r = (r01 + r12 * cmath.exp(2j * beta)) / (1 + r01 * r12 * cmath.exp(2j * beta))
    t = (t01 * t12 * cmath.exp(1j * beta)) / (1 + r01 * r12 * cmath.exp(2j * beta))
    R = abs(r) ** 2
    T = (abs(t) ** 2) * (n2 / n0).real
    return R, T


def check(label: str, got: tuple[float, float], expected: tuple[float, float], tol: float = 1e-9) -> None:
    dR = abs(got[0] - expected[0])
    dT = abs(got[1] - expected[1])
    status = "OK" if (dR < tol and dT < tol) else "FAIL"
    print(f"[{status}] {label}: tmm=({got[0]:.6f},{got[1]:.6f})  expected=({expected[0]:.6f},{expected[1]:.6f})  "
          f"diff=({dR:.2e},{dT:.2e})")
    assert dR < tol and dT < tol, label


def main() -> None:
    lam = 600.0

    # 1) Одна граница раздела, нормальное падение: воздух(1.0) -> стекло(1.5)
    layers = [Layer(1.0, 0.0), Layer(1.5, 0.0)]
    got = power_rt(layers, lam, 0.0, "s")
    exp = fresnel_normal(1.0, 1.5)
    check("одна граница, нормальное падение, s-pol", got, exp)
    got_p = power_rt(layers, lam, 0.0, "p")
    check("одна граница, нормальное падение, p-pol", got_p, exp)

    # 2) Одна граница, угол 30 градусов, обе поляризации
    got_s30 = power_rt(layers, lam, 30.0, "s")
    exp_s30 = fresnel_oblique(1.0, 1.5, 30.0, "s")
    check("одна граница, 30 град, s-pol", got_s30, exp_s30)
    got_p30 = power_rt(layers, lam, 30.0, "p")
    exp_p30 = fresnel_oblique(1.0, 1.5, 30.0, "p")
    check("одна граница, 30 град, p-pol", got_p30, exp_p30)

    # 3) Тонкая беспотерьная плёнка n=2.0, d=100nm между воздухом и стеклом (n=1.5), нормальное падение
    n0, n1, n2 = 1.0, 2.0, 1.5
    d = 100.0
    layers_film = [Layer(n0, 0.0), Layer(n1, d), Layer(n2, 0.0)]
    got_film = power_rt(layers_film, lam, 0.0, "s")
    exp_film = airy_thin_film(n0, n1, n2, d, lam)
    check("тонкая плёнка (формула Эйри), нормальное падение", got_film, exp_film, tol=1e-6)

    # 4) Проверка энергетического баланса R+T=1 для беспотерьного случая на сетке длин волн
    import numpy as np
    max_imbalance = 0.0
    for lam_scan in np.arange(400.0, 1000.0, 17.0):
        R, T = power_rt(layers_film, float(lam_scan), 0.0, "s")
        max_imbalance = max(max_imbalance, abs(R + T - 1.0))
    print(f"[OK] максимальное |R+T-1| для беспотерьной плёнки по сетке длин волн: {max_imbalance:.2e}")
    assert max_imbalance < 1e-9

    print("\nВсе проверки TMM против независимых формул Френеля/Эйри прошли.\n")

    # 5) Нестрогий первый взгляд: однородный (НЕ периодический) слой с
    # eps_eff из composite_ema (Au, C=10%, n_host=1.77) между подложкой
    # n=1.45 и воздухом - имитация "если бы решётку размазать в сплошную
    # плёнку". Это НЕ решётка и не должно совпадать со статьёй, только
    # sanity-check, что резонанс формулы (1) вообще проявляется в R/T/P
    # плоского слоя качественно на своём месте.
    print("--- Нестрогий прогон: однородный слой с eps_eff(Au, C=10%) между n=1.45 и воздухом ---")
    n_host = 1.77
    eps_h = complex(n_host ** 2, 0.0)
    C = 0.10
    d_layer = 100.0  # нм, толщина слоя из статьи
    results = []
    for lam_scan in np.arange(400.0, 900.0, 5.0):
        eps_p = ema.eps_au(float(lam_scan))
        eps_eff = ema.composite_from_eps(eps_h, eps_p, C)
        # ВАЖНО: composite_ema хранит eps в конвенции eps=eps'-i*eps''
        # (см. его докстринг), а tmm.py/_kz() написан в стандартной
        # фотонной конвенции n=n'+ik, eps=eps'+i*eps'' (k>0 <=> потери).
        # Это одна и та же физика, разная бухгалтерия знака (e^{+iwt} vs
        # e^{-iwt}) - конвертируем сопряжением ПЕРЕД тем, как отдать eps в
        # tmm.py, иначе знак мнимой части перепутается и появится
        # нефизичное усиление вместо поглощения (см. README.md).
        n_eff = cmath.sqrt(eps_eff.conjugate())
        assert n_eff.imag >= -1e-9, n_eff  # лоссовая среда: k=Im(n)>=0 в стандартной конвенции
        layers_demo = [Layer(1.45, 0.0), Layer(n_eff, d_layer), Layer(1.0, 0.0)]
        R, T = power_rt(layers_demo, float(lam_scan), 0.0, "s")
        P = 1.0 - R - T
        assert P > -1e-6, f"P={P} < 0 at {lam_scan} nm - нефизичное усиление, конвенция знака сломана"
        results.append((lam_scan, R, T, P))

    peak = max(results, key=lambda row: row[3])
    print(f"Максимум потерь P={peak[3]:.3f} на lambda={peak[0]:.0f} nm (R={peak[1]:.3f}, T={peak[2]:.3f})")
    print("Не путать с результатом решётки: это план-волна на сплошной плёнке без периодичности/дифракции,")
    print("используется только как первичный sanity-check самой TMM-машины перед добавлением решётки.")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    lam_arr = [row[0] for row in results]
    R_arr = [row[1] for row in results]
    T_arr = [row[2] for row in results]
    P_arr = [row[3] for row in results]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(lam_arr, R_arr, label="R")
    ax.plot(lam_arr, T_arr, label="T")
    ax.plot(lam_arr, P_arr, label="P=1-R-T")
    ax.set_xlabel("wavelength, nm")
    ax.set_ylabel("power fraction")
    ax.set_title("TMM (не решётка): однородный слой eps_eff(Au,C=10%), d=100nm, n_sub=1.45")
    ax.legend()
    fig.tight_layout()
    out_dir = Path(__file__).resolve().parents[1] / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    png_path = out_dir / "tmm_homogeneous_layer_demo.png"
    fig.savefig(png_path, dpi=150)
    print(f"PNG записан: {png_path}")


if __name__ == "__main__":
    main()
