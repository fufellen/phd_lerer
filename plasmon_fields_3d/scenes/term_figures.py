"""Иллюстрации для заметок-терминов ваута, посчитанные тем же провалидированным ядром.

Две картинки, каждая показывает ровно тот эффект, ради которого заметка написана.

1. Разложение по собственным модам (EME). Основная мода узкого диэлектрического
   слоя раскладывается по модам широкого. Видно и распределение мощности по
   модам, и то, что набор направляемых мод неполон: сумма долей до единицы не
   дотягивает, и недостача - это излучение.

2. Биортогональное перекрытие. Две моды плёнки золота ортогональны в
   несопряжённой форме и НЕ ортогональны в сопряжённой. Развёртка по мнимой
   части проницаемости металла показывает, что расхождение появляется ровно
   тогда, когда включаются потери.

Обе картинки проверяются числами: скрипт печатает контрольные значения и падает,
если они разошлись с ожидаемыми.

Запуск:
    python plasmon_fields_3d/scenes/term_figures.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT.parent / "lrspp_coupling"))

from fields3d import materials as mat  # noqa: E402
from slabmodes import (  # noqa: E402
    Stack,
    find_modes,
    mode_fields,
    overlap_power,
    solve_mode,
    trapz,
)

OUT = ROOT / "results"
OUT.mkdir(parents=True, exist_ok=True)

LAMBDA_UM = 1.55
K0 = mat.k0_from_lambda(LAMBDA_UM)


# ============================================================ EME

N_CLAD = 1.45
N_CORE = 1.60
EPS_CLAD_D = complex(N_CLAD**2, 0.0)
EPS_CORE_D = complex(N_CORE**2, 0.0)

T_NARROW = 0.8   # мкм, одномодовый слой
T_WIDE = 4.0     # мкм, многомодовый слой


def slab(t_um: float) -> Stack:
    return Stack(eps=(EPS_CLAD_D, EPS_CORE_D, EPS_CLAD_D), thickness=(t_um,),
                 names=("обкладка", "сердцевина", "обкладка"))


def v_parameter(t_um: float) -> float:
    return float(K0 * t_um * np.sqrt(N_CORE**2 - N_CLAD**2))


def eme_figure(path_png: Path) -> dict:
    st_n, st_w = slab(T_NARROW), slab(T_WIDE)

    n_narrow = solve_mode(st_n, K0, complex(0.5 * (N_CORE + N_CLAD), 0.0))
    guesses = [complex(x, 0.0) for x in np.linspace(N_CLAD + 1e-4, N_CORE - 1e-4, 200)]
    modes_wide = find_modes(st_w, K0, guesses, n_min=N_CLAD, n_max=N_CORE)

    # общая сетка; обе структуры центрируются по своей середине
    z = np.linspace(-14.0, 14.0, 240001)
    h_n = mode_fields(n_narrow, st_n, K0, z + 0.5 * T_NARROW)[0]
    eps_n = st_n.eps_at(z + 0.5 * T_NARROW)

    etas, profiles = [], []
    for nw in modes_wide:
        h_w = mode_fields(nw, st_w, K0, z + 0.5 * T_WIDE)[0]
        eps_w = st_w.eps_at(z + 0.5 * T_WIDE)
        etas.append(overlap_power((h_n, n_narrow, eps_n), (h_w, nw, eps_w), z))
        profiles.append(h_w)
    etas = np.array(etas)
    captured = float(etas.sum())

    fig, axes = plt.subplots(1, 3, figsize=(16.0, 4.6))

    ax = axes[0]
    ax.plot(np.real(h_n) / np.abs(h_n).max(), z, lw=2.2, color="k",
            label=f"мода узкого слоя\n{T_NARROW:g} мкм, V = {v_parameter(T_NARROW):.2f}")
    for i, (nw, h_w) in enumerate(zip(modes_wide, profiles)):
        ax.plot(np.real(h_w) / np.abs(h_w).max(), z, lw=1.2,
                label=f"мода {i} широкого, n = {nw.real:.4f}")
    ax.axhspan(-0.5 * T_WIDE, 0.5 * T_WIDE, color="tab:blue", alpha=0.12)
    ax.axhspan(-0.5 * T_NARROW, 0.5 * T_NARROW, color="tab:blue", alpha=0.22)
    ax.set_ylim(-4.5, 4.5)
    ax.set_xlabel(r"$H_y$, нормировано")
    ax.set_ylabel("z, мкм")
    ax.set_title(f"Что во что раскладывается\nширокий слой {T_WIDE:g} мкм, "
                 f"V = {v_parameter(T_WIDE):.2f}")
    ax.legend(fontsize=7, loc="upper right")
    ax.grid(alpha=0.3)

    ax = axes[1]
    idx = np.arange(len(etas))
    ax.bar(idx, 100 * etas, color="tab:blue", label="направляемые моды")
    ax.bar([len(etas)], [100 * (1.0 - captured)], color="tab:red",
           label="не покрыто набором")
    ax.set_xticks(list(idx) + [len(etas)])
    ax.set_xticklabels([f"мода {i}" for i in idx] + ["излучение"], fontsize=8)
    ax.set_ylabel("доля мощности, %")
    ax.set_title("Куда уходит мощность на стыке")
    for i, e in enumerate(etas):
        ax.text(i, 100 * e + 1.5, f"{100 * e:.1f}", ha="center", fontsize=8)
    ax.text(len(etas), 100 * (1 - captured) + 1.5, f"{100 * (1 - captured):.1f}",
            ha="center", fontsize=8)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")

    ax = axes[2]
    cum = np.cumsum(etas)
    ax.plot(np.arange(1, len(cum) + 1), 100 * cum, "o-", color="tab:blue")
    ax.axhline(100.0, color="k", ls="--", lw=1, label="полный баланс")
    ax.axhline(100 * captured, color="tab:red", ls=":", lw=1.2,
               label=f"предел набора {100 * captured:.1f} %")
    ax.set_xlabel("сколько мод включено в набор")
    ax.set_ylabel("накопленная доля мощности, %")
    ax.set_title("Сходимость по числу мод\nи почему одних направляемых мало")
    ax.set_xticks(np.arange(1, len(cum) + 1))
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    fig.suptitle(
        f"Разложение по собственным модам: стык слоёв {T_NARROW:g} и {T_WIDE:g} мкм, "
        f"n = {N_CORE} в обкладке {N_CLAD}, {LAMBDA_UM * 1000:.0f} нм",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(path_png, dpi=150)
    plt.close(fig)

    return {
        "n_narrow": n_narrow,
        "modes_wide": modes_wide,
        "etas": etas,
        "captured": captured,
    }


# ================================================ биортогональность

EPS_D = mat.ZPU450
EPS_AU = mat.AU_PARK_1550
T_AU = 0.014


DN_MAX = 0.002   # предельная асимметрия обкладок, при которой LR-мода ещё связана


def film(loss_scale: float = 1.0, dn: float = 0.0) -> Stack:
    """Плёнка золота с управляемыми потерями металла и асимметрией обкладок.

    `loss_scale` умножает мнимую часть проницаемости золота: единица - реальный
    металл, ноль - искусственный металл без поглощения. `dn` добавляется к
    показателю верхней обкладки и нарушает зеркальную симметрию структуры.

    Оба параметра нужны, и это главный смысл картинки. В симметричной структуре
    LR- и SR-моды имеют противоположную чётность, поэтому подынтегральное
    выражение нечётно и ОБЕ формы перекрытия обращаются в нуль по симметрии -
    различить их на таком примере невозможно. Разница между формами проявляется
    только когда сняты обе защиты: симметрия и эрмитовость.
    """
    eps_m = complex(EPS_AU.real, EPS_AU.imag * loss_scale)
    eps_top = complex((float(np.sqrt(EPS_D).real) + dn) ** 2, 0.0)
    return Stack(eps=(EPS_D, eps_m, eps_top), thickness=(T_AU,),
                 names=("ZPU450", "Au", "ZPU450 + dn"))


def overlap_grid(half_um: float = 30.0) -> np.ndarray:
    """Сетка со сгущением у плёнки: мелкий шаг в металле, широкое окно в обкладках.

    Равномерная сетка здесь не годится. Окно должно покрывать много длин
    спадания слабо связанной LR-моды, иначе обрезанный хвост сам даёт вклад
    порядка 1e-6 и маскирует настоящую ортогональность; при этом внутри 14 нм
    металла нужен шаг в доли нанометра. Совместить это на равномерной сетке -
    миллионы узлов без всякой нужды.
    """
    lo = np.linspace(-half_um, -0.6, 60000)
    mid = np.linspace(-0.6, 0.6 + T_AU, 30000)
    hi = np.linspace(0.6 + T_AU, half_um, 60000)
    return np.unique(np.concatenate([lo, mid, hi]))


def residuals(h1: np.ndarray, h2: np.ndarray, eps_z: np.ndarray, z: np.ndarray) -> tuple[float, float]:
    """Несопряжённая и сопряжённая формы перекрытия, нормированные на свои нормы.

    Несопряжённая - та, что обращается в нуль у неэрмитовой задачи:
        int H1 H2 / eps dz.
    Сопряжённая - привычная из эрмитова случая:
        int H1 conj(H2) / eps dz.
    """
    unc = abs(trapz(h1 * h2 / eps_z, z))
    unc_d = np.sqrt(abs(trapz(h1 * h1 / eps_z, z)) * abs(trapz(h2 * h2 / eps_z, z)))
    con = abs(trapz(h1 * np.conj(h2) / eps_z, z))
    con_d = np.sqrt(abs(trapz(h1 * np.conj(h1) / eps_z, z)) * abs(trapz(h2 * np.conj(h2) / eps_z, z)))
    return float(unc / unc_d), float(con / con_d)


def _pair(st: Stack, seed_lr: complex, seed_sr: complex) -> tuple[complex, complex, bool]:
    """Пара мод плёнки и признак того, что обе ещё связаны."""
    lr = solve_mode(st, K0, seed_lr)
    sr = solve_mode(st, K0, seed_sr)
    n_max = max(float(np.sqrt(st.eps[0]).real), float(np.sqrt(st.eps[-1]).real))
    return lr, sr, bool(lr.real > n_max and sr.real > n_max)


def biorthogonal_figure(path_png: Path) -> dict:
    z = overlap_grid(30.0)

    # опорный случай: асимметрия есть, потери есть
    st = film(1.0, DN_MAX)
    n_lr, n_sr, bound = _pair(st, complex(1.4512, 1e-5), complex(1.5599, 2e-2))
    h_lr = mode_fields(n_lr, st, K0, z)[0]
    h_sr = mode_fields(n_sr, st, K0, z)[0]
    r_unc, r_con = residuals(h_lr, h_sr, st.eps_at(z), z)

    # симметричный случай для сравнения: обе формы гибнут по чётности
    st_sym = film(1.0, 0.0)
    s_lr, s_sr, _ = _pair(st_sym, complex(1.4512, 1e-5), complex(1.5599, 2e-2))
    sym_unc, sym_con = residuals(
        mode_fields(s_lr, st_sym, K0, z)[0], mode_fields(s_sr, st_sym, K0, z)[0],
        st_sym.eps_at(z), z,
    )

    # развёртка по асимметрии при полных потерях
    dn_arr = np.linspace(0.0, DN_MAX, 21)
    seed = (complex(1.4512, 1e-5), complex(1.5599, 2e-2))
    sweep_dn = []
    for dn in dn_arr:
        st_d = film(1.0, float(dn))
        lr, sr, ok = _pair(st_d, *seed)
        seed = (lr, sr)
        u, c = residuals(mode_fields(lr, st_d, K0, z)[0], mode_fields(sr, st_d, K0, z)[0],
                         st_d.eps_at(z), z)
        sweep_dn.append((float(dn), u, c, ok))

    # развёртка по потерям при фиксированной асимметрии
    scales = np.linspace(1.0, 0.0, 21)
    seed = (n_lr, n_sr)
    sweep_loss = []
    for s in scales:
        st_s = film(float(s), DN_MAX)
        lr, sr, ok = _pair(st_s, *seed)
        seed = (lr, sr)
        u, c = residuals(mode_fields(lr, st_s, K0, z)[0], mode_fields(sr, st_s, K0, z)[0],
                         st_s.eps_at(z), z)
        sweep_loss.append((float(s), u, c, ok))

    fig, axes = plt.subplots(1, 3, figsize=(16.0, 4.7))

    # Первая панель показывает, что именно значит «обращается в нуль»: накопленный
    # интеграл по сечению. У несопряжённой формы он набирает величину у плёнки и
    # затем возвращается к нулю, у сопряжённой - выходит на постоянное значение.
    from scipy.integrate import cumulative_trapezoid

    eps_z = st.eps_at(z)
    run_u = cumulative_trapezoid(h_lr * h_sr / eps_z, z, initial=0.0)
    run_c = cumulative_trapezoid(h_lr * np.conj(h_sr) / eps_z, z, initial=0.0)

    ax = axes[0]
    ax.plot(z, np.abs(run_u) / max(np.abs(run_u).max(), 1e-300), color="tab:green",
            label=r"несопряжённая $\int H_1H_2/\varepsilon$")
    ax.plot(z, np.abs(run_c) / max(np.abs(run_c).max(), 1e-300), color="tab:red",
            label=r"сопряжённая $\int H_1H_2^*/\varepsilon$")
    ax.axvline(0.0, color="goldenrod", lw=2.5, alpha=0.8)
    # логарифмическая ось обязательна: на линейной хвосты 0.00002 и 0.02
    # неразличимы, а разница между ними и есть весь смысл картинки
    ax.set_yscale("log")
    ax.set_xlim(-8.0, 8.0)
    ax.set_ylim(1e-7, 3.0)
    ax.set_xlabel("верхний предел интегрирования z, мкм")
    ax.set_ylabel("накопленный интеграл, нормирован на свой максимум")
    ax.set_title("Что значит «обращается в нуль»\n"
                 f"плёнка Au в жёлтой линии, dn = {DN_MAX}")
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(alpha=0.3, which="both")

    ax = axes[1]
    d_arr = np.array([p[0] for p in sweep_dn])
    u_arr = np.array([max(p[1], 1e-13) for p in sweep_dn])
    c_arr = np.array([max(p[2], 1e-13) for p in sweep_dn])
    ax.semilogy(d_arr, c_arr, "s-", ms=3.5, color="tab:red",
                label=r"сопряжённая $\int H_1H_2^*/\varepsilon\,dz$")
    ax.semilogy(d_arr, u_arr, "o-", ms=3.5, color="tab:green",
                label=r"несопряжённая $\int H_1H_2/\varepsilon\,dz$")
    ax.set_xlabel("асимметрия обкладок dn")
    ax.set_ylabel("перекрытие, нормировано")
    ax.set_title("Симметрия скрывает различие\n(при dn = 0 гибнут обе формы)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    ax = axes[2]
    s_arr = np.array([p[0] for p in sweep_loss])
    u2 = np.array([max(p[1], 1e-13) for p in sweep_loss])
    c2 = np.array([max(p[2], 1e-13) for p in sweep_loss])
    ax.semilogy(s_arr, c2, "s-", ms=3.5, color="tab:red", label="сопряжённая")
    ax.semilogy(s_arr, u2, "o-", ms=3.5, color="tab:green", label="несопряжённая")
    ax.set_xlabel(r"множитель при $\mathrm{Im}\,\varepsilon$ золота")
    ax.set_ylabel("перекрытие, нормировано")
    ax.set_title(f"При dn = {DN_MAX} различие задают потери\n(без потерь задача эрмитова)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    fig.suptitle(
        "Биортогональность: в нуль обращается несопряжённое перекрытие, а не привычное "
        f"сопряжённое, и увидеть это можно только без зеркальной симметрии ({LAMBDA_UM * 1000:.0f} нм)",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(path_png, dpi=150)
    plt.close(fig)

    return {"n_lr": n_lr, "n_sr": n_sr, "bound": bound,
            "unconjugated": r_unc, "conjugated": r_con,
            "sym_unconjugated": sym_unc, "sym_conjugated": sym_con,
            "sweep_dn": sweep_dn, "sweep_loss": sweep_loss}


def main() -> int:
    failures = 0

    print("Иллюстрация 1: разложение по собственным модам")
    eme = eme_figure(OUT / "eme_mode_expansion.png")
    print(f"  узкий слой {T_NARROW:g} мкм: n_eff = {eme['n_narrow'].real:.6f}, "
          f"V = {v_parameter(T_NARROW):.3f}")
    print(f"  широкий слой {T_WIDE:g} мкм: найдено мод {len(eme['modes_wide'])}, "
          f"V = {v_parameter(T_WIDE):.3f}")
    for i, (nw, e) in enumerate(zip(eme["modes_wide"], eme["etas"])):
        print(f"    мода {i}: n_eff = {nw.real:.6f}, доля мощности {100 * e:.2f} %")
    print(f"  сумма по направляемым модам {100 * eme['captured']:.2f} %, "
          f"остальное - излучение")

    # число мод у симметричного слоя: floor(V/pi) + 1
    expected = int(np.floor(v_parameter(T_WIDE) / np.pi)) + 1
    ok = len(eme["modes_wide"]) == expected
    print(f"  [{'OK' if ok else 'СБОЙ'}] число мод широкого слоя совпадает с V/pi: "
          f"найдено {len(eme['modes_wide'])}, ожидалось {expected}")
    failures += 0 if ok else 1

    # чётный стык возбуждает только чётные моды: нечётные должны быть подавлены
    odd = eme["etas"][1::2]
    ok = float(np.max(odd)) < 1e-6 if odd.size else True
    print(f"  [{'OK' if ok else 'СБОЙ'}] симметричный стык не возбуждает нечётные моды: "
          f"максимум {100 * float(np.max(odd)) if odd.size else 0:.2e} %")
    failures += 0 if ok else 1

    ok = 0.5 < eme["captured"] < 1.0
    print(f"  [{'OK' if ok else 'СБОЙ'}] набор направляемых мод неполон, но основную часть "
          f"мощности покрывает: {100 * eme['captured']:.2f} %")
    failures += 0 if ok else 1

    print()
    print("Иллюстрация 2: биортогональное перекрытие")
    bio = biorthogonal_figure(OUT / "biorthogonal_overlap.png")
    print(f"  n_LR = {bio['n_lr'].real:.6f} + {bio['n_lr'].imag:.3e}i")
    print(f"  n_SR = {bio['n_sr'].real:.6f} + {bio['n_sr'].imag:.3e}i")
    print(f"  симметричная структура: несопряжённое {bio['sym_unconjugated']:.2e}, "
          f"сопряжённое {bio['sym_conjugated']:.2e}")
    print(f"  асимметричная (dn = {DN_MAX}): несопряжённое {bio['unconjugated']:.2e}, "
          f"сопряжённое {bio['conjugated']:.2e}")

    ok = bio["bound"]
    print(f"  [{'OK' if ok else 'СБОЙ'}] обе моды при dn = {DN_MAX} ещё связаны")
    failures += 0 if ok else 1

    # у симметричной структуры обе формы гибнут по чётности - на ней различие не видно
    ok = bio["sym_conjugated"] < 1e-8 and bio["sym_unconjugated"] < 1e-8
    print(f"  [{'OK' if ok else 'СБОЙ'}] в симметричной структуре обе формы обращаются в нуль "
          f"по чётности, и различить их нельзя")
    failures += 0 if ok else 1

    ok = bio["unconjugated"] < 1e-6
    print(f"  [{'OK' if ok else 'СБОЙ'}] без симметрии несопряжённая форма всё равно "
          f"обращается в нуль: {bio['unconjugated']:.2e}")
    failures += 0 if ok else 1

    ratio = bio["conjugated"] / max(bio["unconjugated"], 1e-300)
    ok = ratio > 1e3
    print(f"  [{'OK' if ok else 'СБОЙ'}] сопряжённая форма в нуль НЕ обращается: "
          f"она больше несопряжённой в {ratio:.0e} раз")
    failures += 0 if ok else 1

    lossless = bio["sweep_loss"][-1]
    ok = lossless[2] < 1e-6
    print(f"  [{'OK' if ok else 'СБОЙ'}] без потерь задача эрмитова и сопряжённая форма "
          f"тоже обращается в нуль: {lossless[2]:.2e} при множителе {lossless[0]:.2f}")
    failures += 0 if ok else 1

    print()
    print(f"Картинки записаны в {OUT}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
