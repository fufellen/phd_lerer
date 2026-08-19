"""Сцена 2. Связанный плазмон: две границы тонкой плёнки работают как одна система.

Что показывается. Если металл сделать тоньше глубины проникновения поля, плазмоны
двух его границ перестают быть независимыми. Они связываются и дают пару мод с
общим полем на всю плёнку:

  - длиннопробежная (LR-SPP): E_z по разные стороны плёнки в противофазе, поле в
    металле почти гасится, поглощение падает, и мода живёт на миллиметры;
  - короткопробежная (SR-SPP): E_z синфазна, поле в металл вдавливается, и мода
    гибнет за единицы микрометров.

Отсюда главный вывод для техники: выигрыш в дальности покупается не материалом, а
геометрией, и платой служит слабое удержание. У LR-моды поле расплывается в
диэлектрик на десятки микрометров, поэтому она плохо совместима с компактными
элементами - и именно это делает нужным преобразователь размера моды из сцены 3.

Чем проверяется:
  - оба корня против уравнений Майера (2.29a,b), решённых независимо;
  - предел толстой плёнки: обе моды сходятся к плазмону одной границы;
  - несопряжённая ортогональность двух мод одной структуры;
  - невязка уравнений Максвелла и её порядок по шагу сетки.

Запуск:
    python plasmon_fields_3d/scenes/scene2_coupled_spp.py
"""

from __future__ import annotations

import csv
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
from fields3d import render  # noqa: E402
from fields3d.section import FieldSection, residual_convergence  # noqa: E402
from slabmodes import (  # noqa: E402
    Stack,
    orthogonality_residual,
    propagation_length_um,
    propagation_loss_db_per_cm,
    solve_mode,
)
from slabmodes.tmm import decaying_sqrt  # noqa: E402

OUT = ROOT / "results"
OUT.mkdir(parents=True, exist_ok=True)

LAMBDA_UM = 1.55
K0 = mat.k0_from_lambda(LAMBDA_UM)
EPS_D = mat.ZPU450
EPS_AU = mat.AU_PARK_1550
N_D = float(np.sqrt(EPS_D).real)
T_FILM_UM = 0.014          # Park 2011: плёнка золота 14 нм
PAPER_LR_NEFF = 1.4512     # Park 2011, рабочая мода ДМД


def film(t_um: float) -> Stack:
    return Stack(eps=(EPS_D, EPS_AU, EPS_D), thickness=(t_um,), names=("ZPU450", "Au", "ZPU450"))


def single_interface_neff() -> complex:
    n = np.sqrt(EPS_AU * EPS_D / (EPS_AU + EPS_D))
    return complex(n if n.imag >= 0 else np.conj(n))


def _safe_tanh(w: complex) -> complex:
    """tanh без переполнения: при больших Re аргумента он насыщается на +-1.

    У толстой плёнки gamma_m * a достигает сотен, и прямой вызов np.tanh даёт
    переполнение и NaN, из-за чего ньютоновский шаг уходит в бесконечность.
    """
    if w.real > 20.0:
        return complex(1.0, 0.0)
    if w.real < -20.0:
        return complex(-1.0, 0.0)
    return complex(np.tanh(w))


def maier_residual(kind: str, neff: complex, t_um: float) -> complex:
    """Уравнения Майера (2.29a,b) - независимая от матрицы переноса запись.

    Для симметричной структуры дисперсионное уравнение распадается по чётности:

        tanh(gamma_m a) = -gamma_d eps_m / (gamma_m eps_d)   - ветвь a,
        tanh(gamma_m a) = -gamma_m eps_d / (gamma_d eps_a)   - ветвь b,

    где a - половина толщины плёнки. Совпадение их корней с корнями общей
    матрицы переноса проверяет ядро, а не повторяет его.
    """
    half = 0.5 * t_um
    gd = K0 * decaying_sqrt(neff * neff - EPS_D)
    gm = K0 * decaying_sqrt(neff * neff - EPS_AU)
    th = _safe_tanh(gm * half)
    if kind == "a":
        return th + gd * EPS_AU / (gm * EPS_D)
    return th + gm * EPS_D / (gd * EPS_AU)


def maier_root(kind: str, t_um: float, guess: complex) -> complex:
    z = complex(guess)
    for _ in range(300):
        v = maier_residual(kind, z, t_um)
        if not np.isfinite(v.real) or not np.isfinite(v.imag):
            return complex("nan")
        h = max(1e-10, abs(z) * 1e-8)
        d = (maier_residual(kind, z + h, t_um) - v) / h
        if d == 0:
            break
        step = v / d
        z -= step
        if abs(step) < 1e-15:
            break
    return z.conjugate() if z.imag < 0 else z


# Лестница толщин, по которой ведётся продолжение ветвей. Начинать надо с
# тонкой плёнки: там моды разведены далеко и приближение для каждой очевидно.
# У толстой плёнки они вырождаются к плазмону одной границы, и начатый оттуда
# поиск обеих ветвей давал бы один и тот же корень.
LADDER = np.geomspace(0.004, 1.0, 90)


def trace_branches() -> list[dict]:
    """Обе ветви прослеживаются продолжением по толщине плёнки.

    На каждом шаге приближением служит корень с предыдущего шага, поэтому ветвь
    не перескакивает на соседнюю: шаг по толщине мал, и корень смещается мало.
    Найденный корень затем уточняется независимой матрицей переноса, и оба
    результата сохраняются - расхождение между ними и есть проверка.
    """
    rows: list[dict] = []
    guess = {"a": complex(1.4501, 1e-7), "b": complex(2.35, 0.12)}
    for t in LADDER:
        row = {"t_um": float(t)}
        st = film(float(t))
        for kind, key in (("a", "lr"), ("b", "sr")):
            r_maier = maier_root(kind, float(t), guess[kind])
            if np.isfinite(r_maier.real):
                guess[kind] = r_maier
            r_tmm = solve_mode(st, K0, r_maier)
            row[f"{key}_maier"] = r_maier
            row[f"{key}_tmm"] = r_tmm
        rows.append(row)
    return rows


def at_thickness(rows: list[dict], t_um: float) -> dict:
    """Ближайшая к заданной толщине запись лестницы, уточнённая точно на ней."""
    i = int(np.argmin([abs(r["t_um"] - t_um) for r in rows]))
    st = film(t_um)
    out = {"t_um": t_um}
    for kind, key in (("a", "lr"), ("b", "sr")):
        seed = rows[i][f"{key}_maier"]
        r_maier = maier_root(kind, t_um, seed)
        out[f"{key}_maier"] = r_maier
        out[key] = solve_mode(st, K0, r_maier)
    return out


def parity(stack: Stack, neff: complex) -> str:
    """Чётность H_y относительно середины плёнки."""
    edges = stack.interfaces()
    c = 0.5 * (edges[0] + edges[-1])
    d = 0.3 * stack.total_thickness
    from slabmodes import mode_fields

    h = mode_fields(neff, stack, K0, np.array([c - d, c + d]))[0]
    return "чётная" if np.real(h[0] * np.conj(h[1])) > 0 else "нечётная"


def thickness_sweep(traced: list[dict], path_csv: Path, path_png: Path, single: complex) -> None:
    rows = [
        {
            "t_nm": r["t_um"] * 1000.0,
            "lr_re": r["lr_tmm"].real,
            "lr_im": r["lr_tmm"].imag,
            "sr_re": r["sr_tmm"].real,
            "sr_im": r["sr_tmm"].imag,
            "lr_L_um": propagation_length_um(r["lr_tmm"], LAMBDA_UM),
            "sr_L_um": propagation_length_um(r["sr_tmm"], LAMBDA_UM),
            "maier_vs_tmm_lr": abs(r["lr_tmm"] - r["lr_maier"]),
            "maier_vs_tmm_sr": abs(r["sr_tmm"] - r["sr_maier"]),
        }
        for r in traced
    ]

    with path_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    t_nm = np.array([r["t_nm"] for r in rows])
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3))

    ax = axes[0]
    ax.semilogx(t_nm, [r["lr_re"] for r in rows], "o-", ms=3, label="LR-SPP")
    ax.semilogx(t_nm, [r["sr_re"] for r in rows], "s-", ms=3, label="SR-SPP")
    ax.axhline(single.real, color="k", ls="--", lw=1, label="плазмон одной границы")
    ax.axhline(N_D, color="0.6", ls=":", lw=1, label="показатель обкладки")
    ax.set_xlabel("толщина плёнки, нм")
    ax.set_ylabel(r"$\mathrm{Re}\,n_\mathrm{eff}$")
    ax.set_title("расщепление при связывании границ")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    ax = axes[1]
    ax.loglog(t_nm, [r["lr_L_um"] for r in rows], "o-", ms=3, label="LR-SPP")
    ax.loglog(t_nm, [r["sr_L_um"] for r in rows], "s-", ms=3, label="SR-SPP")
    ax.axhline(propagation_length_um(single, LAMBDA_UM), color="k", ls="--", lw=1,
               label="одна граница")
    ax.axvline(T_FILM_UM * 1000, color="tab:red", ls=":", lw=1.2, label="плёнка сцены, 14 нм")
    ax.set_xlabel("толщина плёнки, нм")
    ax.set_ylabel("длина распространения, мкм")
    ax.set_title("дальность против толщины")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    fig.suptitle("Связанные плазмоны тонкой плёнки Au в ZPU450, 1550 нм", fontsize=10)
    fig.tight_layout()
    fig.savefig(path_png, dpi=150)
    plt.close(fig)


def build_section(neff: complex, nx: int, nz: int, x_max: float, z_half: float) -> FieldSection:
    st = film(T_FILM_UM)
    z = np.linspace(-z_half, z_half + T_FILM_UM, nz)
    x = np.linspace(0.0, x_max, nx)
    return FieldSection.from_mode(st, neff, LAMBDA_UM, x, z)


def plot_profiles(pair: dict, path: Path) -> None:
    st = film(T_FILM_UM)
    z = np.linspace(-6.0, 6.0 + T_FILM_UM, 4001)
    from slabmodes import mode_fields

    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.4))
    for ax, key, name in ((axes[0], "lr", "LR-SPP"), (axes[1], "sr", "SR-SPP")):
        hy, ex, ez = mode_fields(pair[key], st, K0, z)
        scale = max(np.abs(ez).max(), 1e-300)
        ax.plot(np.real(ez) / scale, z, label=r"$\mathrm{Re}\,E_z$")
        ax.plot(np.abs(hy) / max(np.abs(hy).max(), 1e-300), z, ls="--", label=r"$|H_y|$")
        ax.axhspan(0.0, T_FILM_UM, color="goldenrod", alpha=0.5)
        ax.set_ylim(-4.0, 4.0)
        ax.set_xlabel("нормированная амплитуда")
        ax.set_ylabel("z, мкм")
        ax.set_title(
            f"{name}\n"
            f"n_eff = {pair[key].real:.6f} + {pair[key].imag:.2e}i\n"
            f"L = {propagation_length_um(pair[key], LAMBDA_UM):.3g} мкм"
        )
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    ax = axes[2]
    z2 = np.linspace(-0.02, 0.02 + T_FILM_UM, 2001)
    for key, name in (("lr", "LR-SPP"), ("sr", "SR-SPP")):
        ez = mode_fields(pair[key], st, K0, z2)[2]
        ax.plot(np.real(ez) / max(np.abs(ez).max(), 1e-300), z2 * 1000, label=name)
    ax.axhspan(0.0, T_FILM_UM * 1000, color="goldenrod", alpha=0.5)
    ax.set_xlabel(r"$\mathrm{Re}\,E_z$, норм.")
    ax.set_ylabel("z, нм")
    ax.set_title("поле внутри плёнки 14 нм\n(причина разной дальности)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def render_3d(pair: dict) -> None:
    x_max = 12.0
    z_half = 3.0
    y_half = 2.4

    for key, name in (("lr", "LR-SPP"), ("sr", "SR-SPP")):
        neff = pair[key]
        sec = build_section(neff, 561, 401, x_max, z_half)
        grid = render.build_grid(sec, y_half_um=y_half, ny=14)
        lp = propagation_length_um(neff, LAMBDA_UM)
        p = render.new_plotter(
            f"{name}: связанная мода плёнки Au 14 нм в ZPU450, 1550 нм\n"
            f"изоповерхности мгновенной E_z, задняя стенка - её срез, "
            f"жёлтая плита - золото\n"
            f"n_eff = {neff.real:.6f} + {neff.imag:.2e}i, "
            f"длина распространения {lp:.3g} мкм, "
            f"потери {propagation_loss_db_per_cm(neff, LAMBDA_UM):.3g} дБ/см"
        )
        render.add_cut_plane(p, grid, "Ez", y_um=-0.98 * y_half, title="Ez", clim_frac=0.3)
        render.add_signed_isosurfaces(p, grid, "Ez", levels=(0.4,), opacity=0.45)
        render.add_layer_boxes(
            p, [("Au", 0.0, T_FILM_UM, "goldenrod")], (0.0, x_max), y_half, opacity=0.95
        )
        render.show_axes(p)
        render.finish(
            p,
            OUT / f"scene2_{key}_mode3d.png",
            camera="iso",
            html=OUT / f"scene2_{key}_mode3d.html",
            zoom=1.3,
        )

    # ---- сравнение дальности: обе моды на одной длине трассы
    lp_lr = propagation_length_um(pair["lr"], LAMBDA_UM)
    x_far = 3.0 * propagation_length_um(pair["sr"], LAMBDA_UM)
    p = render.new_plotter(
        f"Почему длиннопробежная мода называется длиннопробежной\n"
        f"обе моды одной и той же плёнки на одной трассе {x_far:.1f} мкм\n"
        f"сзади LR-SPP (L = {lp_lr:.0f} мкм), спереди SR-SPP "
        f"(L = {propagation_length_um(pair['sr'], LAMBDA_UM):.2f} мкм)\n"
        f"цвет - огибающая |E|, нормирована на своё начальное значение"
    )
    for key, y_c in (("lr", 4.0), ("sr", -4.0)):
        sec = build_section(pair[key], 601, 241, x_far, 2.0)
        norm = sec.magnitude_e()[0].max()
        sec.ex = sec.ex / norm
        sec.ez = sec.ez / norm
        sec.hy = sec.hy / norm
        grid = render.build_grid(sec, y_half_um=1.6, ny=8, y_center_um=y_c)
        render.add_field_volume(p, grid, scalars="absE", cmap=render.SEQUENTIAL,
                                clim=(0.0, 1.0), title="|E| норм.")
    render.show_axes(p)
    render.finish(p, OUT / "scene2_range_comparison3d.png", camera="iso",
                  scale=(0.6, 1.0, 1.0), zoom=1.25)


def main() -> int:
    lines: list[str] = []
    add = lines.append
    failures = 0

    single = single_interface_neff()
    traced = trace_branches()
    pair = at_thickness(traced, T_FILM_UM)
    st = film(T_FILM_UM)

    add("Сцена 2. Связанные плазмоны тонкой металлической плёнки")
    add(f"Плёнка Au {T_FILM_UM * 1000:.0f} нм между одинаковыми слоями ZPU450 (n = {N_D}), "
        f"{LAMBDA_UM * 1000:.0f} нм")
    add("")
    add(f"  плазмон одной границы     n_eff = {single.real:.9f} + {single.imag:.6e}i, "
        f"L = {propagation_length_um(single, LAMBDA_UM):.2f} мкм")
    for key, name in (("lr", "LR-SPP"), ("sr", "SR-SPP")):
        n = pair[key]
        add(f"  {name:8s}                 n_eff = {n.real:.9f} + {n.imag:.6e}i")
        add(f"            H_y {parity(st, n):9s}   L = {propagation_length_um(n, LAMBDA_UM):.4g} мкм, "
            f"потери {propagation_loss_db_per_cm(n, LAMBDA_UM):.4g} дБ/см")
        gd = K0 * decaying_sqrt(n * n - EPS_D)
        add(f"            спад в обкладку {1.0 / abs(gd.real):.3f} мкм")
    ratio = propagation_length_um(pair["lr"], LAMBDA_UM) / propagation_length_um(pair["sr"], LAMBDA_UM)
    add(f"  отношение дальностей      {ratio:.0f}")
    add("")

    # --- проверки
    d_lr = abs(pair["lr"] - pair["lr_maier"]) / abs(pair["lr_maier"])
    ok = d_lr < 1e-9
    add(f"  [{'OK' if ok else 'СБОЙ'}] LR-мода против уравнения Майера 2.29a: отклонение {d_lr:.2e}")
    failures += 0 if ok else 1

    d_sr = abs(pair["sr"] - pair["sr_maier"]) / abs(pair["sr_maier"])
    ok = d_sr < 1e-9
    add(f"  [{'OK' if ok else 'СБОЙ'}] SR-мода против уравнения Майера 2.29b: отклонение {d_sr:.2e}")
    failures += 0 if ok else 1

    worst = max(
        max(abs(r["lr_tmm"] - r["lr_maier"]), abs(r["sr_tmm"] - r["sr_maier"])) for r in traced
    )
    ok = worst < 1e-8
    add(f"  [{'OK' if ok else 'СБОЙ'}] обе ветви совпадают с матрицей переноса на всей "
        f"лестнице из {len(traced)} толщин: худшее расхождение {worst:.2e}")
    failures += 0 if ok else 1

    thick = traced[-1]  # плёнка 1 мкм - границы уже практически независимы
    d_thick = max(abs(thick["lr_tmm"] - single), abs(thick["sr_tmm"] - single)) / abs(single)
    ok = d_thick < 1e-6
    add(f"  [{'OK' if ok else 'СБОЙ'}] предел толстой плёнки {thick['t_um'] * 1000:.0f} нм: обе моды "
        f"сходятся к плазмону одной границы, расхождение {d_thick:.2e}")
    failures += 0 if ok else 1

    from slabmodes import mode_fields

    zz = np.linspace(-30.0, 30.0 + T_FILM_UM, 400001)
    h_lr = mode_fields(pair["lr"], st, K0, zz)[0]
    h_sr = mode_fields(pair["sr"], st, K0, zz)[0]
    resid = orthogonality_residual(h_lr, h_sr, st.eps_at(zz), zz)
    ok = resid < 1e-6
    add(f"  [{'OK' if ok else 'СБОЙ'}] моды несопряжённо ортогональны: невязка {resid:.2e}")
    failures += 0 if ok else 1

    # Расчётная сетка проверки обязана разрешать саму плёнку, а не только
    # обкладку. Если внутри 14 нм металла не помещается ни одна точка с полным
    # трёхточечным шаблоном, металл целиком выпадает из проверки; при удвоении
    # сетки он в неё возвращается уже со своей крупной ошибкой, и наблюдаемый
    # порядок получается отрицательным. Шаг 0.5 нм даёт около 28 точек на плёнку.
    def build(refine: int) -> FieldSection:
        nz = 2400 * refine + 1
        nx = 20 * refine + 1
        return build_section(pair["lr"], nx, nz, 0.1, 0.6)

    conv = residual_convergence(build)
    ok = abs(conv["observed_order"] - 2.0) < 0.2
    add(f"  [{'OK' if ok else 'СБОЙ'}] невязка Максвелла падает как шаг в квадрате: "
        f"{conv['residual_coarse']:.2e} -> {conv['residual_fine']:.2e}, "
        f"порядок {conv['observed_order']:.2f}")
    failures += 0 if ok else 1

    d_paper = abs(pair["lr"].real - PAPER_LR_NEFF)
    ok = d_paper < 5e-4
    add(f"  [{'OK' if ok else 'СБОЙ'}] LR-мода против значения Park 2011: "
        f"{pair['lr'].real:.6f} против {PAPER_LR_NEFF}, разность {d_paper:.1e}")
    failures += 0 if ok else 1
    add("")

    add("Физический итог")
    add("  Связывание расщепляет один плазмон границы на пару мод. У LR-моды поле")
    add("  внутри металла почти обращается в нуль, поэтому поглощается мало и")
    add(f"  дальность растёт в {ratio:.0f} раз против SR-моды. Но та же слабая связь с")
    add(f"  металлом означает слабое удержание: поле уходит в обкладку на")
    add(f"  {1.0 / abs((K0 * decaying_sqrt(pair['lr'] ** 2 - EPS_D)).real):.1f} мкм, то есть мода в тысячи раз шире плёнки,")
    add("  которая её ведёт. Совместить такую моду с компактным диэлектрическим")
    add("  волноводом напрямую нельзя - этим занимается сцена 3.")

    plot_profiles(pair, OUT / "scene2_mode_profiles.png")
    thickness_sweep(traced, OUT / "scene2_thickness_sweep.csv",
                    OUT / "scene2_thickness_sweep.png", single)
    render_3d(pair)

    (OUT / "scene2_coupled_spp.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
