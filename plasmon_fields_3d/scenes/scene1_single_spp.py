"""Сцена 1. Обычный плазмон: поверхностный плазмон-поляритон одной границы.

Что показывается. Простейший плазмон - волна на единственной границе металла и
диэлектрика. У неё три отличительных признака, и все три видны на сцене:

1. поле спадает по обе стороны границы экспоненциально, но резко несимметрично:
   в металл на десятки нанометров, в диэлектрик на микрометры;
2. продольная компонента E_z много больше поперечной E_x, поэтому вектор E почти
   лежит вдоль направления распространения и рисует характерные петли;
3. на границе колеблется поверхностный заряд - именно он отличает плазмон от
   обычной направляемой волны и делает возможным n_eff выше показателя
   диэлектрика.

Чем проверяется:
  - n_eff против замкнутой формулы sqrt(eps_m eps_d / (eps_m + eps_d));
  - невязка уравнений Максвелла на самой сетке сцены;
  - спад мощности вдоль x против 2 Im(beta);
  - отношение |E_z / E_x| в диэлектрике против beta / kappa_d.

Запуск:
    python plasmon_fields_3d/scenes/scene1_single_spp.py
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
from fields3d import render  # noqa: E402
from fields3d.section import FieldSection, residual_convergence  # noqa: E402
from slabmodes import Stack, propagation_length_um, solve_mode  # noqa: E402
from slabmodes.tmm import decaying_sqrt  # noqa: E402

OUT = ROOT / "results"
OUT.mkdir(parents=True, exist_ok=True)

CASES = {
    "au_polymer_1550": {
        "title": "Au / полимер ZPU450, 1550 нм",
        "eps_metal": mat.AU_PARK_1550,
        "eps_diel": mat.ZPU450,
        "lambda_um": 1.55,
        "metal": "Au",
        "diel": "ZPU450",
    },
    "ag_air_633": {
        "title": "Ag / воздух, 633 нм",
        "eps_metal": mat.AG_633,
        "eps_diel": mat.AIR,
        "lambda_um": 0.633,
        "metal": "Ag",
        "diel": "воздух",
    },
}


def analytic_neff(eps_m: complex, eps_d: complex) -> complex:
    """Замкнутая формула ППП одной границы."""
    n = np.sqrt(eps_m * eps_d / (eps_m + eps_d))
    return complex(n if n.imag >= 0 else np.conj(n))


def solve_case(case: dict) -> dict:
    k0 = mat.k0_from_lambda(case["lambda_um"])
    stack = Stack(eps=(case["eps_metal"], case["eps_diel"]), thickness=())
    exact = analytic_neff(case["eps_metal"], case["eps_diel"])
    neff = solve_mode(stack, k0, exact * (1.0 + 1e-3))

    gamma_d = k0 * decaying_sqrt(neff * neff - case["eps_diel"])
    gamma_m = k0 * decaying_sqrt(neff * neff - case["eps_metal"])
    return {
        "stack": stack,
        "k0": k0,
        "neff": neff,
        "exact": exact,
        "delta_d_um": 1.0 / abs(gamma_d.real),
        "delta_m_um": 1.0 / abs(gamma_m.real),
        "l_prop_um": propagation_length_um(neff, case["lambda_um"]),
        "lambda_spp_um": case["lambda_um"] / neff.real,
        "ez_ex_ratio": abs(k0 * neff / gamma_d),
    }


def make_sections(case: dict, sol: dict) -> tuple[FieldSection, FieldSection]:
    """Ближняя сцена со структурой волны и дальняя со спадом вдоль трассы."""
    lam = case["lambda_um"]
    # сетку по z смещаем на полшага, чтобы узел не попал точно на границу:
    # там E_z разрывна, и значение в самой точке разрыва не определено
    z_up = 2.2 * sol["delta_d_um"]
    z_dn = 5.0 * sol["delta_m_um"]
    nz = 361
    z = np.linspace(-z_dn, z_up, nz)
    z = z + 0.5 * (z[1] - z[0]) - z[np.argmin(np.abs(z))]

    x_near = np.linspace(0.0, 5.0 * lam / 1.0, 481)
    near = FieldSection.from_mode(
        sol["stack"], sol["neff"], lam, x_near, z, title=case["title"] + ", структура волны"
    )

    x_far = np.linspace(0.0, 3.0 * sol["l_prop_um"], 601)
    far = FieldSection.from_mode(
        sol["stack"], sol["neff"], lam, x_far, z, title=case["title"] + ", затухание вдоль трассы"
    )
    return near, far


def plot_profiles(case: dict, sol: dict, near: FieldSection, path: Path) -> None:
    z = near.z
    i0 = 0
    hy = np.abs(near.hy[i0])
    ex = np.abs(near.ex[i0])
    ez = np.abs(near.ez[i0])
    scale = max(hy.max(), 1e-300)

    fig, axes = plt.subplots(1, 4, figsize=(17.5, 4.4))

    ax = axes[0]
    # H_y рисуется пунктиром поверх E_z намеренно: внутри одного слоя обе
    # величины имеют один и тот же профиль по z, поскольку E_z = -(k0 n/eps) H_y,
    # и сплошная линия просто закрыла бы вторую. Различаются они только скачком
    # на границе, где меняется eps
    ax.plot(ez / max(ez.max(), 1e-300), z, label=r"$|E_z|$", color="tab:green")
    ax.plot(hy / scale, z, label=r"$|H_y|$", color="tab:blue", ls="--", lw=1.6)
    ax.plot(ex / max(ez.max(), 1e-300), z, label=r"$|E_x|$", color="tab:orange")
    ax.axhline(0.0, color="k", lw=1.2)
    ax.axhspan(z.min(), 0.0, color="0.75", alpha=0.6)
    ax.set_xlabel("нормированная амплитуда")
    ax.set_ylabel("z, мкм")
    ax.set_title(f"профиль поля\n{case['metal']} снизу, {case['diel']} сверху")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1]
    ex_t, ez_t, _ = near.instantaneous(0.0)
    m = ax.pcolormesh(near.x, z, ez_t.T, cmap="RdBu_r", shading="auto")
    ax.axhline(0.0, color="k", lw=1.0)
    ax.set_xlabel("x, мкм")
    ax.set_ylabel("z, мкм")
    ax.set_title(r"мгновенная $E_z$")
    fig.colorbar(m, ax=ax, fraction=0.046)

    # вектор E: направление показывает петли, характерные для плазмона
    ax = axes[2]
    sx, sz = 26, 9
    xq = near.x[::sx]
    zsel = (z > -1.2 * sol["delta_m_um"]) & (z < 1.1 * sol["delta_d_um"])
    zq = z[zsel][::sz]
    exq = ex_t[::sx][:, zsel][:, ::sz]
    ezq = ez_t[::sx][:, zsel][:, ::sz]
    norm = np.hypot(exq, ezq)
    norm[norm == 0] = 1.0
    ax.quiver(
        xq, zq, (exq / norm).T, (ezq / norm).T,
        np.log10(np.hypot(exq, ezq) / np.hypot(exq, ezq).max()).T,
        cmap="viridis", pivot="mid", scale=26, width=0.005,
    )
    ax.axhline(0.0, color="k", lw=1.2)
    ax.axhspan(z.min(), 0.0, color="0.75", alpha=0.6)
    ax.set_xlim(near.x.min(), near.x.min() + 3.0 * sol["lambda_spp_um"])
    ax.set_xlabel("x, мкм")
    ax.set_ylabel("z, мкм")
    ax.set_title("направление вектора E\n(цвет — lg|E|)")

    ax = axes[3]
    sigma = near.surface_charge(0.0)
    ax.plot(near.x, sigma / max(abs(sigma).max(), 1e-300))
    ax.set_xlabel("x, мкм")
    ax.set_ylabel(r"$\sigma$, норм.")
    ax.set_title("поверхностный заряд на границе")
    ax.grid(alpha=0.3)

    fig.suptitle(
        f"{case['title']}: n_eff = {sol['neff'].real:.6f} + {sol['neff'].imag:.3e}i, "
        f"L = {sol['l_prop_um']:.1f} мкм",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def render_3d(case: dict, sol: dict, near: FieldSection, far: FieldSection, tag: str) -> None:
    lam = case["lambda_um"]
    y_half = 1.5 * lam
    x_rng = (near.x.min(), near.x.max())

    # ---- ближняя сцена: срез поля, гребни волны, стрелки направления E
    grid = render.build_grid(near, y_half_um=y_half, ny=14, phase_rad=0.0)
    p = render.new_plotter(
        f"{case['title']}: плазмон одной границы\n"
        f"изоповерхности - гребни и впадины мгновенной E_z, задняя стенка - её срез\n"
        f"серая плита - {case['metal']}, выше - {case['diel']}\n"
        f"вдоль оси y структура однородна, показан отрезок {2 * y_half:.1f} мкм"
    )
    # срез ставится задней стенкой, иначе он закрывает объём; ровно на границе
    # области срез вырождается в пустой набор, поэтому берётся чуть внутри
    render.add_cut_plane(p, grid, "Ez", y_um=-0.98 * y_half, title="Ez", clim_frac=0.22)
    render.add_signed_isosurfaces(p, grid, "Ez", levels=(0.4,), opacity=0.45)
    render.add_layer_boxes(p, [("металл", near.z.min(), 0.0, "dimgray")], x_rng, y_half, opacity=0.55)
    render.show_axes(p)
    render.finish(
        p,
        OUT / f"scene1_{tag}_wave3d.png",
        camera="iso",
        html=OUT / f"scene1_{tag}_wave3d.html",
        zoom=1.3,
    )

    # ---- поперечная локализация: объём |E|
    p = render.new_plotter(
        f"{case['title']}: локализация у границы\n"
        f"объём |E|; спад {sol['delta_d_um'] * 1000:.0f} нм в {case['diel']} "
        f"и {sol['delta_m_um'] * 1000:.1f} нм в {case['metal']}, отношение {sol['delta_d_um'] / sol['delta_m_um']:.0f}"
    )
    render.add_field_volume(p, grid, scalars="absE", cmap=render.SEQUENTIAL, title="|E|")
    render.show_axes(p)
    render.finish(p, OUT / f"scene1_{tag}_confinement3d.png", camera="iso", zoom=1.25)

    # ---- дальняя сцена: спад вдоль трассы, ось x сжата так, чтобы трасса была
    # примерно втрое длиннее сечения; срез читается лучше объёма, поскольку
    # огибающая почти всюду мала и объёмный рендер даёт тёмный кадр
    grid3 = render.build_grid(far, y_half_um=y_half, ny=8)
    z_span = float(far.z.max() - far.z.min())
    x_span = float(far.x.max() - far.x.min())
    xscale = 3.0 * z_span / x_span
    p = render.new_plotter(
        f"{case['title']}: затухание вдоль трассы\n"
        f"цвет - огибающая |E| на длине {far.x.max():.0f} мкм = 3 L, "
        f"L = {sol['l_prop_um']:.1f} мкм\n"
        f"ось x сжата в {1.0 / xscale:.0f} раз: трасса {x_span:.0f} мкм против сечения "
        f"{z_span:.1f} мкм"
    )
    render.add_cut_plane(p, grid3, "absE", y_um=0.0, cmap=render.SEQUENTIAL,
                         symmetric=False, title="|E|")
    render.add_layer_boxes(p, [("металл", far.z.min(), 0.0, "dimgray")],
                          (far.x.min(), far.x.max()), y_half, opacity=0.55)
    render.show_axes(p)
    render.finish(p, OUT / f"scene1_{tag}_decay3d.png", camera="iso",
                  scale=(xscale, 1.0, 1.0), zoom=1.1)

    # ---- анимация бегущей волны
    def frame(phase: float):
        g = render.build_grid(near, y_half_um=y_half, ny=10, phase_rad=phase)
        pp = render.new_plotter(
            f"{case['title']}: бегущая волна\nцвет - мгновенная E_z, фаза меняется от 0 до 2 pi",
            window_size=(1100, 700),
        )
        render.add_cut_plane(pp, g, "Ez", y_um=0.0, title="Ez")
        render.add_layer_boxes(pp, [("металл", near.z.min(), 0.0, "dimgray")], x_rng, y_half, 0.4)
        render.show_axes(pp)
        pp.set_scale(1.0, 1.0, 1.0)
        pp.view_isometric()
        pp.camera.zoom(1.3)
        return pp

    render.animate(frame, OUT / f"scene1_{tag}_wave.gif", frames=20, fps=10)
    grid.save(str(OUT / f"scene1_{tag}_field.vti"))


def main() -> int:
    lines: list[str] = []
    add = lines.append
    add("Сцена 1. Плазмон одной границы металл-диэлектрик")
    add("")
    failures = 0

    for tag, case in CASES.items():
        sol = solve_case(case)
        near, far = make_sections(case, sol)
        neff = sol["neff"]

        add(f"--- {case['title']}")
        add(f"  n_eff                       {neff.real:.9f} + {neff.imag:.6e}i")
        add(f"  замкнутая формула           {sol['exact'].real:.9f} + {sol['exact'].imag:.6e}i")
        add(f"  длина волны плазмона        {sol['lambda_spp_um'] * 1000:.1f} нм "
            f"против {case['lambda_um'] * 1000:.0f} нм в свободном пространстве")
        add(f"  спад в диэлектрик           {sol['delta_d_um'] * 1000:.1f} нм")
        add(f"  спад в металл               {sol['delta_m_um'] * 1000:.2f} нм")
        add(f"  отношение спадов            {sol['delta_d_um'] / sol['delta_m_um']:.0f}")
        add(f"  длина распространения       {sol['l_prop_um']:.2f} мкм = "
            f"{sol['l_prop_um'] / case['lambda_um']:.0f} длин волн")
        add(f"  |E_z / E_x| в диэлектрике   {sol['ez_ex_ratio']:.2f}")

        # --- проверки
        d_neff = abs(neff - sol["exact"]) / abs(sol["exact"])
        ok1 = d_neff < 1e-10
        add(f"  [{'OK' if ok1 else 'СБОЙ'}] корень против замкнутой формулы: отклонение {d_neff:.2e}")
        failures += 0 if ok1 else 1

        def build(refine: int, _case=case, _sol=sol) -> FieldSection:
            nz = 1000 * refine + 1
            nx = 8 * refine + 1
            z_local = np.linspace(-5.0 * _sol["delta_m_um"], 3.0 * _sol["delta_d_um"], nz)
            z_local = z_local + 0.5 * (z_local[1] - z_local[0]) - z_local[np.argmin(np.abs(z_local))]
            x_local = np.linspace(0.0, 0.2 * _case["lambda_um"], nx)
            return FieldSection.from_mode(_sol["stack"], _sol["neff"], _case["lambda_um"], x_local, z_local)

        # критерий - именно порядок. Абсолютный уровень невязки задаётся шагом
        # сетки и потому произволен; второй порядок означает, что не осталось
        # ничего, кроме ошибки центральной разности
        conv = residual_convergence(build)
        ok2 = abs(conv["observed_order"] - 2.0) < 0.2
        add(f"  [{'OK' if ok2 else 'СБОЙ'}] невязка Максвелла падает как шаг в квадрате: "
            f"{conv['residual_coarse']:.2e} -> {conv['residual_fine']:.2e} при удвоении сетки, "
            f"порядок {conv['observed_order']:.2f}")

        res = near.maxwell_residual()
        add(f"       на сетке отрисовки невязка {max(res['faraday'], res['ampere_x']):.1e} "
            f"по {res['points_checked']} точкам: сетка картинки грубее расчётной, "
            f"её шаг задан наглядностью")
        failures += 0 if ok2 else 1

        dec = far.energy_decay_check(neff)
        ok3 = dec["rel_error"] < 1e-6
        add(f"  [{'OK' if ok3 else 'СБОЙ'}] спад мощности вдоль x против 2 Im(beta): "
            f"отклонение {dec['rel_error']:.2e}")
        failures += 0 if ok3 else 1

        ratio_num = float(np.abs(near.ez[0, -1] / near.ex[0, -1]))
        ok4 = abs(ratio_num - sol["ez_ex_ratio"]) / sol["ez_ex_ratio"] < 1e-9
        add(f"  [{'OK' if ok4 else 'СБОЙ'}] |E_z/E_x| поля против beta/kappa: "
            f"{ratio_num:.6f} и {sol['ez_ex_ratio']:.6f}")
        failures += 0 if ok4 else 1
        add("")

        plot_profiles(case, sol, near, OUT / f"scene1_{tag}_profiles.png")
        render_3d(case, sol, near, far, tag)

    add("Физический итог")
    add("  Плазмон существует потому, что у металла Re(eps) < 0: только тогда")
    add("  поперечные постоянные по обе стороны границы одновременно вещественны")
    add("  и поле остаётся связанным. Плата за связанность - поглощение, поэтому")
    add("  сильная локализация и большая длина распространения несовместимы:")
    add("  у Ag/воздух 633 нм поле прижато к границе на сотни нанометров, но")
    add("  живёт десятки микрометров, у Au/полимер 1550 нм - наоборот.")

    (OUT / "scene1_single_spp.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
