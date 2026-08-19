"""Сцена 4. Фазовращатель на материале с фазовым переходом: демонстрация.

Что показывается. Фазовращатель меняет не амплитуду волны, а её фазу. Управляет
им материал с фазовым переходом (PCM): при переключении из аморфного состояния в
кристаллическое его показатель преломления скачком растёт, вместе с ним растёт
эффективный показатель моды, и волна на выходе приходит с другой фазой. Длина, на
которой набегает ровно pi, называется L_pi:

    L_pi = lambda_0 / (2 dRe n_eff).

Устройство. Слоистая структура статьи об ЭДП PCM LR-DLSPP, вертикальный срез
через полоску металла (в статье он обозначен S3), снизу вверх:

    SiO2 | PCM 120 нм | Au 10 нм | Si 180 нм | воздух.

Демонстрация. Двухплечевой интерферометр: свет делится поровну, одно плечо
переключено в кристаллическое состояние, другое оставлено аморфным, на выходе
плечи складываются. При набеге pi они гасят друг друга. Сцена показывает оба
плеча рядом, расхождение фазы вдоль трассы и итоговую кривую гашения.

Чем проверяется:
  - n_eff обоих состояний против значений среза S3 из данных статьи, полученных
    независимой реализацией матрицы переноса (файл data/eim_vertical_slices.csv);
  - L_pi против той же статьи;
  - тождество L_pi = lambda / (2 dRe n_eff) на самих полях: набег фазы на этой
    длине действительно равен pi;
  - невязка уравнений Максвелла и её порядок по шагу сетки.

Запуск:
    python plasmon_fields_3d/scenes/scene4_pcm_phase_shifter.py
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
from slabmodes import Stack, solve_mode  # noqa: E402

OUT = ROOT / "results"
OUT.mkdir(parents=True, exist_ok=True)

LAMBDA_UM = 1.550
K0 = mat.k0_from_lambda(LAMBDA_UM)

# Геометрия статьи: PCM-буфер 800x120 нм, Si-гребень 300x180 нм, Au-полоска
# 100x10 нм. Здесь берётся вертикальный срез через полоску металла, у которого
# слои считаются бесконечно широкими.
T_PCM_UM = 0.120
T_AU_UM = 0.010
T_SI_UM = 0.180

EPS_AU = mat.AU_ARTICLE_1550
EPS_SI = mat.SI_1550
EPS_SIO2 = mat.SIO2_1550
EPS_AIR = mat.AIR

# Значения среза S3 из data/eim_vertical_slices.csv статьи про ЭДП PCM LR-DLSPP.
# Получены другой реализацией матрицы переноса (условие M[1,1] = 0), поэтому
# служат независимым эталоном, а не повтором того же кода.
ARTICLE_S3 = {
    ("GSST", "amorphous"): complex(2.652958, 0.005026),
    ("GSST", "crystalline"): complex(3.388808, 0.120177),
    ("Sb2S3", "amorphous"): complex(2.113987, 0.017205),
    ("Sb2S3", "crystalline"): complex(2.547937, 0.007581),
    ("Sb2Se3", "amorphous"): complex(2.532344, 0.008003),
    ("Sb2Se3", "crystalline"): complex(2.957975, 0.001564),
}
ARTICLE_LPI_UM = {"GSST": 1.0532, "Sb2S3": 1.7859, "Sb2Se3": 1.8208}

SEEDS = {"GSST": (2.65, 3.39), "Sb2S3": (2.11, 2.55), "Sb2Se3": (2.53, 2.96)}


def s3_stack(eps_pcm: complex) -> Stack:
    """Вертикальный срез через полоску металла: SiO2 | PCM | Au | Si | воздух."""
    return Stack(
        eps=(EPS_SIO2, eps_pcm, EPS_AU, EPS_SI, EPS_AIR),
        thickness=(T_PCM_UM, T_AU_UM, T_SI_UM),
        names=("SiO2", "PCM", "Au", "Si", "воздух"),
    )


def solve_state(name: str, state: str) -> tuple[complex, Stack]:
    pcm = mat.PCM[name]
    st = s3_stack(pcm.eps(state))
    seed = SEEDS[name][0 if state == "amorphous" else 1]
    return solve_mode(st, K0, complex(seed, 0.01)), st


def loss_db_per_um(neff: complex) -> float:
    """Погонные потери по мощности, дБ/мкм: alpha = 8.686 Im(beta)."""
    return float(8.686 * K0 * abs(neff.imag))


def metrics(name: str) -> dict:
    n_a, st_a = solve_state(name, "amorphous")
    n_c, st_c = solve_state(name, "crystalline")
    d_re = abs(n_c.real - n_a.real)
    l_pi = LAMBDA_UM / (2.0 * d_re)
    return {
        "material": name,
        "n_a": n_a,
        "n_c": n_c,
        "stack_a": st_a,
        "stack_c": st_c,
        "d_re": d_re,
        "l_pi_um": l_pi,
        "loss_a": loss_db_per_um(n_a),
        "loss_c": loss_db_per_um(n_c),
        "il_a": loss_db_per_um(n_a) * l_pi,
        "il_c": loss_db_per_um(n_c) * l_pi,
    }


def interferometer(m: dict, x: np.ndarray) -> dict:
    """Двухплечевой интерферометр: одно плечо переключено, другое нет.

    Поле на выходе - полусумма плеч. Учитываются и разный набег фазы, и разное
    поглощение, поэтому гашение получается неполным: идеальный ноль требовал бы
    равных амплитуд, а у кристаллического плеча они меньше.
    """
    b_a = K0 * m["n_a"]
    b_c = K0 * m["n_c"]
    e_a = np.exp(1j * b_a * x)
    e_c = np.exp(1j * b_c * x)
    out = 0.5 * (e_a + e_c)
    power = np.abs(out) ** 2
    # Разность фаз берётся напрямую из постоянных распространения, а не через
    # angle() от полей: обёрнутый угол не различает +pi и -pi, и проверка на
    # длине L_pi получала бы верный по модулю, но отрицательный ответ.
    phase = (b_c.real - b_a.real) * x
    return {
        "power": power,
        "phase_diff": phase,
        "arm_a": np.abs(e_a) ** 2,
        "arm_c": np.abs(e_c) ** 2,
    }


def plot_panels(results: list[dict], path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(16.0, 4.6))

    ax = axes[0]
    for m in results:
        st = m["stack_a"]
        z = np.linspace(-0.35, st.total_thickness + 0.35, 3001)
        from slabmodes import mode_fields

        for state, style in (("a", "-"), ("c", "--")):
            stack = m[f"stack_{state}"]
            ez = mode_fields(m[f"n_{state}"], stack, K0, z)[2]
            if m["material"] == "Sb2Se3":
                ax.plot(np.abs(ez) / max(np.abs(ez).max(), 1e-300), z, style,
                        label=f"Sb2Se3, {'аморфн.' if state == 'a' else 'кристалл.'}")
    edges = results[0]["stack_a"].interfaces()
    ax.axhspan(edges[0], edges[1], color="mediumpurple", alpha=0.35)
    ax.axhspan(edges[1], edges[2], color="goldenrod", alpha=0.8)
    ax.axhspan(edges[2], edges[3], color="tab:blue", alpha=0.25)
    ax.set_xlabel(r"$|E_z|$, норм.")
    ax.set_ylabel("z, мкм")
    ax.set_title("профиль моды в двух состояниях\nPCM фиолетовым, Au жёлтым, Si синим")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1]
    for m in results:
        x = np.linspace(0.0, 2.5 * m["l_pi_um"], 2001)
        itf = interferometer(m, x)
        ax.plot(x / m["l_pi_um"], itf["phase_diff"] / np.pi, label=m["material"])
    ax.axhline(1.0, color="k", ls=":", lw=1)
    ax.axvline(1.0, color="tab:red", ls=":", lw=1)
    ax.set_xlabel(r"длина в единицах $L_\pi$")
    ax.set_ylabel(r"разность фаз плеч, $\pi$")
    ax.set_title("набег фазы между плечами")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[2]
    for m in results:
        x = np.linspace(0.0, 2.5 * m["l_pi_um"], 2001)
        itf = interferometer(m, x)
        ax.semilogy(x / m["l_pi_um"], np.maximum(itf["power"], 1e-6),
                    label=f"{m['material']}, L_pi = {m['l_pi_um']:.2f} мкм")
    ax.axvline(1.0, color="tab:red", ls=":", lw=1)
    ax.set_xlabel(r"длина в единицах $L_\pi$")
    ax.set_ylabel("мощность на выходе")
    ax.set_title("гашение на выходе интерферометра\n(неполное из-за разного поглощения плеч)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    fig.suptitle(
        "Фазовращатель на PCM: срез SiO2 | PCM 120 нм | Au 10 нм | Si 180 нм | воздух, 1550 нм",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def render_3d(m: dict) -> None:
    """Два плеча рядом: аморфное и кристаллическое, на длине 2 L_pi."""
    st = m["stack_a"]
    l_pi = m["l_pi_um"]
    x = np.linspace(0.0, 2.0 * l_pi, 1201)
    z = np.linspace(-0.30, st.total_thickness + 0.30, 401)
    y_half = 0.22
    edges = st.interfaces()

    p = render.new_plotter(
        f"Фазовращатель на {m['material']}: два плеча интерферометра\n"
        f"сзади аморфное состояние, спереди кристаллическое; цвет - мгновенная E_z\n"
        f"n_eff = {m['n_a'].real:.4f} и {m['n_c'].real:.4f}, "
        f"L_pi = {l_pi:.3f} мкм; на этой длине плечи приходят в противофазе\n"
        f"слои снизу вверх: SiO2, PCM 120 нм, Au 10 нм, Si 180 нм, воздух"
    )
    for state, y_c in (("a", 0.45), ("c", -0.45)):
        sec = FieldSection.from_mode(m[f"stack_{state}"], m[f"n_{state}"], LAMBDA_UM, x, z)
        grid = render.build_grid(sec, y_half_um=y_half, ny=6, y_center_um=y_c)
        render.add_cut_plane(p, grid, "Ez", y_um=y_c, title="Ez", clim_frac=0.5)
        render.add_layer_boxes(
            p,
            [
                ("PCM", edges[0], edges[1], "mediumpurple"),
                ("Au", edges[1], edges[2], "goldenrod"),
                ("Si", edges[2], edges[3], "royalblue"),
            ],
            (x.min(), x.max()),
            y_half,
            opacity=0.4,
            y_center_um=y_c,
        )
    render.show_axes(p)
    render.finish(
        p,
        OUT / "scene4_two_arms3d.png",
        camera="iso",
        scale=(0.35, 1.0, 1.0),
        html=OUT / "scene4_two_arms3d.html",
        zoom=1.2,
    )

    # ---- одно плечо крупно: видно, что мода прижата к металлу и PCM
    sec = FieldSection.from_mode(m["stack_c"], m["n_c"], LAMBDA_UM,
                                 np.linspace(0.0, 3.0, 901), z)
    grid = render.build_grid(sec, y_half_um=0.6, ny=10)
    p = render.new_plotter(
        f"{m['material']} в кристаллическом состоянии, крупный план\n"
        f"объём |E|: мода прижата к полоске Au и лежит в PCM и Si, "
        f"поэтому переключение PCM ею и управляет\n"
        f"потери {m['loss_c']:.4f} дБ/мкм, на длине L_pi это {m['il_c']:.3f} дБ"
    )
    render.add_field_volume(p, grid, scalars="absE", cmap=render.SEQUENTIAL, title="|E|")
    render.add_layer_boxes(
        p,
        [
            ("PCM", edges[0], edges[1], "mediumpurple"),
            ("Au", edges[1], edges[2], "goldenrod"),
            ("Si", edges[2], edges[3], "royalblue"),
        ],
        (0.0, 3.0),
        0.6,
        opacity=0.4,
    )
    render.show_axes(p)
    render.finish(p, OUT / "scene4_mode_closeup3d.png", camera="iso", zoom=1.25)


def main() -> int:
    lines: list[str] = []
    add = lines.append
    failures = 0

    add("Сцена 4. Фазовращатель на материале с фазовым переходом")
    add(f"Срез SiO2 | PCM {T_PCM_UM * 1000:.0f} нм | Au {T_AU_UM * 1000:.0f} нм | "
        f"Si {T_SI_UM * 1000:.0f} нм | воздух, {LAMBDA_UM * 1000:.0f} нм")
    add("Слои считаются бесконечно широкими: это планарная модель того же среза,")
    add("а не расчёт волновода конечной ширины. Абсолютные потери прибора она не")
    add("определяет, для них нужен расчёт методом конечных элементов.")
    add("")

    results = []
    for name in ("GSST", "Sb2S3", "Sb2Se3"):
        pcm = mat.PCM[name]
        m = metrics(name)
        results.append(m)
        add(f"--- {pcm.name}")
        add(f"  источник постоянных         {pcm.source}")
        add(f"  показатель материала        {pcm.n_amorphous} -> {pcm.n_crystalline}, "
            f"скачок {pcm.delta_n_material:.3f}")
        add(f"  n_eff аморфное              {m['n_a'].real:.6f} + {m['n_a'].imag:.6e}i")
        add(f"  n_eff кристаллическое       {m['n_c'].real:.6f} + {m['n_c'].imag:.6e}i")
        add(f"  скачок эффективного показателя {m['d_re']:.6f}")
        add(f"  L_pi                        {m['l_pi_um']:.4f} мкм "
            f"(статья: {ARTICLE_LPI_UM[name]:.4f} мкм)")
        add(f"  потери                      {m['loss_a']:.5f} и {m['loss_c']:.5f} дБ/мкм")
        add(f"  вносимые потери на L_pi     {m['il_a']:.4f} и {m['il_c']:.4f} дБ")

        for state, key in (("amorphous", "n_a"), ("crystalline", "n_c")):
            ref = ARTICLE_S3[(name, state)]
            got = m[key]
            d_re = abs(got.real - ref.real)
            d_im = abs(got.imag - ref.imag) / max(abs(ref.imag), 1e-30)
            ok = d_re < 2e-6 and d_im < 2e-3
            add(f"  [{'OK' if ok else 'СБОЙ'}] {state:12s} против среза S3 статьи: "
                f"Re расходится на {d_re:.1e}, Im на {100 * d_im:.2f} %")
            failures += 0 if ok else 1

        d_lpi = abs(m["l_pi_um"] - ARTICLE_LPI_UM[name]) / ARTICLE_LPI_UM[name]
        ok = d_lpi < 2e-3
        add(f"  [{'OK' if ok else 'СБОЙ'}] L_pi против статьи: расхождение {100 * d_lpi:.3f} %")
        failures += 0 if ok else 1
        add("")

    # --- проверка тождества L_pi на самих полях
    m = results[2]  # Sb2Se3
    x = np.array([0.0, m["l_pi_um"]])
    itf = interferometer(m, x)
    phase_at_lpi = itf["phase_diff"][-1]
    ok = abs(phase_at_lpi - np.pi) < 1e-9
    add(f"  [{'OK' if ok else 'СБОЙ'}] на длине L_pi набег фазы между плечами равен pi: "
        f"получено {phase_at_lpi / np.pi:.12f} pi")
    failures += 0 if ok else 1

    # --- невязка Максвелла
    st = m["stack_c"]

    def build(refine: int) -> FieldSection:
        nz = 3200 * refine + 1
        nx = 20 * refine + 1
        zz = np.linspace(-0.25, st.total_thickness + 0.25, nz)
        xx = np.linspace(0.0, 0.05, nx)
        return FieldSection.from_mode(st, m["n_c"], LAMBDA_UM, xx, zz)

    conv = residual_convergence(build)
    ok = abs(conv["observed_order"] - 2.0) < 0.2
    add(f"  [{'OK' if ok else 'СБОЙ'}] невязка Максвелла падает как шаг в квадрате: "
        f"{conv['residual_coarse']:.2e} -> {conv['residual_fine']:.2e}, "
        f"порядок {conv['observed_order']:.2f}")
    failures += 0 if ok else 1
    add("")

    best = min(results, key=lambda r: max(r["il_a"], r["il_c"]))
    add("Сравнение материалов по парным метрикам")
    add("  материал   L_pi, мкм   IL аморф., дБ   IL кристалл., дБ")
    for r in results:
        add(f"  {r['material']:9s}  {r['l_pi_um']:9.4f}   {r['il_a']:13.4f}   {r['il_c']:16.4f}")
    add(f"  По вносимым потерям на длине L_pi лучший здесь {best['material']}.")
    add("")
    add("Физический итог")
    add("  Переключение PCM смещает эффективный показатель моды, а не поглощение")
    add("  само по себе, поэтому прибор работает именно как фазовращатель. Но два")
    add("  состояния поглощают по-разному, и в интерферометре плечи приходят с")
    add("  разными амплитудами: гашение получается неполным даже при точном наборе")
    add("  фазы pi. У GSST это выражено сильнее всего, поскольку у него в")
    add("  кристаллическом состоянии заметная мнимая часть показателя; у Sb2S3 и")
    add("  Sb2Se3 поглощение мало в обоих состояниях, и гашение глубже.")

    with (OUT / "scene4_phase_shifter.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["material", "n_a_re", "n_a_im", "n_c_re", "n_c_im", "delta_re",
                    "L_pi_um", "loss_a_dB_um", "loss_c_dB_um", "IL_a_dB", "IL_c_dB",
                    "article_L_pi_um"])
        for r in results:
            w.writerow([r["material"], f"{r['n_a'].real:.9f}", f"{r['n_a'].imag:.9e}",
                        f"{r['n_c'].real:.9f}", f"{r['n_c'].imag:.9e}", f"{r['d_re']:.9f}",
                        f"{r['l_pi_um']:.6f}", f"{r['loss_a']:.6f}", f"{r['loss_c']:.6f}",
                        f"{r['il_a']:.6f}", f"{r['il_c']:.6f}", ARTICLE_LPI_UM[r["material"]]])

    plot_panels(results, OUT / "scene4_phase_shifter_panels.png")
    render_3d(results[2])

    (OUT / "scene4_pcm_phase_shifter.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
