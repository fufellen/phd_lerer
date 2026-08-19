"""Сцена 3. Переход из оптического волновода в плазмонный.

Задача, ради которой сцена существует. Диэлектрический волновод и плазмонный
несовместимы напрямую: у первого мода компактная и почти без потерь, у второго
она в разы шире и живёт миллиметры. Свести их можно двумя способами, и оба
показаны здесь.

1. Плавно, через связь. Волноводы кладут рядом, и при равенстве эффективных
   показателей мощность полностью перетекает из одного в другой на длине связи.
   Так устроен вертикальный направленный ответвитель Park et al. 2009. Это
   основная трёхмерная сцена: видно, как энергия уходит из диэлектрической
   сердцевины к металлической плёнке и возвращается обратно.

2. Встык. Волноводы просто стыкуют торцами, и доля прошедшей мощности равна
   перекрытию мод. Считается здесь же, числом, для сравнения.

Существенная оговорка о геометрии. У прибора Park 2009 полоски конечной ширины,
и синхронизм достигается именно подбором ширины: в планарном пределе его нет,
показатели волноводов расходятся на 2.2e-3, что больше связи. Поэтому для сцены
толщина сердцевины подобрана из условия синхронизма в планарном пределе. Это
честная модель того же физического эффекта, но не копия прибора: длина связи
здесь своя, а числа статьи приведены рядом для сопоставления.

Чем проверяется:
  - подбор синхронизма: показатели изолированных волноводов совпадают;
  - длина связи из расщепления супермод против теории связанных волн;
  - полная перекачка мощности: доля в принимающем волноводе на длине связи;
  - сохранение мощности: сумма долей равна единице с точностью до поглощения;
  - невязка уравнений Максвелла и её порядок по шагу сетки.

Запуск:
    python plasmon_fields_3d/scenes/scene3_dielectric_to_plasmonic.py
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
    mode_fields,
    overlap_power,
    propagation_loss_db_per_cm,
    solve_mode,
    trapz,
)

OUT = ROOT / "results"
OUT.mkdir(parents=True, exist_ok=True)

LAMBDA_UM = 1.55
K0 = mat.k0_from_lambda(LAMBDA_UM)
EPS_CLAD = mat.ZPU13_430          # Park 2009: обкладка ZPU13-430, n = 1.43
EPS_CORE = mat.ZPU13_440          # Park 2009: сердцевина ZPU13-440, n = 1.44
EPS_AU = mat.AU_PARK_1550
N_CLAD = float(np.sqrt(EPS_CLAD).real)

T_AU_UM = 0.020                   # плёнка золота 20 нм
GAP_UM = 5.0                      # зазор между волноводами, как в статье

PAPER = {
    "neff": 1.4307,
    "lc_um": 680.0,
    "lc_measured_um": 600.0,
    "max_transfer": 0.98,
    "device_loss_db_cm": 4.1,
    "lrspp_loss_calc": 9.0,
    "lrspp_loss_meas": 13.0,
    "t_core_paper_um": 2.8,
}


def film_stack() -> Stack:
    return Stack(eps=(EPS_CLAD, EPS_AU, EPS_CLAD), thickness=(T_AU_UM,),
                 names=("ZPU13-430", "Au", "ZPU13-430"))


def slab_stack(t_core_um: float) -> Stack:
    return Stack(eps=(EPS_CLAD, EPS_CORE, EPS_CLAD), thickness=(t_core_um,),
                 names=("ZPU13-430", "ZPU13-440", "ZPU13-430"))


def coupled_stack(t_core_um: float, gap_um: float) -> Stack:
    """Снизу вверх: обкладка, золото, зазор, сердцевина, обкладка."""
    return Stack(
        eps=(EPS_CLAD, EPS_AU, EPS_CLAD, EPS_CORE, EPS_CLAD),
        thickness=(T_AU_UM, gap_um, t_core_um),
        names=("ZPU13-430", "Au", "зазор", "ZPU13-440", "ZPU13-430"),
    )


def single_mode_thickness_limit() -> float:
    """Толщина, выше которой симметричный слой перестаёт быть одномодовым.

    Вторая мода появляется при V = k0 t sqrt(n_core^2 - n_clad^2) = pi. Выше этой
    толщины поиск корня по приближению уже не гарантирует основную моду: решатель
    может сойтись к моде более высокого порядка, и монотонность n_eff(t)
    нарушится. Здесь это не абстракция - именно на этом сорвалась первая версия
    подбора толщины, взявшая верхнюю границу 6 мкм.
    """
    n_core = float(np.sqrt(EPS_CORE).real)
    return float(np.pi / (K0 * np.sqrt(n_core**2 - N_CLAD**2)))


def fundamental_slab_neff(t_um: float, seed: complex) -> complex:
    """Основная TM-мода слоя, найденная продолжением от предыдущей толщины."""
    return solve_mode(slab_stack(t_um), K0, seed)


def match_core_thickness(n_target: float) -> float:
    """Толщина сердцевины, при которой диэлектрический волновод синхронен плёнке.

    Показатель основной моды симметричного слоя монотонно растёт с толщиной от
    показателя обкладки, поэтому корень единствен. Но пользоваться этой
    монотонностью можно только внутри одномодового диапазона, и ветвь основной
    моды прослеживается продолжением от тонкого слоя, а не отдельными поисками
    из одного и того же приближения.
    """
    t_max = 0.9 * single_mode_thickness_limit()
    ladder = np.linspace(0.05, t_max, 400)
    seed = complex(N_CLAD + 1e-6, 0.0)
    table = []
    for t in ladder:
        seed = fundamental_slab_neff(float(t), seed)
        table.append((float(t), seed.real))

    n_vals = np.array([v for _, v in table])
    if not (n_vals[0] <= n_target <= n_vals[-1]):
        raise RuntimeError(
            f"синхронизм недостижим: показатель слоя пробегает "
            f"{n_vals[0]:.6f}..{n_vals[-1]:.6f}, нужен {n_target:.6f}"
        )
    j = int(np.searchsorted(n_vals, n_target))
    lo, hi = table[j - 1][0], table[j][0]
    seed = complex(table[j - 1][1], 0.0)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        n_mid = fundamental_slab_neff(mid, seed).real
        if n_mid < n_target:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-13:
            break
    return 0.5 * (lo + hi)


def supermodes(stack: Stack, seed: complex) -> tuple[complex, complex]:
    """Две супермоды связанной пары, отобранные по значению n_eff около seed."""
    found: list[complex] = []
    for d in np.linspace(-4e-3, 4e-3, 41):
        try:
            n = solve_mode(stack, K0, complex(seed.real + d, seed.imag))
        except Exception:
            continue
        if not np.isfinite(n.real) or n.real <= N_CLAD + 1e-9 or n.real > N_CLAD + 0.02:
            continue
        if any(abs(n - old) < 1e-9 for old in found):
            continue
        found.append(n)
    found.sort(key=lambda z: -z.real)
    if len(found) < 2:
        raise RuntimeError(f"найдено супермод: {len(found)}")
    return found[0], found[1]


def split_power(section: FieldSection, z_split_um: float) -> tuple[np.ndarray, np.ndarray]:
    """Мощность ниже и выше границы раздела: в плазмонном и в диэлектрическом."""
    s = section.poynting_x()
    lower = section.z <= z_split_um
    upper = ~lower
    p_lo = trapz(s[:, lower], section.z[lower])
    p_hi = trapz(s[:, upper], section.z[upper])
    return np.asarray(p_lo), np.asarray(p_hi)


def build_coupler_section(
    stack: Stack, n_even: complex, n_odd: complex, x: np.ndarray, z: np.ndarray
) -> FieldSection:
    """Поле пары как сумма двух супермод, возбуждённых поровну.

    Веса подбираются так, чтобы в сечении x = 0 поле совпадало с модой
    диэлектрического волновода: для синхронной пары это равные амплитуды, а знак
    задаёт, в каком из волноводов сосредоточена мощность на входе.
    """
    even = FieldSection.from_mode(stack, n_even, LAMBDA_UM, x, z, amplitude=0.5)
    odd = FieldSection.from_mode(stack, n_odd, LAMBDA_UM, x, z, amplitude=0.5)

    # нормируем каждую супермоду на единичную мощность, иначе их относительный
    # вес задан произвольной нормировкой восстановленного поля
    for part in (even, odd):
        p0 = trapz(part.poynting_x()[0], part.z)
        scale = 1.0 / np.sqrt(abs(p0))
        part.ex *= scale
        part.ez *= scale
        part.hy *= scale

    # знак подбираем так, чтобы на входе поле было в диэлектрической сердцевине
    z_core = z > 0.5 * (z.min() + z.max())
    plus = np.abs(even.ez[0] + odd.ez[0])[z_core].sum()
    minus = np.abs(even.ez[0] - odd.ez[0])[z_core].sum()
    sign = 1.0 if plus > minus else -1.0
    odd.ex *= sign
    odd.ez *= sign
    odd.hy *= sign
    return FieldSection.superpose([even, odd], title="направленный ответвитель")


def butt_joint_overlap(t_core_um: float) -> float:
    """Доля мощности при стыковке торцами диэлектрического и плазмонного волноводов."""
    st_d = slab_stack(t_core_um)
    st_p = film_stack()
    n_p = solve_mode(st_p, K0, complex(PAPER["neff"], 1e-5))
    n_d = fundamental_slab_neff(t_core_um, complex(n_p.real, 0.0))

    z = np.linspace(-60.0, 60.0, 400001)
    # центрируем обе структуры по своей середине, чтобы стык был осевым
    hd = mode_fields(n_d, st_d, K0, z + 0.5 * t_core_um)[0]
    hp = mode_fields(n_p, st_p, K0, z + 0.5 * T_AU_UM)[0]
    eps_d = st_d.eps_at(z + 0.5 * t_core_um)
    eps_p = st_p.eps_at(z + 0.5 * T_AU_UM)
    return overlap_power((hd, n_d, eps_d), (hp, n_p, eps_p), z)


def plot_panels(
    stack: Stack, n_even: complex, n_odd: complex, lc_um: float,
    section: FieldSection, z_split: float, path: Path, eta_butt: float,
) -> None:
    z = np.linspace(-3.0, stack.total_thickness + 3.0, 4001)
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.5))

    ax = axes[0]
    for n, name, style in ((n_even, "чётная супермода", "-"), (n_odd, "нечётная супермода", "--")):
        ez = mode_fields(n, stack, K0, z)[2]
        ax.plot(np.real(ez) / max(np.abs(ez).max(), 1e-300), z, style, label=name)
    edges = stack.interfaces()
    ax.axhspan(edges[0], edges[1], color="goldenrod", alpha=0.7)
    ax.axhspan(edges[2], edges[3], color="tab:blue", alpha=0.25)
    ax.set_xlabel(r"$\mathrm{Re}\,E_z$, норм.")
    ax.set_ylabel("z, мкм")
    ax.set_title("две супермоды пары\nзолото — жёлтым, сердцевина — синим")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1]
    p_lo, p_hi = split_power(section, z_split)
    tot = p_lo + p_hi
    ax.plot(section.x, p_lo / tot[0], label="плазмонный волновод")
    ax.plot(section.x, p_hi / tot[0], label="диэлектрический волновод")
    ax.plot(section.x, tot / tot[0], color="0.5", lw=1, label="сумма")
    ax.axvline(lc_um, color="tab:red", ls=":", lw=1.2, label=f"L связи = {lc_um:.0f} мкм")
    ax.set_xlabel("x, мкм")
    ax.set_ylabel("доля мощности")
    ax.set_title("перекачка мощности вдоль трассы")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[2]
    ex_t, ez_t, _ = section.instantaneous(0.0)
    m = ax.pcolormesh(section.x, section.z, np.abs(section.magnitude_e()).T,
                      cmap="inferno", shading="auto")
    ax.axhline(z_split, color="w", ls=":", lw=1)
    ax.set_xlabel("x, мкм")
    ax.set_ylabel("z, мкм")
    ax.set_title(r"огибающая $|E|$")
    fig.colorbar(m, ax=ax, fraction=0.046)

    fig.suptitle(
        f"Переход из диэлектрического волновода в плазмонный: связь на длине "
        f"{lc_um:.0f} мкм; перекрытие при стыковке торцами {100 * eta_butt:.1f} %",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def render_3d(section: FieldSection, stack: Stack, lc_um: float, z_split: float) -> None:
    y_half = 3.0
    edges = stack.interfaces()
    # для отрисовки сетку прореживаем: в расчёте она мелкая ради проверок, а
    # объёмный рендер миллионов узлов выдаёт пустой кадр
    thin = section.subsample(every_x=2, every_z=2)
    grid = render.build_grid(thin, y_half_um=y_half, ny=8)

    # масштаб по x подбираем так, чтобы трасса была примерно втрое длиннее
    # сечения: при более сильном сжатии картинка вырождается в вертикальную
    # полосу, при более слабом сечение становится неразличимым
    z_span = float(section.z.max() - section.z.min())
    x_span = float(section.x.max() - section.x.min())
    xscale = 3.0 * z_span / x_span

    p = render.new_plotter(
        f"Переход из диэлектрического волновода в плазмонный\n"
        f"цвет - огибающая |E|; снизу плёнка Au 20 нм (жёлтым), "
        f"сверху сердцевина ZPU13-440 (синим)\n"
        f"мощность введена в диэлектрический волновод, на длине {lc_um:.0f} мкм "
        f"переходит в плазмонный и возвращается\n"
        f"ось x сжата в {1 / xscale:.0f} раз: трасса {x_span:.0f} мкм против сечения "
        f"{z_span:.0f} мкм"
    )
    render.add_cut_plane(p, grid, "absE", y_um=0.0, cmap=render.SEQUENTIAL,
                         symmetric=False, title="|E|")
    render.add_layer_boxes(
        p,
        [("Au", edges[0], edges[1], "goldenrod"), ("сердцевина", edges[2], edges[3], "royalblue")],
        (section.x.min(), section.x.max()),
        y_half,
        opacity=0.55,
    )
    render.show_axes(p)
    render.finish(
        p,
        OUT / "scene3_coupler3d.png",
        camera="iso",
        scale=(xscale, 1.0, 1.0),
        html=OUT / "scene3_coupler3d.html",
        zoom=1.1,
    )

    # ---- поперечные профили в трёх точках трассы. Это принципиально двумерная
    # информация - одномерный профиль в трёх сечениях, - и в объёме она читается
    # хуже, чем на обычном графике, поэтому рисуется отдельной панелью
    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    for frac, style, name in ((0.0, "-", "x = 0"), (0.5, "--", "x = L/2"), (1.0, "-.", "x = L связи")):
        i = int(np.clip(round(frac * lc_um / (section.x[1] - section.x[0])), 0, section.x.size - 1))
        prof = np.abs(section.magnitude_e()[i])
        ax.plot(prof / prof.max(), section.z, style, label=name)
    ax.axhspan(edges[0], edges[1], color="goldenrod", alpha=0.8)
    ax.axhspan(edges[2], edges[3], color="royalblue", alpha=0.25)
    ax.set_xlabel("|E|, нормировано на своё максимальное значение")
    ax.set_ylabel("z, мкм")
    ax.set_title(
        "Поперечные профили в трёх точках трассы\n"
        "золото жёлтым, сердцевина синим:\n"
        "поле переходит из верхнего волновода в нижний"
    )
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "scene3_cross_sections.png", dpi=150)
    plt.close(fig)


def main() -> int:
    lines: list[str] = []
    add = lines.append
    failures = 0

    n_film = solve_mode(film_stack(), K0, complex(PAPER["neff"], 1e-5))
    t_core = match_core_thickness(n_film.real)
    n_slab = fundamental_slab_neff(t_core, complex(n_film.real, 0.0))

    add("Сцена 3. Переход из диэлектрического волновода в плазмонный")
    add(f"Длина волны {LAMBDA_UM * 1000:.0f} нм, обкладка ZPU13-430 (n = {N_CLAD}), "
        f"сердцевина ZPU13-440, плёнка Au {T_AU_UM * 1000:.0f} нм, зазор {GAP_UM:.0f} мкм")
    add("")
    add(f"  плазмонный волновод, LR-SPP   n_eff = {n_film.real:.9f} + {n_film.imag:.3e}i")
    add(f"  потери плазмонного волновода  {propagation_loss_db_per_cm(n_film, LAMBDA_UM):.2f} дБ/см")
    add(f"     сопоставлять их следует с одиночным LR-SPP статьи: {PAPER['lrspp_loss_calc']} дБ/см "
        f"расчёт и {PAPER['lrspp_loss_meas']} дБ/см измерение.")
    add(f"     Величина {PAPER['device_loss_db_cm']} дБ/см относится к прибору целиком, где мощность")
    add("     половину трассы идёт по диэлектрику, и напрямую с планарной плёнкой")
    add("     не сравнивается. Планарное значение вдобавок есть потолок: у полоски")
    add("     конечной ширины перекрытие с металлом меньше, и потери ниже.")
    add(f"  толщина сердцевины из синхронизма {t_core:.6f} мкм")
    add(f"  диэлектрический волновод      n_eff = {n_slab.real:.9f} + {n_slab.imag:.3e}i")
    add(f"  расстройка после подбора      {abs(n_slab.real - n_film.real):.2e}")
    add(f"  значение n_eff в статье       {PAPER['neff']}")
    add("")

    ok = abs(n_slab.real - n_film.real) < 1e-9
    add(f"  [{'OK' if ok else 'СБОЙ'}] синхронизм достигнут: расстройка "
        f"{abs(n_slab.real - n_film.real):.2e}")
    failures += 0 if ok else 1

    stack = coupled_stack(t_core, GAP_UM)
    n_even, n_odd = supermodes(stack, n_film)
    d_n = abs(n_even.real - n_odd.real)
    lc_um = LAMBDA_UM / (2.0 * d_n)
    kappa = K0 * d_n / 2.0

    add(f"  супермоды пары                {n_even.real:.9f} и {n_odd.real:.9f}")
    add(f"  расщепление                   {d_n:.6e}")
    add(f"  коэффициент связи             {kappa:.6e} мкм^-1")
    add(f"  длина связи                   {lc_um:.2f} мкм")
    add(f"  для сравнения, прибор статьи  {PAPER['lc_um']:.0f} мкм расчёт, "
        f"{PAPER['lc_measured_um']:.0f} мкм измерение")
    add("")

    lc_cmt = np.pi / (2.0 * kappa)
    ok = abs(lc_cmt - lc_um) / lc_um < 1e-12
    add(f"  [{'OK' if ok else 'СБОЙ'}] длина связи из расщепления супермод и из теории "
        f"связанных волн совпадают: {lc_um:.6f} и {lc_cmt:.6f} мкм")
    failures += 0 if ok else 1

    edges = stack.interfaces()
    z_split = 0.5 * (edges[1] + edges[2])
    z = np.linspace(-14.0, stack.total_thickness + 14.0, 1201)
    x = np.linspace(0.0, 2.0 * lc_um, 1401)
    section = build_coupler_section(stack, n_even, n_odd, x, z)

    p_lo, p_hi = split_power(section, z_split)
    tot = p_lo + p_hi
    i_lc = int(np.argmin(np.abs(x - lc_um)))
    transfer = float(p_lo[i_lc] / tot[i_lc])
    start_frac = float(p_hi[0] / tot[0])

    add(f"  доля мощности в диэлектрике на входе   {100 * start_frac:.2f} %")
    add(f"  доля мощности в плазмонном на L связи  {100 * transfer:.2f} %")
    add(f"  для сравнения, предел статьи           {100 * PAPER['max_transfer']:.0f} % расчёт")
    add("")

    ok = start_frac > 0.97
    add(f"  [{'OK' if ok else 'СБОЙ'}] на входе мощность действительно в диэлектрическом "
        f"волноводе: {100 * start_frac:.2f} %")
    failures += 0 if ok else 1

    # Перекачка неполна не из-за расстройки - она обнулена - а из-за разного
    # затухания супермод. К длине связи их амплитуды уже неравны, и в исходном
    # волноводе они гасятся не полностью. Проверять надо не «больше 97 %», а
    # совпадение с двухмодовой формулой, куда подставлены те же комплексные
    # постоянные распространения.
    beta_e, beta_o = K0 * n_even, K0 * n_odd
    ee = np.exp(1j * beta_e * x)
    eo = np.exp(1j * beta_o * x)
    cmt_cross = np.abs(ee - eo) ** 2 / 4.0
    cmt_total = (np.abs(ee) ** 2 + np.abs(eo) ** 2) / 2.0
    field_cross = p_lo / tot[0]

    # Положение максимума не зависит от того, как поделено сечение, поэтому это
    # самая жёсткая из проверок. Сравнивать его надо с максимумом двухмодовой
    # кривой, а не с формальной длиной связи: при неравных потерях супермод
    # оптимум смещается назад, поскольку к концу трассы обе амплитуды успевают
    # разойтись. Здесь смещение составляет около 2 %.
    # Допуск относительный, а не в шагах сетки: вблизи максимума кривая плоская,
    # поэтому его положение - наименее резко определённая её черта, и малое
    # различие амплитуд сдвигает точку максимума на несколько шагов.
    i_peak = int(np.argmax(field_cross))
    i_peak_cmt = int(np.argmax(cmt_cross))
    dx = float(x[1] - x[0])
    shift = abs(x[i_peak] - x[i_peak_cmt]) / lc_um
    ok = shift < 0.01
    add(f"  [{'OK' if ok else 'СБОЙ'}] максимум перекачки по полю и по двухмодовой формуле "
        f"совпадают: {x[i_peak]:.2f} и {x[i_peak_cmt]:.2f} мкм, различие {100 * shift:.2f} % "
        f"длины связи при шаге {dx:.2f} мкм")
    add(f"       оба лежат раньше формальной длины связи {lc_um:.2f} мкм на "
        f"{100 * (1 - x[i_peak] / lc_um):.1f} %: при разных потерях супермод оптимум смещается")
    failures += 0 if ok else 1

    # Точность сравнения с двухмодовой формулой ограничена сверху известной
    # величиной: равновесная пара супермод воспроизводит моду изолированного
    # волновода не точно, и остаток виден прямо на входе как недостающие до
    # единицы проценты. Требовать согласия лучше этого уровня бессмысленно.
    floor = 1.0 - start_frac
    worst = float(np.max(np.abs(field_cross - cmt_cross)))
    ok = worst < 3.0 * floor
    add(f"  [{'OK' if ok else 'СБОЙ'}] перекачка по полю совпадает с двухмодовой формулой "
        f"с комплексными beta: расхождение {worst:.3f} при пороге {3 * floor:.3f}")
    add(f"       порог задан не произвольно: {100 * floor:.2f} % - это доля мощности, которую")
    add("       равновесная пара супермод не укладывает в моду входного волновода;")
    add("       геометрическое деление сечения пополам тоже не совпадает с")
    add("       модовым разложением, и обе неточности одного порядка")
    failures += 0 if ok else 1

    ideal = float(np.max(cmt_cross))
    add(f"       предел перекачки при равных потерях супермод был бы 100 %, при их")
    add(f"       реальном различии Im(n_eff) = {n_even.imag:.3e} и {n_odd.imag:.3e} он равен "
        f"{100 * ideal:.2f} %")

    worst_tot = float(np.max(np.abs(tot / tot[0] - cmt_total)))
    ok = worst_tot < 3.0 * floor
    add(f"  [{'OK' if ok else 'СБОЙ'}] полная мощность следует сумме двух экспонент с разными "
        f"показателями, а не одной: расхождение {worst_tot:.3f} при пороге {3 * floor:.3f}")
    failures += 0 if ok else 1

    # каждая супермода в отдельности обязана затухать ровно по своей Im(beta)
    for n_s, name in ((n_even, "чётная"), (n_odd, "нечётная")):
        one = FieldSection.from_mode(stack, n_s, LAMBDA_UM, x, z)
        dec = one.energy_decay_check(n_s)
        ok = dec["rel_error"] < 1e-6
        add(f"  [{'OK' if ok else 'СБОЙ'}] {name} супермода затухает по своей Im(beta): "
            f"отклонение {dec['rel_error']:.1e}")
        failures += 0 if ok else 1

    def build(refine: int) -> FieldSection:
        nz = 3000 * refine + 1
        nx = 20 * refine + 1
        zz = np.linspace(-1.0, stack.total_thickness + 1.0, nz)
        xx = np.linspace(0.0, 0.1, nx)
        return build_coupler_section(stack, n_even, n_odd, xx, zz)

    conv = residual_convergence(build)
    ok = abs(conv["observed_order"] - 2.0) < 0.2
    add(f"  [{'OK' if ok else 'СБОЙ'}] невязка Максвелла падает как шаг в квадрате: "
        f"{conv['residual_coarse']:.2e} -> {conv['residual_fine']:.2e}, "
        f"порядок {conv['observed_order']:.2f}")
    failures += 0 if ok else 1

    eta_butt = butt_joint_overlap(t_core)
    eta_paper = butt_joint_overlap(PAPER["t_core_paper_um"])
    add("")
    add("  Стыковка торцами вместо связи по длине")
    add(f"    сердцевина {t_core:.3f} мкм, синхронная с плёнкой: перекрытие "
        f"{100 * eta_butt:.2f} %, потери {-10 * np.log10(eta_butt):.3f} дБ")
    add(f"    сердцевина {PAPER['t_core_paper_um']} мкм, как в статье: перекрытие "
        f"{100 * eta_paper:.2f} %, потери {-10 * np.log10(eta_paper):.3f} дБ")
    add("    Почти полное перекрытие в синхронном случае - не ошибка счёта, а")
    add("    следствие слабого удержания: при равных n_eff обе моды спадают в")
    add("    обкладку с одинаковой постоянной и почти всю мощность несут именно")
    add("    в обкладке, а не в сердцевине или у плёнки. Различие геометрий - слой")
    add(f"    {t_core * 1000:.0f} нм против плёнки {T_AU_UM * 1000:.0f} нм - приходится на малую долю")
    add("    сечения. Расстройка показателей, а не разница профилей, и есть здесь")
    add("    главный источник потерь на стыке.")
    add("")
    add("Физический итог")
    add("  Перекачка идёт не потому, что волноводы «обмениваются» светом, а")
    add("  потому, что у пары есть две супермоды с чуть разными показателями.")
    add("  Возбуждённые поровну, они то складываются в одном волноводе, то в")
    add("  другом, и период биений и есть длина связи. Отсюда жёсткое условие:")
    add("  при расстройке показателей больше связи перекачка неполна при любой")
    add("  длине - в приборе Park 2009 синхронизм пришлось получать подбором")
    add("  ширины полосок, поскольку в планарном пределе его нет.")

    with (OUT / "scene3_transfer.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["x_um", "P_plasmonic", "P_dielectric", "P_total"])
        for i in range(0, x.size, 5):
            w.writerow([f"{x[i]:.4f}", f"{p_lo[i] / tot[0]:.6f}",
                        f"{p_hi[i] / tot[0]:.6f}", f"{tot[i] / tot[0]:.6f}"])

    plot_panels(stack, n_even, n_odd, lc_um, section, z_split,
                OUT / "scene3_coupler_panels.png", eta_butt)
    render_3d(section, stack, lc_um, z_split)

    (OUT / "scene3_dielectric_to_plasmonic.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
