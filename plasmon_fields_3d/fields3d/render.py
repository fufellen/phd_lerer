"""Сборка трёхмерных сцен из плоского сечения поля.

Сечение (x, z) вытягивается вдоль y, потому что структура вдоль y однородна.
Отрезок по y выбирается только ради наглядности и всегда подписывается: это
кусок бесконечно широкой структуры, а не волновод конечной ширины.

Масштабы. У плазмонных задач разброс размеров огромен: плёнка золота 10-20 нм
против сотен микрометров распространения. Данные хранятся в микрометрах без
искажения, а непропорциональность вводится только в камере через set_scale, и
коэффициент выводится в подпись. Иначе картинка врёт о геометрии.

Шрифты. Кириллицу отрисовывает только текстовый слой, и лишь при явно заданном
файле шрифта: собственный шрифт VTK её не содержит и молча выбрасывает символы.
Подписи осей рисует CubeAxesActor, который файл шрифта в этой сборке VTK не
принимает, поэтому названия осей заданы латиницей. Это проверено, а не
предположено: с кириллическим названием ось выводится как "x,".
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pyvista as pv
import vtk

from .section import FieldSection

pv.OFF_SCREEN = True

FONT = "C:/Windows/Fonts/arial.ttf" if os.path.exists("C:/Windows/Fonts/arial.ttf") else None

DIVERGING = "RdBu_r"
SEQUENTIAL = "inferno"

AXES = {"xtitle": "x, um", "ytitle": "y, um", "ztitle": "z, um"}


def build_grid(
    section: FieldSection,
    y_half_um: float = 1.0,
    ny: int = 2,
    phase_rad: float = 0.0,
    y_center_um: float = 0.0,
) -> pv.ImageData:
    """Трёхмерная сетка с мгновенными и усреднёнными величинами.

    Точечные массивы:
      Ez, Ex, Hy   - мгновенные вещественные поля при заданной фазе;
      absE         - |E| по модулю фазора, огибающая;
      Sx           - продольная плотность потока мощности;
      eps_re       - действительная часть проницаемости, для показа геометрии;
      E            - вектор мгновенного поля (E_x, 0, E_z).
    """
    ex_t, ez_t, hy_t = section.instantaneous(phase_rad)
    abs_e = section.magnitude_e()
    sx = section.poynting_x()
    eps_re = np.real(section.eps)

    nx, nz = section.x.size, section.z.size
    ny = max(int(ny), 2)
    y = np.linspace(y_center_um - y_half_um, y_center_um + y_half_um, ny)

    grid = pv.ImageData(
        dimensions=(nx, ny, nz),
        spacing=(
            float(section.x[1] - section.x[0]),
            float(y[1] - y[0]),
            float(section.z[1] - section.z[0]),
        ),
        origin=(float(section.x[0]), float(y[0]), float(section.z[0])),
    )

    def spread(a2d: np.ndarray) -> np.ndarray:
        return np.repeat(a2d[:, None, :], ny, axis=1).flatten(order="F")

    grid["Ez"] = spread(ez_t)
    grid["Ex"] = spread(ex_t)
    grid["Hy"] = spread(hy_t)
    grid["absE"] = spread(abs_e)
    grid["Sx"] = spread(sx)
    grid["eps_re"] = spread(eps_re)

    vec = np.zeros((grid.n_points, 3))
    vec[:, 0] = grid["Ex"]
    vec[:, 2] = grid["Ez"]
    grid["E"] = vec
    return grid


def new_plotter(title: str, window_size: tuple[int, int] = (1700, 1050)) -> pv.Plotter:
    p = pv.Plotter(off_screen=True, window_size=window_size)
    p.set_background("white")
    p.add_text(title, font_size=12, color="black", position="upper_left", font_file=FONT)
    return p


def show_axes(plotter: pv.Plotter) -> None:
    plotter.show_grid(color="black", **AXES)


def add_layer_boxes(
    plotter: pv.Plotter,
    layers: list[tuple[str, float, float, str]],
    x_range: tuple[float, float],
    y_half_um: float,
    opacity: float = 0.3,
    y_center_um: float = 0.0,
) -> None:
    """Полупрозрачные плиты материалов: подпись, z_min, z_max, цвет."""
    for name, z0, z1, color in layers:
        box = pv.Box(
            bounds=(x_range[0], x_range[1],
                    y_center_um - y_half_um, y_center_um + y_half_um, z0, z1)
        )
        plotter.add_mesh(box, color=color, opacity=opacity, lighting=False, label=name)


def add_cut_plane(
    plotter: pv.Plotter,
    grid: pv.ImageData,
    scalars: str = "Ez",
    y_um: float = 0.0,
    cmap: str = DIVERGING,
    symmetric: bool = True,
    title: str | None = None,
    opacity: float = 1.0,
    clim_frac: float = 1.0,
) -> None:
    """Срез y = const, раскрашенный по полю. Основной носитель информации.

    `clim_frac` обрезает шкалу до доли от максимума. Это нужно почти всегда:
    у плазмона |E| у самой границы на порядок больше, чем в объёме диэлектрика,
    и при полной шкале вся интересная область выглядит белой.
    """
    sl = grid.slice(normal="y", origin=(0.0, y_um, 0.0))
    if sl.n_points == 0:
        return
    data = sl[scalars]
    peak = float(np.abs(data).max()) * float(clim_frac)
    clim = (-peak, peak) if symmetric else (float(data.min()), float(data.max()))
    plotter.add_mesh(
        sl,
        scalars=scalars,
        cmap=cmap,
        clim=clim,
        opacity=opacity,
        show_scalar_bar=True,
        scalar_bar_args={"title": title or scalars, "vertical": True, "title_font_size": 14},
    )


def add_field_volume(
    plotter: pv.Plotter,
    grid: pv.ImageData,
    scalars: str = "absE",
    clim: tuple[float, float] | None = None,
    cmap: str = SEQUENTIAL,
    opacity: str | list = "linear",
    title: str | None = None,
) -> None:
    plotter.add_volume(
        grid,
        scalars=scalars,
        cmap=cmap,
        clim=clim,
        opacity=opacity,
        scalar_bar_args={"title": title or scalars, "vertical": True, "title_font_size": 14},
    )


def add_signed_isosurfaces(
    plotter: pv.Plotter,
    grid: pv.ImageData,
    scalars: str = "Ez",
    levels: tuple[float, ...] = (0.45,),
    cmap: str = DIVERGING,
    opacity: float = 0.5,
) -> float:
    """Пары изоповерхностей +-level*max: гребни и впадины бегущей волны."""
    data = grid[scalars]
    peak = float(np.abs(data).max())
    if peak <= 0:
        return 0.0
    values = sorted({s * f * peak for f in levels for s in (-1.0, 1.0)})
    surf = grid.contour(values, scalars=scalars)
    if surf.n_points:
        plotter.add_mesh(
            surf,
            scalars=scalars,
            cmap=cmap,
            clim=(-peak, peak),
            opacity=opacity,
            smooth_shading=True,
            show_scalar_bar=False,
        )
    return peak


def add_direction_arrows(
    plotter: pv.Plotter,
    section: FieldSection,
    phase_rad: float = 0.0,
    every_x: int = 16,
    every_z: int = 10,
    y_um: float = 0.0,
    length_um: float = 0.25,
    z_range: tuple[float, float] | None = None,
    cmap: str = "Greys",
) -> None:
    """Стрелки направления мгновенного вектора E в плоскости y = const.

    Длина у всех стрелок одинаковая, а модуль передан цветом. Так сделано
    намеренно: у плазмона |E| падает на порядки за десятки нанометров, и стрелки,
    масштабированные по модулю, вырождаются в иглы у границы и в невидимые точки
    всюду остальном. Направление же информативно везде - именно по нему видно,
    что вектор E почти лежит вдоль распространения и заворачивается в петли.
    """
    ex_t, ez_t, _ = section.instantaneous(phase_rad)
    zmask = np.ones(section.z.size, dtype=bool)
    if z_range is not None:
        zmask = (section.z >= z_range[0]) & (section.z <= z_range[1])

    xs = section.x[::every_x]
    zs = section.z[zmask][::every_z]
    ex_s = ex_t[::every_x][:, zmask][:, ::every_z]
    ez_s = ez_t[::every_x][:, zmask][:, ::every_z]

    xx, zz = np.meshgrid(xs, zs, indexing="ij")
    pts = np.column_stack([xx.ravel(), np.full(xx.size, y_um), zz.ravel()])
    vec = np.zeros((pts.shape[0], 3))
    vec[:, 0] = ex_s.ravel()
    vec[:, 2] = ez_s.ravel()
    mag = np.linalg.norm(vec, axis=1)
    good = mag > 0
    if not np.any(good):
        return
    unit = np.zeros_like(vec)
    unit[good] = vec[good] / mag[good, None]

    cloud = pv.PolyData(pts[good])
    cloud["dir"] = unit[good]
    cloud["mag"] = mag[good] / mag[good].max()
    glyphs = cloud.glyph(orient="dir", scale=False, factor=length_um, geom=pv.Arrow())
    plotter.add_mesh(glyphs, scalars="mag", cmap=cmap, clim=(0.0, 1.0), show_scalar_bar=False)


def frame_isometric(
    plotter: pv.Plotter,
    scale: tuple[float, float, float] = (1.0, 1.0, 1.0),
    direction=(1.0, -1.2, 0.55),
    padding: float = 1.25,
) -> None:
    """Ставит камеру по фактическим границам сцены, а не через view_isometric.

    Зачем понадобилось. `set_scale` растягивает актёры, но камеру под новые
    размеры не подстраивает, и `view_isometric` после него кадр не исправляет:
    у вытянутой сцены объект занимает около процента кадра. Проверено прямым
    измерением доли непустых пикселей - 0.012 против 0.92 у правильного кадра.
    Поэтому положение камеры считается здесь явно, из границ сцены после
    масштабирования, а не запрашивается у готового пресета.
    """
    # plotter.bounds отдаёт границы ДО применения set_scale, поэтому масштаб
    # приходится учитывать здесь вручную. Именно на этом ломался кадр: радиус
    # считался по нерастянутым сотням микрометров, а актёр к тому моменту был
    # сжат в полсотни раз.
    b = plotter.bounds
    sx, sy, sz = scale
    center = np.array([
        sx * (b[0] + b[1]) / 2,
        sy * (b[2] + b[3]) / 2,
        sz * (b[4] + b[5]) / 2,
    ])
    extents = np.array([sx * (b[1] - b[0]), sy * (b[3] - b[2]), sz * (b[5] - b[4])])
    radius = 0.5 * float(np.linalg.norm(extents))
    if radius <= 0:
        return
    view_angle = float(plotter.camera.view_angle)
    distance = padding * radius / np.tan(np.radians(0.5 * view_angle))
    d = np.array(direction, dtype=float)
    d /= np.linalg.norm(d)
    plotter.camera_position = [tuple(center + d * distance), tuple(center), (0.0, 0.0, 1.0)]


def finish(
    plotter: pv.Plotter,
    out_png: Path,
    camera: str | tuple = "iso",
    scale: tuple[float, float, float] = (1.0, 1.0, 1.0),
    window_size: tuple[int, int] = (1700, 1050),
    html: Path | None = None,
    zoom: float = 1.0,
) -> None:
    """Ставит камеру, пишет PNG и, если нужно, интерактивный HTML."""
    plotter.set_scale(xscale=scale[0], yscale=scale[1], zscale=scale[2])
    if camera == "iso":
        frame_isometric(plotter, scale=scale)
    elif camera == "side":
        plotter.camera.tight(padding=0.12, view="xz")
    elif camera == "xz":
        plotter.view_xz(negative=True)
    elif isinstance(camera, tuple):
        plotter.camera_position = camera
    plotter.camera.zoom(zoom)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plotter.screenshot(str(out_png), window_size=window_size)
    if html is not None:
        try:
            plotter.export_html(str(html))
        except Exception as exc:  # pragma: no cover - зависит от окружения trame
            print(f"    HTML не записан: {type(exc).__name__}: {exc}")
    plotter.close()


def animate(
    build_plotter,
    out_gif: Path,
    frames: int = 24,
    fps: int = 12,
) -> None:
    """Кадры по фазе omega t от 0 до 2 pi, склеенные в GIF.

    `build_plotter(phase_rad)` должна вернуть готовый Plotter; кадр снимается
    и плоттер закрывается. Анимация показывает именно бегущую волну, а не
    вращение камеры: движется поле, геометрия неподвижна.
    """
    import imageio.v2 as imageio

    out_gif.parent.mkdir(parents=True, exist_ok=True)
    images = []
    for i in range(frames):
        phase = 2.0 * np.pi * i / frames
        p = build_plotter(phase)
        img = p.screenshot(None, return_img=True)
        p.close()
        images.append(img)
    imageio.mimsave(str(out_gif), images, fps=fps, loop=0)


def patch_axes_font(actor) -> None:
    """Пытается подсунуть CubeAxesActor файл шрифта. Оставлено для наглядности.

    В текущей сборке VTK не срабатывает, поэтому названия осей задаются
    латиницей; функция сохранена, чтобы попытка была видна в коде, а не
    повторялась заново при каждом возврате к сценам.
    """
    if FONT is None:
        return
    for i in range(3):
        for getter in (actor.GetTitleTextProperty, actor.GetLabelTextProperty):
            tp = getter(i)
            tp.SetFontFamily(vtk.VTK_FONT_FILE)
            tp.SetFontFile(FONT)
