# -*- coding: utf-8 -*-
"""Запуск собранной VIBR_T из Python: одна длина волны - один процесс.

Проницаемость столбика и слоя задаётся напрямую как комплексное eps (Im eps < 0, конвенция программы), поэтому
модель эффективной среды и таблицу материала выбирает вызывающий код, а COMPOSITE и COMPOSITE_3 не участвуют.
Геометрия программы - решётка эллиптических цилиндров; прямоугольных столбиков в этой версии нет.

  from vibr_run import Grating, solve
  r = solve(Grating(period=400, ax=250, ay=250, h_pillar=100, h_layer=2000, h_buffer=325), lam=550, eps_pillar=e, eps_layer=e)
"""
import cmath
import math
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from read_output import read

EXE = Path(__file__).resolve().parent / "build" / "vibr_t.exe"


@dataclass
class Grating:
    period: float = 400.0
    ax: float = 250.0
    ay: float = 250.0
    h_pillar: float = 100.0
    h_layer: float = 2000.0
    h_buffer: float = 325.0
    n_buffer: float = 1.45
    n_sub: float = 1.45


@dataclass
class Basis:
    """Числа базисных функций и гармоник. N_g ограничено сверху, см. max_n_g."""
    n_g: int = 20
    m_z: int = 10
    m_fi: int = 3
    m_r: int = 3


def max_n_g(g):
    """Наибольшее N_g без переполнения: RECUR_ZERO.c считает exp(gamma*H), при gamma*H > 709 результат программы - NaN.

    Осевые гармоники (p, 0) и (0, q) в F_df.c не отсекаются, поэтому gamma доходит до 2*pi*N_g/период.
    """
    return int(709.0 * g.period / (2.0 * math.pi * max(g.h_pillar, g.h_layer, g.h_buffer)))


def _nk(eps):
    """n + ik при зависимости exp(iwt): Im eps = -2nk."""
    s = cmath.sqrt(eps.conjugate())
    if s.real < 0:
        s = -s
    return s.real, s.imag


def input_text(g, b, lam, eps_pillar, eps_layer):
    np_, kp = _nk(eps_pillar)
    nl, kl = _nk(eps_layer)
    return "\n".join([
        "/ lambda / %.6f 1000000 %.6f" % (lam, lam + 0.5),
        "/ dx dy tg_pci / %.6f %.6f 0.0" % (g.period, g.period),
        "/ Nsl N_hol N_ell / 4 1 1",
        "/ A_t B_t A_b B_b N_l fi_0 / %.6f %.6f %.6f %.6f 1 0.00" % (g.ax, g.ay, g.ax, g.ay),
        "/ liner / 0 0 0.0",
        "/ n_l k_l / %.8f %.8f" % (np_, kp),
        "/ XC YC / 0.0 0.0",
        "/ top / 0 0 0.0 1.0 0.0",
        "/ pillar layer / 0 0 0.0 %.6f 1. 0.0" % g.h_pillar,
        "/ composite layer / 0 0 0.0 %.6f %.8f %.8f" % (g.h_layer, nl, kl),
        "/ buffer / 0 0 0.0 %.6f %.6f 0" % (g.h_buffer, g.n_buffer),
        "/ substrate / 0 0 0.0 %.6f 0.0" % g.n_sub,
        "/ Azimuth tet fi / 0.00 0.00",
        "/ AM_z / %d 1000 %d" % (b.m_z, b.m_z),
        "/ AM_fi / %d 1000 %d" % (b.m_fi, b.m_fi),
        "/ AM_r / %d 1000 %d" % (b.m_r, b.m_r),
        "/ AN_g / %d 100000 %d" % (b.n_g, b.n_g),
        "/ N_az N_rad / 50 200",
        "/ f_to / 0.0001",
        "/ Rad_ro / 1.00",
        "", ""])


def run_text(text, timeout=1800):
    """Запуск по готовому тексту входного файла; возвращает блоки вывода (см. read_output.read)."""
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        (d / "in.txt").write_text(text, encoding="ascii")
        r = subprocess.run([str(EXE)], cwd=d, input="1\nin.txt\nout.txt\nn\n", capture_output=True, text=True,
                           encoding="cp1251", errors="replace", timeout=timeout)
        if r.returncode or not (d / "out.txt").exists():
            raise RuntimeError("VIBR_T, код %s: %s" % (r.returncode, r.stdout[-400:]))
        return read(str(d / "out.txt"))


def solve(g, lam, eps_pillar, eps_layer, basis=None, pol="s", timeout=1800):
    """R, T, A на длине волны lam. A = delta P из вывода программы."""
    basis = basis or Basis()
    if basis.n_g > max_n_g(g):
        raise ValueError("N_g = %d больше предела %d для этой геометрии: программа выдаст NaN" % (basis.n_g, max_n_g(g)))
    blocks = run_text(input_text(g, basis, lam, eps_pillar, eps_layer), timeout)
    rows = blocks[-1]["rows"]
    if not rows:
        raise RuntimeError("в выводе программы нет чисел (NaN?), lam = %g" % lam)
    key = min(rows, key=lambda k: (abs(k[0] - lam), k[1] != pol))
    return rows[key]
