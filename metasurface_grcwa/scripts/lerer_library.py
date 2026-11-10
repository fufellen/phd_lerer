# -*- coding: utf-8 -*-
"""Материальная библиотека программы А. М. Лерера (письмо 10.11.2026, архив newLIb+mat.rar).

Таблицы *.c содержат Re eps и Im eps на сетке 1 нм; функция программы возвращает Im со знаком минус (конвенция Лерера).
Здесь они переписаны в data/lerer_lib/<имя>.txt (lambda_nm Re Im, Im > 0) и читаются в оптической конвенции Im eps >= 0.
За пределами таблицы берётся крайнее значение; диапазон каждой таблицы - в RANGE.

  python lerer_library.py <папка newLib>   переписать таблицы из исходников
"""
import re
import sys
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parent.parent / "data" / "lerer_lib"
FILES = {"AZO_2": "10_AZO.c", "AZO_wt_0_8": "11_AZO_wt_0_8.c", "AZO_wt_2_0": "12_AZO_wt_2_0.c", "AZO_wt_3_2": "13_AZO_wt_3_2.c", "GZO": "14_GZO.c",
         "HfN": "15_HfN.c", "ITO": "16_ITO.c", "TaN": "17_TaN.c", "TiN": "18_TiN.c", "ZrN": "19_ZrN.c"}
_ROW = re.compile(r"Lam\[\d+\]\s*=\s*([-\d.eE+]+);\s*reEps\[\d+\]\s*=\s*([-\d.eE+]+);\s*imEps\[\d+\]\s*=\s*([-\d.eE+]+);")
_CACHE = {}


def convert(src):
    DATA.mkdir(parents=True, exist_ok=True)
    for name, fn in FILES.items():
        text = (Path(src) / fn).read_text(encoding="cp1251", errors="ignore")
        rows = sorted({float(a): (float(b), float(c)) for a, b, c in _ROW.findall(text)}.items())
        with open(DATA / (name + ".txt"), "w", encoding="utf-8") as f:
            f.write("# %s из библиотеки программы А. М. Лерера (%s, письмо 10.11.2026): lambda_nm Re_eps Im_eps (Im > 0)\n" % (name, fn))
            for lam, (re_, im_) in rows:
                f.write("%g %g %g\n" % (lam, re_, abs(im_)))
        print(name, len(rows), rows[0][0], rows[-1][0])


def table(name):
    if name not in _CACHE:
        a = np.loadtxt(DATA / (name + ".txt"), comments="#")
        _CACHE[name] = a
    return _CACHE[name]


def RANGE(name):
    a = table(name)
    return float(a[0, 0]), float(a[-1, 0])


def eps(name, lam):
    a = table(name)
    return complex(float(np.interp(lam, a[:, 0], a[:, 1])), float(np.interp(lam, a[:, 0], a[:, 2])))


if __name__ == "__main__":
    convert(sys.argv[1])
