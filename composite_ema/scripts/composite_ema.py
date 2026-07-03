"""
Python-порт COMPOSITE()/COMPOSITE_3() из composite_ema/source_c/COMPOSIT.c
(формула Клаузиуса-Моссотти / Maxwell-Garnett для эффективной диэл.
проницаемости диэлектрика с плазмонными наночастицами) и базы данных
Ag()/Au()/Cu() из composite_ema/source_c/EPS_MET_AuAgCu.c.

Источник статьи: Лерер А.М. и др., "Электродинамический анализ оптических
поглощающих метаповерхностей, содержащих диэлектрические слои с плазмонными
наночастицами" (см. Obsidian: PhD/Scientific_study/Scientific_study/
paper_analysis/from_Lerer/от Лерера/Лерер АМ_1.md, формула (1)).

Что воспроизводится: только эффективная диэлектрическая проницаемость
eps_eff(lambda) диэлектрического слоя с ПНЧ (формула (1) статьи), то есть
физика "резонанс Фрёлиха внутри знаменателя Клаузиуса-Моссотти".

Что НЕ воспроизводится: полные коэффициенты R, T, P (рис. 2-5 статьи).
Они считаются другой, не найденной в доступных материалах программой -
решением объёмного интегрального уравнения для всей двухпериодической
решётки (метод из [1] в статье). В шаблоне рядом (папка
"нановолноводы от Лерера/N_W_Plan - real") есть родственный, но другой
(старый, для нановолноводов) солвер - не то же самое, что использовалось
для этой конкретной статьи, и не портируется здесь.

Материальные данные металлов (чистовая версия):
- Au, Ag, Cu берутся из Johnson & Christy 1972 (Phys. Rev. B 6, 4370),
  таблицы n,k с refractiveindex.info (main/{Au,Ag,Cu}/nk/Johnson.yml, CC0,
  сырые .yml лежат в composite_ema/data/). Единый диапазон 187.9-1937 нм,
  49 узлов, один источник на все три металла - без прежней склейки из
  C-кода Лерера и УФ-графта. n,k -> eps=eps'-i*eps'' (Re=n^2-k^2, Im=-2nk).
- индекс материала в libmat() (1-Cu, 2-Ag, 4-Au) - по комментарию в
  COMPOSIT.c; сам диспетчер восстановлен, значения теперь из J&C, а не из
  EPS_MET_AuAgCu.c (тот файл оставлен в репозитории как исторический);
- исправлена off-by-one граница интерполяционного окна (в оригинале
  использовалась объявленная ёмкость массива вместо реального числа
  заполненных точек - см. комментарии в EPS_MET_AuAgCu.c);
- HfN, ZrN, TiN (рис. 5 статьи) добавлены отдельно, не из оригинального C:
  TiN - открытые табличные данные (Pflüger et al. 1984, refractiveindex.info,
  n,k -> eps сконвертированы в ту же конвенцию eps=eps'-i*eps'', что и
  Ag/Au/Cu); HfN и ZrN - оцифрованы вручную (color-tracking по пикселям) из
  графика Fig. 4 в Naik, Kim, Boltasseva, "Oxides and nitrides as alternative
  plasmonic materials in the optical range", Opt. Mater. Express 1(6), 1090
  (2011), препринт arXiv:1108.0993 (открытый). Это НЕ тот образец, что
  использовал Лерер - см. README.md, раздел про точность.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"

# ---------------------------------------------------------------------------
# Материальная база: Johnson & Christy 1972 для Au/Ag/Cu.
# Источник: P. B. Johnson, R. W. Christy, "Optical constants of the noble
# metals", Phys. Rev. B 6, 4370 (1972), таблицы n,k с refractiveindex.info
# (main/{Au,Ag,Cu}/nk/Johnson.yml, лицензия CC0, сырые .yml сохранены рядом
# в composite_ema/data/). Единый диапазон 187.9-1937 нм, 49 узлов на металл.
#
# n,k сконвертированы в eps = (n+ik)^2 и ЗАПИСАНЫ в конвенции eps=eps'-i*eps''
# (Im хранится отрицательным, как в исходном C-коде Лерера):
#     Re(eps) = n^2 - k^2 ,   Im_stored = -2*n*k .
# Это чистовая замена прежних таблиц, снятых с C-кода Лерера неясного
# происхождения (Au/Ag) и с УФ-графтом (Cu); см. историю в README.md.
# Столбец 0 = Re(eps), столбец 1 = Im_stored(eps), длины волн в нанометрах.
# ---------------------------------------------------------------------------

_AG_XI = [
    187.9, 191.6, 195.3, 199.3, 203.3, 207.3, 211.9, 216.4,
    221.4, 226.2, 231.3, 237.1, 242.6, 249.0, 255.1, 261.6,
    268.9, 276.1, 284.4, 292.4, 300.9, 310.7, 320.4, 331.5,
    342.5, 354.2, 367.9, 381.5, 397.4, 413.3, 430.5, 450.9,
    471.4, 495.9, 520.9, 548.6, 582.1, 616.8, 659.5, 704.5,
    756.0, 821.1, 892.0, 984.0, 1088.0, 1216.0, 1393.0, 1610.0,
    1937.0,
]
_AG_RE = [
    -0.3240, -0.3078, -0.3206, -0.3311, -0.3571, -0.3289, -0.3156, -0.2965,
    -0.2385, -0.2187, -0.2030, -0.2303, -0.2089, -0.2132, -0.1715, -0.1013,
    0.0220, 0.2165, 0.3904, 0.5842, 0.8663, 0.8974, 0.5024, -0.6583,
    -1.2846, -2.0036, -2.7407, -3.4720, -4.2824, -5.1731, -6.0598, -7.0580,
    -8.2287, -9.5641, -11.0465, -12.8558, -14.8817, -17.2355, -20.0948, -23.4046,
    -27.4777, -32.7969, -39.8397, -48.8865, -60.7604, -77.9255, -101.9931, -140.4000,
    -198.1888,
]
_AG_IM = [
    -2.5937, -2.7104, -2.8112, -2.9116, -2.9808, -3.0963, -3.1800, -3.2598,
    -3.3550, -3.3869, -3.4739, -3.4995, -3.5828, -3.6392, -3.7054, -3.7449,
    -3.7867, -3.7534, -3.5645, -3.2276, -2.5835, -1.3922, -0.6350, -0.2819,
    -0.3198, -0.2838, -0.2320, -0.1864, -0.2070, -0.2275, -0.1970, -0.2126,
    -0.2869, -0.3093, -0.3324, -0.4303, -0.3858, -0.4982, -0.4483, -0.3870,
    -0.3145, -0.4582, -0.5050, -0.5594, -0.6236, -1.5890, -2.6260, -3.5550,
    -6.7584,
]

_AU_XI = [
    187.9, 191.6, 195.3, 199.3, 203.3, 207.3, 211.9, 216.4,
    221.4, 226.2, 231.3, 237.1, 242.6, 249.0, 255.1, 261.6,
    268.9, 276.1, 284.4, 292.4, 300.9, 310.7, 320.4, 331.5,
    342.5, 354.2, 367.9, 381.5, 397.4, 413.3, 430.5, 450.9,
    471.4, 495.9, 520.9, 548.6, 582.1, 616.8, 659.5, 704.5,
    756.0, 821.1, 892.0, 984.0, 1088.0, 1216.0, 1393.0, 1610.0,
    1937.0,
]
_AU_RE = [
    0.2271, 0.2952, 0.2925, 0.2039, 0.1382, -0.0104, -0.1325, -0.2338,
    -0.3463, -0.4155, -0.5510, -0.6169, -0.7445, -0.8913, -1.0804, -1.2365,
    -1.3464, -1.3665, -1.3323, -1.3068, -1.2274, -1.2425, -1.2308, -1.3553,
    -1.3102, -1.2320, -1.4006, -1.6049, -1.6494, -1.7022, -1.6922, -1.7590,
    -1.7027, -2.2783, -3.9462, -5.8421, -8.1127, -10.6619, -13.6482, -16.8177,
    -20.6102, -25.8113, -32.0407, -40.2741, -51.0496, -66.2185, -90.4265, -125.3505,
    -189.0420,
]
_AU_IM = [
    -3.0413, -3.1759, -3.2857, -3.3277, -3.3968, -3.3904, -3.5100, -3.6062,
    -3.7102, -3.8252, -3.8922, -4.0550, -4.1633, -4.3385, -4.4901, -4.7223,
    -4.9763, -5.2824, -5.4949, -5.5964, -5.7803, -5.7926, -5.8458, -5.5737,
    -5.5382, -5.5980, -5.6092, -5.6444, -5.7389, -5.7174, -5.6492, -5.2826,
    -4.8444, -3.8126, -2.5804, -2.1113, -1.6605, -1.3742, -1.0352, -1.0668,
    -1.2718, -1.6266, -1.9254, -2.7940, -3.8610, -5.7015, -8.1863, -12.5552,
    -25.3552,
]

# Cu: тот же источник Johnson & Christy 1972 (main/Cu/nk/Johnson.yml), что и
# Au/Ag выше, единый диапазон 187.9-1937 нм. Прежняя таблица Cu собиралась из
# двух кусков (УФ-графт J&C + C-код Лерера от 517 нм) со стыком-интерполяцией
# на 495.9->517 нм; теперь весь диапазон снят из одного набора данных.
_CU_XI = [
    187.9, 191.6, 195.3, 199.3, 203.3, 207.3, 211.9, 216.4,
    221.4, 226.2, 231.3, 237.1, 242.6, 249.0, 255.1, 261.6,
    268.9, 276.1, 284.4, 292.4, 300.9, 310.7, 320.4, 331.5,
    342.5, 354.2, 367.9, 381.5, 397.4, 413.3, 430.5, 450.9,
    471.4, 495.9, 520.9, 548.6, 582.1, 616.8, 659.5, 704.5,
    756.0, 821.1, 892.0, 984.0, 1088.0, 1216.0, 1393.0, 1610.0,
    1937.0,
]
_CU_RE = [
    -0.9040, -1.0240, -1.1327, -1.2686, -1.4224, -1.5367, -1.6442, -1.7202,
    -1.7403, -1.7334, -1.6984, -1.6088, -1.4408, -1.3022, -1.0430, -0.8714,
    -0.6797, -0.5777, -0.5642, -0.6503, -0.8590, -1.0850, -1.2747, -1.5204,
    -1.6249, -1.7942, -2.0510, -2.4131, -2.7351, -3.2324, -3.7505, -4.2080,
    -4.6028, -5.0857, -5.4093, -5.6005, -6.8216, -10.1820, -13.9916, -17.6379,
    -21.7046, -26.7648, -33.1798, -41.1268, -51.9555, -67.7496, -88.7347, -123.0768,
    -179.1768,
]
_CU_IM = [
    -2.5136, -2.6372, -2.7936, -2.9263, -3.0690, -3.2300, -3.4341, -3.6698,
    -3.9256, -4.1725, -4.4083, -4.6131, -4.8213, -4.8854, -4.9096, -4.7686,
    -4.8372, -4.8063, -4.7357, -4.6377, -4.7012, -4.7720, -4.9211, -4.8803,
    -5.0701, -5.2498, -5.3720, -5.4397, -5.5862, -5.6499, -5.7625, -5.9446,
    -6.2075, -6.2562, -6.1549, -5.2571, -3.7856, -1.9230, -1.6487, -1.7661,
    -2.2392, -2.6936, -3.4608, -4.1094, -5.1962, -7.9152, -11.3268, -16.9024,
    -29.2774,
]

assert len(_AG_XI) == len(_AG_RE) == len(_AG_IM) == 49
assert len(_AU_XI) == len(_AU_RE) == len(_AU_IM) == 49
assert len(_CU_XI) == len(_CU_RE) == len(_CU_IM) == 49

# ---------------------------------------------------------------------------
# TiN: Pflüger et al. 1984 (refractiveindex.info, main/TiN/nk/Pfluger.yml,
# открытые данные). Источник даёт n, k (0.031-2.48 um); здесь сконвертировано
# в eps = (n+ik)^2 = (n^2-k^2) + i*2nk и ЗАТЕМ ЗНАК Im ИНВЕРТИРОВАН, чтобы
# попасть в ту же конвенцию eps=eps'-i*eps'' (Im хранится отрицательным), что
# и у существующих таблиц Ag/Au/Cu выше: Im_stored = -2nk.
# Только 50 узлов на весь диапазон - в видимом/ближнем ИК (0.41-0.83 um)
# всего 4 точки, интерполяция там грубее, чем у Ag/Au/Cu.
# ---------------------------------------------------------------------------

_TIN_XI = [
    31.07, 31.79, 32.63, 33.51, 34.44, 35.42, 36.46, 37.57,
    38.74, 39.99, 41.33, 42.75, 44.28, 45.92, 47.68, 49.59,
    51.66, 53.9, 56.35, 59.04, 61.99, 65.25, 68.88, 72.93,
    77.49, 82.65, 88.56, 95.37, 103.3, 112.7, 124.0, 130.5,
    137.8, 145.9, 155.0, 165.3, 177.1, 190.7, 206.6, 225.4,
    248.0, 275.5, 310.0, 354.2, 413.3, 495.9, 619.9, 826.6,
    1239.9, 2479.7,
]
_TIN_RE = [
    0.703, 0.6863, 0.6902, 0.6856, 0.6592, 0.6361, 0.6041, 0.5702,
    0.5328, 0.4958, 0.4653, 0.4342, 0.393, 0.3473, 0.2851, 0.2188,
    0.1612, 0.1474, 0.1479, 0.1482, 0.1741, 0.2357, 0.2633, 0.2237,
    0.0833, -0.0559, -0.135, -0.1591, -0.2216, -0.6916, -0.9044, -0.8568,
    -0.804, -0.5817, -0.4365, -0.3648, -0.2871, 0.1005, 0.4984, 1.554,
    2.6909, 3.4317, 3.64, 3.456, 1.946, -1.2825, -5.4937, -11.2037,
    -18.1655, -47.6321,
]
_TIN_IM = [
    -0.0848, -0.1079, -0.1219, -0.1168, -0.1152, -0.1184, -0.1198, -0.1312,
    -0.1503, -0.1773, -0.2066, -0.2343, -0.2632, -0.2924, -0.3347, -0.4058,
    -0.5077, -0.6242, -0.727, -0.8241, -0.9318, -0.9944, -1.0186, -1.003,
    -1.0294, -1.2068, -1.3862, -1.6306, -1.669, -1.9258, -2.76, -3.1174,
    -3.5462, -3.8144, -4.2228, -4.6136, -5.084, -5.6108, -6.327, -6.7568,
    -6.93, -6.3756, -5.5842, -4.5368, -3.6192, -3.96, -7.1016, -13.8684,
    -27.1152, -74.58,
]

assert len(_TIN_XI) == len(_TIN_RE) == len(_TIN_IM)

# ---------------------------------------------------------------------------
# HfN, ZrN: оцифровано вручную (color-tracking по RGB пикселей, см.
# CODEX-чекпоинт заметки к статье в Obsidian) из Fig. 4 в Naik, Kim,
# Boltasseva, Opt. Mater. Express 1(6), 1090 (2011), arXiv:1108.0993 -
# конкретный образец ellipsometry на сапфире при 800C, N2:Ar=2:8. НЕ тот
# образец, что использовал Лерер - другая шихта/подложка/температура дадут
# другие числа. Ось Im(eps) на исходном графике положительна (конвенция
# eps=eps'+i*eps''); здесь знак Im инвертирован для согласования с
# остальными таблицами этого файла (eps=eps'-i*eps''). Данные HfN Im
# заметно "срезаны" (clipped) выше ~1800 нм - график не показывает полный
# размах в этой области, только нижнюю границу. Проверка по независимому
# факту статьи ("cross-over wavelengths around 430 nm" для HfN и ZrN):
# оцифровка даёт пересечение нуля ~438 нм (HfN) и ~418 нм (ZrN) - совпадает
# в пределах точности оцифровки графика.
# ---------------------------------------------------------------------------

_HFNZRN_XI = [
    300.0, 320.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0,
    460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0,
    620.0, 640.0, 660.0, 680.0, 700.0, 720.0, 740.0, 760.0,
    780.0, 800.0, 820.0, 840.0, 860.0, 880.0, 900.0, 920.0,
    940.0, 960.0, 980.0, 1000.0, 1020.0, 1040.0, 1060.0, 1080.0,
    1100.0, 1120.0, 1140.0, 1160.0, 1180.0, 1200.0, 1220.0, 1240.0,
    1260.0, 1280.0, 1300.0, 1320.0, 1340.0, 1360.0, 1380.0, 1400.0,
    1420.0, 1440.0, 1460.0, 1480.0, 1500.0, 1520.0, 1540.0, 1560.0,
    1580.0, 1600.0, 1620.0, 1640.0, 1660.0, 1680.0, 1700.0, 1720.0,
    1740.0, 1760.0, 1780.0, 1800.0, 1820.0, 1840.0, 1860.0, 1880.0,
    1900.0, 1920.0, 1940.0, 1960.0, 1980.0, 2000.0,
]
_HFN_RE = [
    3.863, 3.843, 3.381, 3.099, 1.982, 1.497, 0.736, -0.054,
    -0.506, -1.611, -1.887, -3.004, -3.685, -4.373, -5.559, -5.823,
    -7.066, -7.723, -8.5, -9.677, -10.016, -11.277, -11.973, -12.796,
    -14.108, -14.315, -15.627, -16.343, -17.008, -18.271, -18.49, -19.839,
    -20.599, -21.151, -22.463, -22.647, -23.913, -24.513, -25.294, -26.329,
    -26.605, -27.598, -28.031, -28.815, -29.295, -29.92, -30.486, -30.886,
    -31.41, -31.715, -32.399, -32.616, -32.908, -33.182, -33.553, -33.798,
    -33.991, -34.269, -34.269, -34.285, -34.338, -34.379, -34.269, -34.338,
    -34.338, -34.131, -34.18, -34.2, -34.062, -33.924, -33.344, -33.172,
    -33.096, -32.549, -32.312, -31.784, -31.631, -31.327, -30.778, -30.61,
    -30.463, -30.106, -29.851, -31.162, -30.679, -31.678,
]
_HFN_IM = [
    -3.387, -3.266, -2.883, -2.682, -2.471, -2.442, -2.517, -2.552,
    -2.7, -2.816, -3.021, -3.251, -3.619, -3.833, -4.348, -4.626,
    -5.126, -5.668, -5.979, -6.819, -7.002, -7.918, -8.594, -9.062,
    -10.114, -10.767, -11.396, -12.449, -13.478, -13.822, -14.92, -16.201,
    -17.483, -17.62, -18.993, -20.229, -21.648, -22.933, -24.399, -25.767,
    -26.174, -27.506, -28.8, -30.266, -31.716, -33.135, -34.554, -36.04,
    -37.507, -39.085, -40.439, -42.014, -43.478, -44.976, -46.442, -47.918,
    -49.466, -50.933, -52.694, -54.001, -55.639, -57.844, -59.725, -61.104,
    -62.654, -64.22, -65.446, -66.924, -68.513, -69.886, -71.506, -72.54,
    -74.088, -75.539, -76.914, -78.444, -79.214, -79.245,
    # начиная с ~1860 нм оцифровка упирается в верхнюю границу графика
    # (eps''=80) - значения ниже помечены как clipped (см. README.md):
    -79.245, -79.245, -79.245, -79.245, -79.245, -79.245, -79.245, -79.245,
]
_ZRN_RE = [
    3.167, 2.946, 2.627, 1.663, 1.665, 0.115, -0.068, -1.473,
    -1.691, -3.13, -3.268, -4.273, -5.05, -6.48, -6.667, -8.308,
    -8.555, -10.256, -10.671, -12.313, -12.842, -14.384, -15.097, -16.525,
    -17.619, -18.596, -20.668, -20.806, -22.877, -22.977, -24.981, -27.089,
    -27.227, -29.298, -31.187, -31.577, -33.821, -35.103, -36.193, -38.274,
    -38.481, -40.622, -40.76, -43.038, -45.167, -45.386, -47.595, -49.116,
    -50.628, -52.205, -53.425, -54.673, -56.709, -56.907, -59.1, -61.335,
    -61.473, -63.614, -64.127, -65.995, -68.171, -69.259, -70.38, -72.727,
    -72.866, -74.868, -76.421, -77.146, -79.287, -80.887, -81.841, -83.607,
    -83.844, -85.777, -86.005, -88.134, -90.196, -90.334, -92.252, -93.499,
    -94.27, -94.391, -95.283, -97.791, -98.815, -99.073,
]
_ZRN_IM = [
    -2.991, -2.4, -2.509, -2.12, -2.048, -2.037, -2.014, -2.086,
    -2.166, -2.298, -2.426, -2.517, -3.204, -3.06, -3.275, -3.469,
    -3.387, -3.893, -3.982, -4.762, -4.851, -4.917, -5.666, -5.767,
    -6.503, -6.709, -7.26, -7.78, -8.252, -8.97, -9.591, -10.252,
    -11.058, -11.617, -12.694, -13.043, -14.281, -15.11, -15.795, -16.203,
    -17.346, -18.59, -19.357, -20.458, -21.984, -22.105, -23.681, -25.309,
    -26.226, -27.048, -28.696, -29.839, -30.512, -32.128, -33.776, -35.469,
    -35.561, -37.386, -39.039, -40.824, -42.517, -44.302, -44.418, -46.133,
    -47.918, -49.657, -51.396, -53.259, -55.011, -56.904, -58.545, -58.87,
    -60.569, -62.308, -65.111, -66.709, -68.378, -70.151, -72.301, -74.223,
    -74.966, -76.156, -77.757, -77.867, -77.081, -76.917,
]

assert len(_HFNZRN_XI) == len(_HFN_RE) == len(_HFN_IM) == len(_ZRN_RE) == len(_ZRN_IM)


def _local_lagrange(x: float, xi: list[float], column: list[float], mint: int = 2) -> float:
    """Порт локальной интерполяции Лагранжа из Ag()/Au()/Cu() (EPS_MET_AuAgCu.c).

    Отличие от оригинала: граница окна считается по фактическому числу
    заполненных точек (len(xi) - 1), а не по объявленной ёмкости массива
    (в оригинале q/qq были на 1-2 больше реального числа точек - см.
    комментарий в шапке EPS_MET_AuAgCu.c). В пределах табличного диапазона
    результат идентичен оригиналу; отличие может проявляться только у
    самого верхнего края таблицы, которого этот скрипт не достигает.
    """
    n = len(xi)
    i = 0
    while i < n and xi[i] < x:
        i += 1
    if i >= n:
        i = n - 1
    q = n - 1
    b1 = i - mint // 2
    if b1 < 0:
        b1 = 0
    if b1 + mint - 1 > q:
        b1 = q - mint + 1
    ed = b1 + mint - 1
    total = 0.0
    for ii in range(b1, ed + 1):
        term = column[ii]
        for jj in range(b1, ed + 1):
            if jj != ii:
                term *= (x - xi[jj]) / (xi[ii] - xi[jj])
        total += term
    return total


def eps_ag(lam_nm: float) -> complex:
    return complex(_local_lagrange(lam_nm, _AG_XI, _AG_RE), _local_lagrange(lam_nm, _AG_XI, _AG_IM))


def eps_au(lam_nm: float) -> complex:
    return complex(_local_lagrange(lam_nm, _AU_XI, _AU_RE), _local_lagrange(lam_nm, _AU_XI, _AU_IM))


def eps_cu(lam_nm: float) -> complex:
    return complex(_local_lagrange(lam_nm, _CU_XI, _CU_RE), _local_lagrange(lam_nm, _CU_XI, _CU_IM))


def eps_tin(lam_nm: float) -> complex:
    return complex(_local_lagrange(lam_nm, _TIN_XI, _TIN_RE), _local_lagrange(lam_nm, _TIN_XI, _TIN_IM))


def eps_hfn(lam_nm: float) -> complex:
    return complex(_local_lagrange(lam_nm, _HFNZRN_XI, _HFN_RE), _local_lagrange(lam_nm, _HFNZRN_XI, _HFN_IM))


def eps_zrn(lam_nm: float) -> complex:
    return complex(_local_lagrange(lam_nm, _HFNZRN_XI, _ZRN_RE), _local_lagrange(lam_nm, _HFNZRN_XI, _ZRN_IM))


MATERIAL_RANGE_NM = {
    1: (_CU_XI[0], _CU_XI[-1]),
    2: (_AG_XI[0], _AG_XI[-1]),
    4: (_AU_XI[0], _AU_XI[-1]),
    5: (_HFNZRN_XI[0], _HFNZRN_XI[-1]),
    6: (_HFNZRN_XI[0], _HFNZRN_XI[-1]),
    7: (_TIN_XI[0], _TIN_XI[-1]),
}
MATERIAL_NAME = {1: "Cu", 2: "Ag", 4: "Au", 5: "HfN", 6: "ZrN", 7: "TiN"}


def libmat(n_mat: int, lam_nm: float) -> complex:
    """Порт диспетчера libmat(n_mat, f, &Re_eps, &Im_eps) из COMPOSIT.c.

    Индексы 1/2/4 взяты из комментария в COMPOSIT.c: "1- Cu, 2-Ag, 4-Au и
    т. д." (сама libmat() не найдена в доступных исходниках - поведение
    восстановлено по единственной найденной материальной базе с таким же
    стилем, EPS_MET_AuAgCu.c). Индексы 5=HfN, 6=ZrN, 7=TiN - это НЕ часть
    оригинального libmat(): своя нумерация для материалов рис. 5 статьи,
    которых не было в найденном фрагменте C (см. README.md).
    """
    if n_mat == 1:
        return eps_cu(lam_nm)
    if n_mat == 2:
        return eps_ag(lam_nm)
    if n_mat == 4:
        return eps_au(lam_nm)
    if n_mat == 5:
        return eps_hfn(lam_nm)
    if n_mat == 6:
        return eps_zrn(lam_nm)
    if n_mat == 7:
        return eps_tin(lam_nm)
    raise ValueError(
        f"libmat: материал n_mat={n_mat} не поддерживается "
        "(1=Cu, 2=Ag, 4=Au из реконструкции COMPOSIT.c; 5=HfN, 6=ZrN, 7=TiN - расширение этого порта)"
    )


# ---------------------------------------------------------------------------
# Порт COMPOSITE()/COMPOSITE_3() (формула Клаузиуса-Моссотти / Maxwell-Garnett).
# divide_()/mult_() из оригинала заменены обычной комплексной арифметикой Python.
# ---------------------------------------------------------------------------


def composite_from_eps(eps_h: complex, eps_p: complex, C: float) -> complex:
    """Чистая формула (1) статьи, без обращения к материальной базе.

    Используется и в composite(), и в самопроверке self_test() ниже.
    """
    K = (eps_p - eps_h) / (eps_p + 2.0 * eps_h)
    ratio = (1.0 + 2.0 * C * K) / (1.0 - C * K)
    return eps_h * ratio


def composite(eps_h: complex, C: float, n_mat: int, lam_nm: float) -> complex:
    """Порт COMPOSITE(): один сорт наночастиц."""
    eps_p = libmat(n_mat, lam_nm)
    return composite_from_eps(eps_h, eps_p, C)


def composite_3_from_eps(eps_h: complex, populations: list[tuple[complex, float]]) -> complex:
    eta = complex(0.0, 0.0)
    for eps_p, C in populations:
        K = (eps_p - eps_h) / (eps_p + 2.0 * eps_h)
        eta += C * K
    ratio = (1.0 + 2.0 * eta) / (1.0 - eta)
    return eps_h * ratio


def composite_3(eps_h: complex, C1: float, n_mat_1: int, C2: float, n_mat_2: int, lam_nm: float) -> complex:
    """Порт COMPOSITE_3(): два сорта наночастиц."""
    eps_p1 = libmat(n_mat_1, lam_nm)
    eps_p2 = libmat(n_mat_2, lam_nm)
    return composite_3_from_eps(eps_h, [(eps_p1, C1), (eps_p2, C2)])


# ---------------------------------------------------------------------------
# Самопроверка: composite_from_eps() должна совпадать с ручной формулой
# Максвелла-Гарнетта на синтетических (не табличных) eps.
# ---------------------------------------------------------------------------


def self_test() -> None:
    eps_h = complex(2.25, 0.0)  # n=1.5 host, как аналитический тест-кейс
    eps_p = complex(-5.0, 2.0)  # произвольная "металлическая" частица
    for C in (0.0, 0.05, 0.1, 0.3):
        got = composite_from_eps(eps_h, eps_p, C)
        K = (eps_p - eps_h) / (eps_p + 2.0 * eps_h)
        expected = eps_h * (1.0 + 2.0 * C * K) / (1.0 - C * K)
        assert abs(got - expected) < 1e-12, (C, got, expected)
    # C=0 должно точно вернуть матрицу без изменений
    assert composite_from_eps(eps_h, eps_p, 0.0) == eps_h
    # опорные узлы таблиц J&C 1972 (точные табличные длины волн):
    # Au(520.9 nm): n=0.62,k=2.081 -> Re=n^2-k^2=-3.9462, Im=-2nk=-2.5804
    got_au = eps_au(520.9)
    assert abs(got_au - complex(-3.9462, -2.5804)) < 1e-9, got_au
    # Cu(520.9 nm): n=1.18,k=2.608 -> Re=-5.4093, Im=-6.1549
    got_cu = eps_cu(520.9)
    assert abs(got_cu - complex(-5.4093, -6.1549)) < 1e-9, got_cu
    # УФ-узел Cu (n=0.94,k=1.337 at 187.9nm -> Re=-0.904, Im=-2.5136)
    got_cu_uv = eps_cu(187.9)
    assert abs(got_cu_uv - complex(-0.904, -2.5136)) < 1e-9, got_cu_uv
    print("[self_test] OK: composite_from_eps совпадает с аналитической формулой MG;")
    print("[self_test] OK: eps_au(520.9nm) и eps_cu(520.9nm) совпадают с узлами J&C 1972.")

    # HfN/TiN/ZrN: проверяем, что знак Im везде отрицателен (та же конвенция
    # eps=eps'-i*eps'', что у Ag/Au/Cu - иначе резонанс дал бы усиление
    # вместо поглощения) и что найденное пересечение Re[eps]=0 попадает
    # рядом с независимо заявленным в статье Naik et al. значением ~430 нм.
    for name, fn in (("HfN", eps_hfn), ("ZrN", eps_zrn), ("TiN", eps_tin)):
        lo, hi = MATERIAL_RANGE_NM[{"HfN": 5, "ZrN": 6, "TiN": 7}[name]]
        lam = lo
        while lam <= hi:
            im = fn(lam).imag
            assert im <= 1e-6, f"{name}: Im(eps)={im} > 0 at {lam} nm - нарушение конвенции знака (не пассивная среда)"
            lam += (hi - lo) / 40.0
    hfn_cross = find_frohlich_crossing(eps_hfn, 0.0, 300.0, 700.0)  # eps_h=0 -> ищем Re[eps]=0 напрямую
    zrn_cross = find_frohlich_crossing(eps_zrn, 0.0, 300.0, 700.0)
    assert hfn_cross.lam_nm is not None and 380.0 < hfn_cross.lam_nm < 480.0, hfn_cross
    assert zrn_cross.lam_nm is not None and 380.0 < zrn_cross.lam_nm < 480.0, zrn_cross
    print(f"[self_test] OK: Im(eps) <= 0 всюду для HfN/ZrN/TiN (пассивная среда, конвенция знака верна).")
    print(f"[self_test] OK: Re[eps]=0 для HfN на {hfn_cross.lam_nm:.0f} нм, для ZrN на {zrn_cross.lam_nm:.0f} нм "
          f"(Naik et al. 2011 заявляют ~430 нм для обоих - согласуется в пределах точности оцифровки графика).")


@dataclass
class FrohlichCrossing:
    lam_nm: float | None
    note: str


def find_frohlich_crossing(material_eps, eps_h_real: float, lam_lo: float, lam_hi: float, step: float = 1.0) -> FrohlichCrossing:
    """Ищет длину волны, где Re[eps_NP] + 2*eps_h пересекает ноль (условие Фрёлиха)."""
    lam = lam_lo
    prev_lam = lam
    prev_val = material_eps(lam).real + 2.0 * eps_h_real
    lam += step
    while lam <= lam_hi:
        val = material_eps(lam).real + 2.0 * eps_h_real
        if prev_val == 0.0:
            return FrohlichCrossing(prev_lam, "точное совпадение узла")
        if (prev_val < 0.0) != (val < 0.0):
            # линейная интерполяция корня между prev_lam и lam
            t = prev_val / (prev_val - val)
            lam_root = prev_lam + t * (lam - prev_lam)
            return FrohlichCrossing(lam_root, "линейная интерполяция знакопеременного интервала")
        prev_lam, prev_val = lam, val
        lam += step
    return FrohlichCrossing(None, f"в диапазоне {lam_lo}-{lam_hi} нм знакопеременного перехода не найдено")


def run_reproduction() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    n_host = 1.77  # показатель преломления слоя с ПНЧ, из статьи (стр. "Результаты")
    eps_h = complex(n_host ** 2, 0.0)
    C = 0.10  # концентрация ПНЧ 10%, из статьи

    lam_grid = np.arange(250.0, 2000.0 + 1e-9, 2.0)  # от 250 нм - захватывает УФ-данные Au и расширенного Cu

    rows = []
    au_im = []
    cu_im = []
    au_re = []
    cu_re = []
    au_lam_valid = []
    cu_lam_valid = []

    cu_lo, cu_hi = MATERIAL_RANGE_NM[1]
    au_lo, au_hi = MATERIAL_RANGE_NM[4]

    for lam in lam_grid:
        row = {"lambda_nm": lam}
        if au_lo <= lam <= au_hi:
            eps_p_au = eps_au(lam)
            eps_eff_au = composite_from_eps(eps_h, eps_p_au, C)
            row["Re_eps_Au"] = eps_p_au.real
            row["Im_eps_Au"] = eps_p_au.imag
            row["Re_eps_eff_Au_C10"] = eps_eff_au.real
            row["Im_eps_eff_Au_C10"] = eps_eff_au.imag
            au_lam_valid.append(lam)
            au_re.append(eps_eff_au.real)
            au_im.append(eps_eff_au.imag)
        if cu_lo <= lam <= cu_hi:
            eps_p_cu = eps_cu(lam)
            eps_eff_cu = composite_from_eps(eps_h, eps_p_cu, C)
            row["Re_eps_Cu"] = eps_p_cu.real
            row["Im_eps_Cu"] = eps_p_cu.imag
            row["Re_eps_eff_Cu_C10"] = eps_eff_cu.real
            row["Im_eps_eff_Cu_C10"] = eps_eff_cu.imag
            cu_lam_valid.append(lam)
            cu_re.append(eps_eff_cu.real)
            cu_im.append(eps_eff_cu.imag)
        rows.append(row)

    fieldnames = [
        "lambda_nm", "Re_eps_Au", "Im_eps_Au", "Re_eps_eff_Au_C10", "Im_eps_eff_Au_C10",
        "Re_eps_Cu", "Im_eps_Cu", "Re_eps_eff_Cu_C10", "Im_eps_eff_Cu_C10",
    ]
    csv_path = RESULTS_DIR / "composite_eps_Au_vs_Cu.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    print(f"[run] CSV записан: {csv_path}")

    au_cross = find_frohlich_crossing(eps_au, eps_h.real, au_lo, min(au_hi, 1200.0))
    cu_cross = find_frohlich_crossing(eps_cu, eps_h.real, cu_lo, min(cu_hi, 1200.0))
    print(f"[run] Условие Фрёлиха Re[eps_NP]+2*eps_h=0, eps_h={eps_h.real:.4f} (n_host={n_host}):")
    print(f"[run]   Au: lambda = {au_cross.lam_nm} nm ({au_cross.note})")
    print(f"[run]   Cu: lambda = {cu_cross.lam_nm} nm ({cu_cross.note})")

    # Табличная база хранит Im(eps) в конвенции eps=eps'-i*eps'' (все значения
    # отрицательны), поэтому резонансный рост потерь проявляется как минимум
    # (большая по модулю отрицательная величина) Im(eps_eff), а не максимум -
    # ищем экстремум по модулю.
    au_peak_idx = int(np.argmax(np.abs(au_im)))
    cu_peak_idx = int(np.argmax(np.abs(cu_im)))
    print(f"[run] Экстремум |Im(eps_eff)| (Au, C=10%): lambda = {au_lam_valid[au_peak_idx]:.1f} nm, "
          f"Im(eps_eff) = {au_im[au_peak_idx]:.3f}")
    print(f"[run] Экстремум |Im(eps_eff)| (Cu, C=10%): lambda = {cu_lam_valid[cu_peak_idx]:.1f} nm, "
          f"Im(eps_eff) = {cu_im[cu_peak_idx]:.3f}")

    # УФ-сравнение Au/Cu (рис. 2 статьи): теперь возможно - таблица Cu
    # расширена в УФ данными Johnson & Christy 1972 (см. комментарий у
    # _CU_UV_XI выше). Сравниваем |Im(eps_eff)| на общем УФ-участке.
    uv_hi = 400.0
    au_uv = [(lam, im) for lam, im in zip(au_lam_valid, au_im) if lam <= uv_hi]
    cu_uv = [(lam, im) for lam, im in zip(cu_lam_valid, cu_im) if lam <= uv_hi]
    if au_uv and cu_uv:
        au_uv_mean = np.mean([abs(im) for _, im in au_uv])
        cu_uv_mean = np.mean([abs(im) for _, im in cu_uv])
        print(f"[run] УФ (<= {uv_hi:.0f} nm), среднее |Im(eps_eff)|: Au = {au_uv_mean:.3f}, Cu = {cu_uv_mean:.3f}")
        if au_uv_mean > cu_uv_mean:
            print("[run] Au даёт больше потерь (сильнее поглощение) в УФ, чем Cu - "
                  "согласуется с рис. 2 статьи ('ДР с золотыми ПНЧ обеспечивает большее поглощение в УФ диапазоне, чем ДР с медными').")
        else:
            print("[run] РАСХОЖДЕНИЕ с рис. 2 статьи: ожидалось, что Au даёт больше потерь в УФ, чем Cu.")
    else:
        print("[run] ВНИМАНИЕ: недостаточно точек в УФ для сравнения Au/Cu.")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax_re, ax_im) = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
    ax_re.plot(au_lam_valid, au_re, label="Au, C=10%", color="#d4a017")
    ax_re.plot(cu_lam_valid, cu_re, label="Cu, C=10%", color="#b5651d")
    ax_re.axhline(0, color="gray", linewidth=0.5)
    ax_re.set_ylabel(r"$\mathrm{Re}(\varepsilon_\mathrm{eff})$")
    ax_re.legend()
    ax_re.set_title(f"eps_eff(lambda) по COMPOSITE(): n_host={n_host}, C={C:.0%}", fontsize=11)

    ax_im.plot(au_lam_valid, au_im, label="Au, C=10%", color="#d4a017")
    ax_im.plot(cu_lam_valid, cu_im, label="Cu, C=10%", color="#b5651d")
    if au_cross.lam_nm:
        ax_im.axvline(au_cross.lam_nm, color="#d4a017", linestyle="--", linewidth=1,
                       label=f"Au: условие Фрёлиха ~{au_cross.lam_nm:.0f} nm")
    if cu_cross.lam_nm:
        ax_im.axvline(cu_cross.lam_nm, color="#b5651d", linestyle="--", linewidth=1,
                       label=f"Cu: условие Фрёлиха ~{cu_cross.lam_nm:.0f} nm")
    ax_im.set_xlabel("wavelength, nm")
    ax_im.set_ylabel(r"$\mathrm{Im}(\varepsilon_\mathrm{eff})$")
    ax_im.legend(fontsize=8)

    fig.tight_layout()
    png_path = RESULTS_DIR / "composite_eps_Au_vs_Cu.png"
    fig.savefig(png_path, dpi=150)
    print(f"[run] PNG записан: {png_path}")


def run_reproduction_hfn_tin() -> None:
    """Аналог run_reproduction(), но для HfN/TiN (материалы рис. 5 статьи).

    HfN/TiN здесь - оцифровка/открытая таблица (см. заголовок файла и
    README.md), не оригинальные данные Лерера - результат ожидаемо совпадёт
    со статьёй только качественно (форма резонанса, порядок величины), не
    количественно.
    """
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    n_host = 1.77
    eps_h = complex(n_host ** 2, 0.0)
    C = 0.10

    lam_grid = np.arange(300.0, 2000.0 + 1e-9, 2.0)
    hfn_lo, hfn_hi = MATERIAL_RANGE_NM[5]
    tin_lo, tin_hi = MATERIAL_RANGE_NM[7]

    rows = []
    hfn_lam, hfn_re, hfn_im = [], [], []
    tin_lam, tin_re, tin_im = [], [], []

    for lam in lam_grid:
        row = {"lambda_nm": lam}
        if hfn_lo <= lam <= hfn_hi:
            eps_eff = composite_from_eps(eps_h, eps_hfn(lam), C)
            row["Re_eps_eff_HfN_C10"] = eps_eff.real
            row["Im_eps_eff_HfN_C10"] = eps_eff.imag
            hfn_lam.append(lam); hfn_re.append(eps_eff.real); hfn_im.append(eps_eff.imag)
        if tin_lo <= lam <= tin_hi:
            eps_eff = composite_from_eps(eps_h, eps_tin(lam), C)
            row["Re_eps_eff_TiN_C10"] = eps_eff.real
            row["Im_eps_eff_TiN_C10"] = eps_eff.imag
            tin_lam.append(lam); tin_re.append(eps_eff.real); tin_im.append(eps_eff.imag)
        rows.append(row)

    fieldnames = ["lambda_nm", "Re_eps_eff_HfN_C10", "Im_eps_eff_HfN_C10", "Re_eps_eff_TiN_C10", "Im_eps_eff_TiN_C10"]
    csv_path = RESULTS_DIR / "composite_eps_HfN_vs_TiN.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    print(f"[run_hfn_tin] CSV записан: {csv_path}")

    hfn_cross = find_frohlich_crossing(eps_hfn, eps_h.real, hfn_lo, min(hfn_hi, 1200.0))
    tin_cross = find_frohlich_crossing(eps_tin, eps_h.real, tin_lo, min(tin_hi, 1200.0))
    print(f"[run_hfn_tin] Условие Фрёлиха, eps_h={eps_h.real:.4f} (n_host={n_host}):")
    print(f"[run_hfn_tin]   HfN: lambda = {hfn_cross.lam_nm} nm ({hfn_cross.note})")
    print(f"[run_hfn_tin]   TiN: lambda = {tin_cross.lam_nm} nm ({tin_cross.note})")

    hfn_peak = int(np.argmax(np.abs(hfn_im)))
    tin_peak = int(np.argmax(np.abs(tin_im)))
    print(f"[run_hfn_tin] Экстремум |Im(eps_eff)| (HfN, C=10%): lambda = {hfn_lam[hfn_peak]:.1f} nm, "
          f"Im(eps_eff) = {hfn_im[hfn_peak]:.3f}")
    print(f"[run_hfn_tin] Экстремум |Im(eps_eff)| (TiN, C=10%): lambda = {tin_lam[tin_peak]:.1f} nm, "
          f"Im(eps_eff) = {tin_im[tin_peak]:.3f}")
    # Качественная проверка утверждения статьи ("Заключение"/Fig.5): TiN
    # должен быть заметно более широкополосным (более плоский/растянутый
    # резонанс), чем HfN - сравниваем ширину полосы |Im(eps_eff)| > half-max.
    def bandwidth_above_half_max(lam_arr, im_arr, peak_idx):
        half = abs(im_arr[peak_idx]) / 2.0
        above = [lam for lam, im in zip(lam_arr, im_arr) if abs(im) >= half]
        return (min(above), max(above)) if above else (None, None)
    hfn_band = bandwidth_above_half_max(hfn_lam, hfn_im, hfn_peak)
    tin_band = bandwidth_above_half_max(tin_lam, tin_im, tin_peak)
    print(f"[run_hfn_tin] Полоса |Im(eps_eff)|>=half-max: HfN {hfn_band[0]:.0f}-{hfn_band[1]:.0f} nm "
          f"({hfn_band[1]-hfn_band[0]:.0f} nm), TiN {tin_band[0]:.0f}-{tin_band[1]:.0f} nm "
          f"({tin_band[1]-tin_band[0]:.0f} nm)")
    if (tin_band[1] - tin_band[0]) > (hfn_band[1] - hfn_band[0]):
        print("[run_hfn_tin] Качественно согласуется со статьёй: полоса TiN шире полосы HfN.")
    else:
        print("[run_hfn_tin] РАСХОЖДЕНИЕ со статьёй: ожидалась более широкая полоса у TiN, чем у HfN.")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax_re, ax_im) = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
    ax_re.plot(hfn_lam, hfn_re, label="HfN, C=10%", color="#8b3a62")
    ax_re.plot(tin_lam, tin_re, label="TiN, C=10%", color="#2f6f8f")
    ax_re.axhline(0, color="gray", linewidth=0.5)
    ax_re.set_ylabel(r"$\mathrm{Re}(\varepsilon_\mathrm{eff})$")
    ax_re.legend()
    ax_re.set_title(f"eps_eff(lambda) по COMPOSITE(): n_host={n_host}, C={C:.0%} - HfN vs TiN (рис. 5 статьи)", fontsize=10)

    ax_im.plot(hfn_lam, hfn_im, label="HfN, C=10%", color="#8b3a62")
    ax_im.plot(tin_lam, tin_im, label="TiN, C=10%", color="#2f6f8f")
    ax_im.set_xlabel("wavelength, nm")
    ax_im.set_ylabel(r"$\mathrm{Im}(\varepsilon_\mathrm{eff})$")
    ax_im.legend(fontsize=8)

    fig.tight_layout()
    png_path = RESULTS_DIR / "composite_eps_HfN_vs_TiN.png"
    fig.savefig(png_path, dpi=150)
    print(f"[run_hfn_tin] PNG записан: {png_path}")


if __name__ == "__main__":
    self_test()
    run_reproduction()
    run_reproduction_hfn_tin()
