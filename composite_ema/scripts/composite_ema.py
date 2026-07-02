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

Известные допущения/отличия от оригинального C (см. README.md рядом):
- индекс материала в libmat(): по комментарию в COMPOSIT.c "1-Cu, 2-Ag, 4-Au";
  сама libmat() не найдена, поведение восстановлено по Ag()/Au()/Cu() из
  EPS_MET_AuAgCu.c (см. README про происхождение этого файла);
- исправлена off-by-one граница интерполяционного окна (в оригинале
  использовалась объявленная ёмкость массива вместо реального числа
  заполненных точек - см. комментарии в EPS_MET_AuAgCu.c);
- отброшена "осиротевшая" последняя точка таблицы Cu (fss[84] без xii[84]);
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
# Материальная база (порт Ag()/Au()/Cu() из EPS_MET_AuAgCu.c).
# Столбец 0 = Re(eps), столбец 1 = Im(eps), длины волн в нанометрах.
# ---------------------------------------------------------------------------

_AG_XI = [
    400.0, 450.0, 500.0, 525.0, 550.0, 575.0, 600.0, 625.0, 650.0, 675.0,
    700.0, 750.0, 800.0, 850.0, 900.0, 950.0, 1000.0, 1050.0, 1100.0, 1150.0,
    1200.0, 1250.0, 1300.0, 1325.0, 1350.0, 1375.0, 1395.0, 1405.0, 1410.0,
    1450.0, 1500.0, 1550.0, 1600.0, 1650.0, 1700.0, 1750.0, 1800.0, 1830.0,
    1880.0, 1900.0, 1930.0, 1980.0, 2000.0,
]
_AG_RE = [
    -3.7736, -6.0784, -8.5095, -9.8427, -11.1412, -12.5172, -13.8985,
    -15.3495, -17.0373, -18.7293, -20.4108, -23.9875, -27.9645, -32.2399,
    -36.5736, -40.92, -45.6535, -50.4992, -55.7433, -61.3927, -67.15,
    -67.2816, -73.4795, -65.296, -67.4084, -69.5544, -69.8805, -69.8805,
    -70.048, -76.2027, -80.9776, -86.6415, -92.6895, -98.7424, -104.987,
    -111.422, -118.922, -123.695, -129.76, -132.916, -137.724, -145.345,
    -148.418,
]
_AG_IM = [
    -0.663, -0.741, -0.7592, -0.8164, -0.8016, -0.8496, -0.8952, -1.0192,
    -1.1564, -1.2124, -1.2656, -1.47, -1.4812, -1.704, -2.057, -2.56,
    -2.8392, -3.2706, -3.5856, -4.2336, -4.92, -5.747, -6.3492, -6.3102,
    -6.576, -6.847, -7.0308, -7.0308, -7.0392, -7.5164, -8.109, -8.7608,
    -9.4472, -10.149, -10.8756, -11.627, -12.4488, -12.9412, -13.6884,
    -14.0618, -14.596, -15.5017, -15.86,
]

_AU_XI = [
    250, 275, 295, 300, 315, 335, 350, 375, 400, 415, 435, 445, 450, 465,
    475, 490, 500, 515, 530, 545, 560, 575, 590, 600, 615, 630, 645, 660,
    675, 690, 700, 715, 730, 745, 760, 775, 790, 800, 850, 900, 950, 1000,
    1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550, 1600,
    1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000, 2050, 2100, 2150, 2200,
    2250, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300,
    3400, 3500, 3600, 3700, 3800, 3900, 4000,
]
_AU_RE = [
    -0.5355, -.7329, -0.518, -.4464, -.2625, -.2196, -.36, -.72, -1.086,
    -1.1847, -1.1583, -1.2384, -1.2844, -1.472, -1.6524, -2.3085, -2.8704,
    -4.0131, -5.3865, -6.72, -7.9065, -8.384, -8.694, -8.8776, -9.1887,
    -9.4464, -9.7645, -10.7328, -12.5769, -14.5668, -15.5769, -17.28,
    -18.984, -20.492, -21.9672, -23.5907, -25.168, -26.182, -31.0964,
    -35.9516, -41.4207, -46.4448, -52.05, -58.8924, -65.5076, -71.6184,
    -72.1056, -68.7219, -73.7664, -84.0515, -84.39, -90.1592, -131.745,
    -102.68, -108.34, -115.322, -122.103, -146.474, -136.185, -143.629,
    -150.182, -158.085, -165.678, -203.473, -180.832, -189.699, -198.313,
    -206.887, -223.672, -242.447, -260.958, -281.245, -300.563, -321.246,
    -343.021, -365.773, -389.75, -412.994, -437.698, -463.97, -489.453,
    -514.149, -541.501, -571.855, -598.587,
]
_AU_IM = [
    -4.9468, -6.068, -6.8352, -6.912, -7.0288, -6.696, -6.475, -6.46,
    -6.5072, -6.3896, -6.1056, -5.852, -5.64, -5.0142, -4.536, -3.7332,
    -3.268, -2.646, -2.2752, -1.9912, -1.8112, -1.6878, -1.5392, -1.495,
    -1.3984, -1.232, -1.1268, -1.0496, -1.136, -1.2224, -1.264, -1.3312,
    -1.3952, -1.5402, -1.5946, -1.6524, -1.8072, -1.8432, -2.232, -2.64,
    -2.9624, -3.5464, -4.0432, -4.608, -5.184, -5.929, -6.46, -6.806,
    -7.568, -8.6292, -9.2, -10.0806, -12.6403, -11.9774, -12.9296,
    -14.0076, -15.1231, -17.6533, -17.4939, -18.7836, -20.0422, -21.4232,
    -22.8365, -27.507, -25.7762, -27.3967, -29.0588, -30.7498, -34.2582,
    -38.1404, -42.1693, -46.5443, -51.0478, -55.8785, -60.986, -66.3659,
    -72.1259, -78.0143, -84.2786, -90.8968, -97.6892, -104.609, -112.1,
    -120.073, -127.989,
]

# УФ-хвост (187.9-495.9 нм) взят из Johnson & Christy 1972 (открытые
# данные, refractiveindex.info: main/Cu/nk/Johnson.yml) - в исходном C
# Лерера (EPS_MET.c) таблица Cu начиналась только с 517 нм, УФ-сравнения
# Au/Cu (рис. 2 статьи) без этого не проверить. n,k Джонсона-Кристи
# сконвертированы в ту же конвенцию eps=eps'-i*eps'' формулой
# Re=n^2-k^2, Im=-2nk. Зазор 495.9->517 нм (21 нм, разные первоисточники,
# область не перекрывается) закрывается линейной интерполяцией.
_CU_UV_XI = [
    187.9, 191.6, 195.3, 199.3, 203.3, 207.3, 211.9, 216.4,
    221.4, 226.2, 231.3, 237.1, 242.6, 249.0, 255.1, 261.6,
    268.9, 276.1, 284.4, 292.4, 300.9, 310.7, 320.4, 331.5,
    342.5, 354.2, 367.9, 381.5, 397.4, 413.3, 430.5, 450.9,
    471.4, 495.9,
]
_CU_UV_RE = [
    -0.904, -1.024, -1.1327, -1.2686, -1.4224, -1.5367, -1.6442, -1.7202,
    -1.7403, -1.7334, -1.6984, -1.6088, -1.4408, -1.3022, -1.043, -0.8714,
    -0.6797, -0.5777, -0.5642, -0.6503, -0.859, -1.085, -1.2747, -1.5204,
    -1.6249, -1.7942, -2.051, -2.4131, -2.7351, -3.2324, -3.7505, -4.208,
    -4.6028, -5.0857,
]
_CU_UV_IM = [
    -2.5136, -2.6372, -2.7936, -2.9263, -3.069, -3.23, -3.4341, -3.6698,
    -3.9256, -4.1725, -4.4083, -4.6131, -4.8213, -4.8854, -4.9096, -4.7686,
    -4.8372, -4.8063, -4.7357, -4.6377, -4.7012, -4.772, -4.9211, -4.8803,
    -5.0701, -5.2498, -5.372, -5.4397, -5.5862, -5.6499, -5.7625, -5.9446,
    -6.2075, -6.2562,
]

_CU_XI = _CU_UV_XI + [
    517, 530, 545, 560, 575, 590, 605, 620, 635, 650, 665, 680, 695, 710,
    725, 740, 755, 770, 785, 800, 815, 830, 860, 890, 920, 950, 980, 1010,
    1050, 1090, 1130, 1170, 1210, 1240, 1270, 1300, 1350, 1400, 1450, 1500,
    1530, 1580, 1600, 1630, 1680, 1700, 1730, 1780, 1800, 1830, 1880, 1900,
    1930, 1980, 2000, 2030, 2080, 2100, 2130, 2180, 2200, 2250, 2300, 2350,
    2400, 2450, 2500, 2550, 2600, 2650, 2700, 2750, 2800, 2900, 3000, 3100,
    3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900,
]
_CU_RE = _CU_UV_RE + [
    -5.624, -6.3525, -7.254, -8.1576, -9.1027, -10.0529, -11.0416,
    -12.1091, -13.144, -14.1264, -15.2685, -16.2425, -17.1375, -18.142,
    -19.7152, -20.1552, -21.252, -22.3689, -23.4171, -24.5979, -25.6965,
    -26.9232, -29.3632, -32.0223, -34.6781, -37.4292, -40.296, -43.2684,
    -47.1696, -51.2388, -55.3152, -59.5595, -64.1088, -67.5132, -71.0056,
    -74.4129, -80.3432, -86.5223, -92.9869, -99.803, -104.012, -111.213,
    -114.157, -118.656, -126.276, -129.284, -133.876, -141.696, -144.905,
    -149.762, -158.027, -161.39, -166.513, -175.248, -178.787, -184.647,
    -194.641, -198.696, -204.355, -213.518, -217.263, -227.209, -237.713,
    -248.452, -259.019, -269.84, -280.848, -291.748, -302.82, -314.131,
    -325.432, -336.859, -348.479, -373.096, -399.763, -427.318, -455.782,
    -485.164, -513.306, -541.296, -570.028, -600.223, -631.739, -664.061,
    # оригинальный fss[84]=-697.188 отброшен: нет пары xii[84] (см. README.md)
]
_CU_IM = _CU_UV_IM + [
    -6.1248, -6.05, -5.8752, -5.719, -5.4636, -5.232, -4.896, -4.602,
    -4.1838, -4.256, -3.3012, -3.24, -3.4112, -3.5952, -3.7464, -3.8786,
    -3.9818, -4.18, -4.374, -4.482, -4.6828, -4.8974, -5.2224, -5.5664,
    -5.91, -6.3856, -6.7522, -7.128, -7.579, -8.0416, -8.6536, -9.1332,
    -9.7966, -10.2176, -10.647, -11.072, -11.9361, -12.7966, -13.7106,
    -14.6849, -15.2981, -16.3263, -16.7551, -17.3878, -18.4328, -18.7416,
    -19.1864, -19.9758, -20.2967, -20.7803, -21.5712, -21.9008, -22.4,
    -23.2173, -23.5572, -24.4035, -25.8667, -26.4457, -27.3073, -28.6749,
    -29.221, -30.8505, -32.6079, -34.444, -36.1379, -37.8741, -39.648,
    -41.7141, -43.8269, -46.0269, -48.2615, -50.5401, -52.9055, -57.6928,
    -62.6518, -67.6823, -72.4717, -77.4244, -82.5972, -87.9993, -93.5682,
    -99.5612, -105.932, -112.495,
    # оригинальный fss[84]=-119.25 отброшен: нет пары xii[84] (см. README.md)
]

assert len(_CU_UV_XI) == len(_CU_UV_RE) == len(_CU_UV_IM) == 34
assert len(_CU_XI) == len(_CU_RE) == len(_CU_IM) == 118
assert len(_AU_XI) == len(_AU_RE) == len(_AU_IM) == 85
assert len(_AG_XI) == len(_AG_RE) == len(_AG_IM) == 43

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
    # известная опорная точка таблицы: Au(530 nm) = -5.3865 - 2.2752j
    got_au_530 = eps_au(530.0)
    assert abs(got_au_530 - complex(-5.3865, -2.2752)) < 1e-9, got_au_530
    got_cu_517 = eps_cu(517.0)
    assert abs(got_cu_517 - complex(-5.624, -6.1248)) < 1e-9, got_cu_517
    # УФ-узел Cu (Johnson & Christy 1972, n=0.94,k=1.337 at 187.9nm ->
    # Re=n^2-k^2=-0.904, Im=-2nk=-2.5136)
    got_cu_uv = eps_cu(187.9)
    assert abs(got_cu_uv - complex(-0.904, -2.5136)) < 1e-9, got_cu_uv
    print("[self_test] OK: composite_from_eps совпадает с аналитической формулой MG;")
    print("[self_test] OK: eps_au(530nm) и eps_cu(517nm) совпадают с табличными узлами.")

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
