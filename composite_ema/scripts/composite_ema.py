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
  HfN/TiN (рис. 5 статьи) в найденных исходниках отсутствуют - не
  воспроизводятся;
- исправлена off-by-one граница интерполяционного окна (в оригинале
  использовалась объявленная ёмкость массива вместо реального числа
  заполненных точек - см. комментарии в EPS_MET_AuAgCu.c);
- отброшена "осиротевшая" последняя точка таблицы Cu (fss[84] без xii[84]).
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

_CU_XI = [
    517, 530, 545, 560, 575, 590, 605, 620, 635, 650, 665, 680, 695, 710,
    725, 740, 755, 770, 785, 800, 815, 830, 860, 890, 920, 950, 980, 1010,
    1050, 1090, 1130, 1170, 1210, 1240, 1270, 1300, 1350, 1400, 1450, 1500,
    1530, 1580, 1600, 1630, 1680, 1700, 1730, 1780, 1800, 1830, 1880, 1900,
    1930, 1980, 2000, 2030, 2080, 2100, 2130, 2180, 2200, 2250, 2300, 2350,
    2400, 2450, 2500, 2550, 2600, 2650, 2700, 2750, 2800, 2900, 3000, 3100,
    3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900,
]
_CU_RE = [
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
_CU_IM = [
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

assert len(_CU_XI) == len(_CU_RE) == len(_CU_IM) == 84
assert len(_AU_XI) == len(_AU_RE) == len(_AU_IM) == 85
assert len(_AG_XI) == len(_AG_RE) == len(_AG_IM) == 43


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


MATERIAL_RANGE_NM = {1: (_CU_XI[0], _CU_XI[-1]), 2: (_AG_XI[0], _AG_XI[-1]), 4: (_AU_XI[0], _AU_XI[-1])}
MATERIAL_NAME = {1: "Cu", 2: "Ag", 4: "Au"}


def libmat(n_mat: int, lam_nm: float) -> complex:
    """Порт диспетчера libmat(n_mat, f, &Re_eps, &Im_eps) из COMPOSIT.c.

    Сопоставление индексов взято из комментария в COMPOSIT.c:
    "1- Cu, 2-Ag, 4-Au и т. д.". Сама libmat() не найдена в доступных
    исходниках - реализована здесь по этому комментарию и по единственной
    найденной материальной базе с таким же стилем (EPS_MET_AuAgCu.c).
    """
    if n_mat == 1:
        return eps_cu(lam_nm)
    if n_mat == 2:
        return eps_ag(lam_nm)
    if n_mat == 4:
        return eps_au(lam_nm)
    raise ValueError(
        f"libmat: материал n_mat={n_mat} не найден в доступных исходниках "
        "(известны только 1=Cu, 2=Ag, 4=Au; HfN/TiN из рис. 5 статьи не восстановлены)"
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
    print("[self_test] OK: composite_from_eps совпадает с аналитической формулой MG;")
    print("[self_test] OK: eps_au(530nm) и eps_cu(517nm) совпадают с табличными узлами.")


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

    lam_grid = np.arange(300.0, 2000.0 + 1e-9, 2.0)

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

    au_uv_mask = [lam for lam in au_lam_valid if lam <= 400.0]
    if au_uv_mask:
        idx = au_lam_valid.index(au_uv_mask[-1])
        print(f"[run] Im(eps_eff) для Au на границе УФ/видимого (~{au_lam_valid[idx]:.0f} nm) = {au_im[idx]:.3f}")
    if cu_lo > 400.0:
        print(f"[run] ВНИМАНИЕ: таблица Cu начинается только с {cu_lo:.0f} nm - "
              "сравнение Au/Cu в УФ (см. рис. 2 статьи, <400 nm) этим файлом не проверяется, "
              "данных для Cu в УФ в найденных исходниках нет.")

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


if __name__ == "__main__":
    self_test()
    run_reproduction()
