"""
Python-порт функций эффективной среды из vibr_t/src_c/VIBR.c:

- COMPOSITE      - формула Клаузиуса-Моссотти / Максвелла-Гарнетта, 1 тип
                   наночастиц: eps_c = eps_m*(1+2CK)/(1-CK),
                   K=(eps_p-eps_m)/(eps_p+2 eps_m);
- COMPOSITE_3    - обобщение на 2 типа наночастиц: S = C1*K1 + C2*K2,
                   eps_c = eps_m*(1+2S)/(1-S);
- COMPOSITE_me   - вариант "линейное смешение + поляризуемость":
                   eps_c = 1 + (1-C)(eps_m-1) + C*3 eps_m (eps_p-eps_m)/(eps_p+2 eps_m).

НАЙДЕННЫЙ БАГ (вероятная цель "Программы попытка отследить баг"):
в COMPOSITE_3 мнимая часть ЧИСЛИТЕЛЯ (1+2S) собрана из действительной:

    u1 = 1 + 2 * Re_sum;
    v1 = 2 * Re_sum;      // <-- должно быть 2 * Im_sum
    rh1 = 1.0 - Re_sum;
    ih1 = -1.0 * Im_sum;

Из-за этого даже при C2=0 COMPOSITE_3 НЕ совпадает с корректной COMPOSITE
(у которой числитель 1+2CK собран верно). Здесь есть обе версии:
composite3_as_written() (bug-for-bug, для воспроизведения расчётов Яковлева)
и composite3_fixed() (исправленная; при C2=0 совпадает с COMPOSITE точно).

Все eps - в конвенции Лерера e^{+i w t}: Im eps <= 0 для пассивной среды.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import vibr_materials as vm  # noqa: E402


# ---------------------------------------------------------------------------
# Чистые формулы (eps подаются напрямую, без libmat)
# ---------------------------------------------------------------------------

def mg_single(eps_m: complex, eps_p: complex, C: float) -> complex:
    """COMPOSITE: Клаузиус-Моссотти для одного типа включений."""
    K = (eps_p - eps_m) / (eps_p + 2.0 * eps_m)
    return eps_m * (1.0 + 2.0 * C * K) / (1.0 - C * K)


def composite3_fixed(eps_m1: complex, eps_m2: complex,
                     eps_p1: complex, eps_p2: complex,
                     C1: float, C2: float) -> complex:
    """COMPOSITE_3 с исправленным числителем (v1 = 2*Im_sum).

    При C2=0 и eps_m1=eps_m2 тождественно равна mg_single(eps_m1, eps_p1, C1).
    """
    K1 = (eps_p1 - eps_m1) / (eps_p1 + 2.0 * eps_m1)
    K2 = (eps_p2 - eps_m2) / (eps_p2 + 2.0 * eps_m2)
    S = C1 * K1 + C2 * K2
    return eps_m1 * (1.0 + 2.0 * S) / (1.0 - S)


def composite3_as_written(eps_m1: complex, eps_m2: complex,
                          eps_p1: complex, eps_p2: complex,
                          C1: float, C2: float) -> complex:
    """COMPOSITE_3 в точности как в C (баг сохранён: Im числителя = 2*Re_sum).

    Нужна для воспроизведения того, что реально считала программа
    (расчёт Яковлева с двумя типами наночастиц).
    """
    K1 = (eps_p1 - eps_m1) / (eps_p1 + 2.0 * eps_m1)
    K2 = (eps_p2 - eps_m2) / (eps_p2 + 2.0 * eps_m2)
    re_sum = (C1 * K1 + C2 * K2).real
    im_sum = (C1 * K1 + C2 * K2).imag
    numerator = complex(1.0 + 2.0 * re_sum, 2.0 * re_sum)  # БАГ: 2*re_sum вместо 2*im_sum
    denominator = complex(1.0 - re_sum, -im_sum)
    return eps_m1 * (numerator / denominator)


def composite_me(eps_m: complex, eps_p: complex, C: float) -> complex:
    """COMPOSITE_me: линейное смешение с поляризуемостью включения."""
    dip = 3.0 * eps_m * (eps_p - eps_m) / (eps_p + 2.0 * eps_m)
    re = 1.0 + (1.0 - C) * (eps_m.real - 1.0) + C * dip.real
    im = (1.0 - C) * eps_m.imag + C * dip.imag
    return complex(re, im)


# ---------------------------------------------------------------------------
# Обёртки, повторяющие сигнатуры C-функций (материал по номеру libmat)
# ---------------------------------------------------------------------------

def COMPOSITE(eps_m: complex, C: float, n_mat: int, lam_nm: float) -> complex:
    """void COMPOSITE(R_eps_m, I_eps_m, C, n_mat, ...) при f=lam_nm."""
    return mg_single(eps_m, vm.libmat(n_mat, lam_nm), C)


def COMPOSITE_3(eps_m1: complex, eps_m2: complex, C1: float, C2: float,
                n_mat_1: int, n_mat_2: int, lam_nm: float,
                fixed: bool = False) -> complex:
    """void COMPOSITE_3(...) при f=lam_nm. fixed=True - исправленный числитель."""
    eps_p1 = vm.libmat(n_mat_1, lam_nm)
    eps_p2 = vm.libmat(n_mat_2, lam_nm)
    fn = composite3_fixed if fixed else composite3_as_written
    return fn(eps_m1, eps_m2, eps_p1, eps_p2, C1, C2)


def COMPOSITE_me(eps_m: complex, C: float, n_mat: int, lam_nm: float) -> complex:
    return composite_me(eps_m, vm.libmat(n_mat, lam_nm), C)
