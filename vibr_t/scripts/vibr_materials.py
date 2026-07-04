"""
Python-порт материальной библиотеки программы VIBR_T А.М. Лерера
(vibr_t/src_c/: Au_.c, Ag_.c, Si_.c, ZnO_.c, EPS_MET.c, AZO_.c, SiO2_HfO.c,
VO2_cold.c, VO2_hot.c, диспетчер BIBL_mat.c).

Табличные данные НЕ перепечатаны руками - читаются из vibr_t/data/*.csv,
которые генерирует scripts/extract_c_tables.py прямо из C-исходников.

Конвенция: у Лерера зависимость от времени e^{+i w t}, поэтому у пассивной
(поглощающей) среды Im eps <= 0; eps возвращается как complex(eps', eps'')
с eps'' <= 0. Для передачи в фотонные солверы с конвенцией e^{-i w t}
(наш tmm/rcwa) берите .conjugate().

Ключевые факты, установленные при порте (важно для статьи):
1. У МЕДИ НЕТ УФ-ДАННЫХ: таблица Cu в EPS_MET.c начинается с 517 нм.
   Алгоритм локальной интерполяции (mint=2) при lam < 517 нм линейно
   ЭКСТРАПОЛИРУЕТ по двум первым точкам (517 и 530 нм); уже при ~415 нм
   Re eps становится ПОЛОЖИТЕЛЬНОЙ - медь в программе Лерера в УФ ведёт
   себя как диэлектрик с потерями, а не как металл. Функции здесь по
   умолчанию воспроизводят это поведение (bug-for-bug), но помечают его:
   используйте strict=True, чтобы получить ошибку вне диапазона таблицы.
2. В "плотных" файлах (Au_.c и т.п.) интерполяционная формула
       n_ = n[m+1]*(Lam[m]-m_) + n[m]*(m_+1.0-Lam[m]);   m=(int)(lam+1e-7)
   при Lam[m]==m тождественно даёт n[floor(lam)] - кусочно-ПОСТОЯННУЮ
   интерполяцию (по-видимому, вместо Lam[m] имелось в виду lam). На целых
   длинах волн (как во всех входных файлах) разницы нет. По умолчанию
   используется линейная ("intended"); mode="c_exact" повторяет C.
3. В EPS_MET.c верхняя граница окна интерполяции клэмпится по ЗАЯВЛЕННОЙ
   ёмкости массива (q=44/86/85), а заполнено на 1 меньше (43/85/84) -
   за правым краем таблицы C-код читал неинициализированную память.
   Порт клэмпит по фактическому числу точек.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parent.parent / "data"


def _load_csv(name: str) -> np.ndarray:
    with (DATA / name).open("r", encoding="utf-8") as fh:
        rows = list(csv.reader(fh))
    return np.array([[float(x) for x in r] for r in rows[1:]], dtype=float)


class _Dense:
    """Плотная n,k-таблица (Au_.c/Ag_.c/Si_.c/ZnO_.c): шаг 1 нм, 200-2000 нм."""

    def __init__(self, csv_name: str, label: str):
        t = _load_csv(csv_name)
        self.lam = t[:, 0]
        self.n = t[:, 1]
        self.k = t[:, 2]
        self.label = label
        self.lo, self.hi = float(self.lam[0]), float(self.lam[-1])

    def eps(self, lam_nm: float, mode: str = "linear") -> complex:
        if lam_nm < self.lo or lam_nm > self.hi:
            # в C: printf("eps not defined"); exit(1);
            raise ValueError(f"{self.label}: lambda={lam_nm} нм вне таблицы "
                             f"[{self.lo}, {self.hi}] нм (в C программа завершалась)")
        if mode == "c_exact":
            m = int(lam_nm + 1.0e-7)
            i = int(np.searchsorted(self.lam, m))
            n_, k_ = self.n[i], self.k[i]
        else:
            n_ = float(np.interp(lam_nm, self.lam, self.n))
            k_ = float(np.interp(lam_nm, self.lam, self.k))
        return complex(n_ * n_ - k_ * k_, -2.0 * n_ * k_)


class _SparseEps:
    """Разреженная eps-таблица EPS_MET.c с локальной интерполяцией Лагранжа.

    Порт цикла:
        while(xi[i]<x)i++;
        b1=i-(mint/2); if(b1<1)b1=0; if(b1+mint-1>q)b1=q-mint+1;
        // далее полином Лагранжа по mint точкам окна
    mint=2 -> локально линейная интерполяция; за краями таблицы -
    линейная ЭКСТРАПОЛЯЦИЯ по двум крайним точкам (важно для Cu в УФ).
    """

    def __init__(self, csv_name: str, label: str, mint: int = 2):
        t = _load_csv(csv_name)
        self.lam = t[:, 0]
        self.re = t[:, 1]
        self.im = t[:, 2]
        self.label = label
        self.mint = mint
        self.lo, self.hi = float(self.lam[0]), float(self.lam[-1])

    def _lagrange(self, x: float, ys: np.ndarray) -> float:
        i = int(np.searchsorted(self.lam, x, side="left"))
        b1 = i - self.mint // 2
        if b1 < 1:
            b1 = 0
        if b1 + self.mint - 1 > len(self.lam) - 1:  # C клэмпил по заявленной ёмкости
            b1 = len(self.lam) - self.mint
        xs = self.lam[b1:b1 + self.mint]
        yw = ys[b1:b1 + self.mint]
        total = 0.0
        for ii in range(self.mint):
            li = yw[ii]
            for jj in range(self.mint):
                if jj != ii:
                    li *= (x - xs[jj]) / (xs[ii] - xs[jj])
            total += li
        return float(total)

    def eps(self, lam_nm: float, strict: bool = False) -> complex:
        if strict and (lam_nm < self.lo or lam_nm > self.hi):
            raise ValueError(f"{self.label}: lambda={lam_nm} нм вне таблицы "
                             f"[{self.lo}, {self.hi}] нм (strict)")
        return complex(self._lagrange(lam_nm, self.re), self._lagrange(lam_nm, self.im))

    def is_extrapolated(self, lam_nm: float) -> bool:
        return lam_nm < self.lo or lam_nm > self.hi


class _VO2:
    """VO2_cold.c/VO2_hot.c: n,k-таблица 1000-6000 нм, шаг 250 нм, линейно."""

    def __init__(self, csv_name: str, label: str):
        t = _load_csv(csv_name)
        self.lam = t[:, 0]
        self.n = t[:, 1]
        self.k = t[:, 2]
        self.label = label

    def eps(self, lam_nm: float) -> complex:
        if lam_nm < self.lam[0] or lam_nm > self.lam[-1]:
            # в C вне диапазона e1,e2 оставались неинициализированными (UB)
            raise ValueError(f"{self.label}: lambda={lam_nm} нм вне таблицы "
                             f"[{self.lam[0]}, {self.lam[-1]}] нм (в C - чтение мусора)")
        n_ = float(np.interp(lam_nm, self.lam, self.n))
        k_ = float(np.interp(lam_nm, self.lam, self.k))
        return complex(n_ * n_ - k_ * k_, -2.0 * n_ * k_)


class _AZO:
    """AZO_.c: eps'(лямбда) и |eps''|(лямбда) на РАЗНЫХ сетках, мкм."""

    def __init__(self):
        tr = _load_csv("azo_eps_re.csv")
        ti = _load_csv("azo_eps_im_abs.csv")
        self.lam_re, self.re = tr[:, 0], tr[:, 1]
        self.lam_im, self.im_abs = ti[:, 0], ti[:, 1]

    def eps(self, lam_nm: float) -> complex:
        lam_um = lam_nm / 1000.0
        if lam_um < self.lam_re[0] or lam_um > self.lam_re[-1]:
            raise ValueError(f"AZO_: lambda={lam_nm} нм вне таблицы Re eps")
        if lam_um < self.lam_im[0] or lam_um > self.lam_im[-1]:
            raise ValueError(f"AZO_: lambda={lam_nm} нм вне таблицы Im eps")
        re = float(np.interp(lam_um, self.lam_re, self.re))
        im = -float(np.interp(lam_um, self.lam_im, self.im_abs))
        return complex(re, im)


def _sio2(lam_nm: float) -> complex:
    """SiO2_ (Селлмейер, SiO2_HfO.c); без потерь."""
    L2 = (lam_nm / 1000.0) ** 2
    re = (1.0 + 0.6961663 * L2 / (L2 - 0.0684043 ** 2)
          + 0.4079426 * L2 / (L2 - 0.1162414 ** 2)
          + 0.8974794 * L2 / (L2 - 9.896161 ** 2))
    return complex(re, 0.0)


def _hfo(lam_nm: float) -> complex:
    """HfO_ (Селлмейер, SiO2_HfO.c); без потерь."""
    L2 = (lam_nm / 1000.0) ** 2
    re = (1.0 + 1.9558 * L2 / (L2 - 0.15494 ** 2)
          + 1.345 * L2 / (L2 - 0.0634 ** 2)
          + 10.41 * L2 / (L2 - 27.12 ** 2))
    return complex(re, 0.0)


# ---------------------------------------------------------------------------
# Ленивая инициализация таблиц
# ---------------------------------------------------------------------------

_cache: dict[str, object] = {}


def _get(name: str):
    if name not in _cache:
        if name == "au_dense":
            _cache[name] = _Dense("au_nk_dense.csv", "Au_ (плотная, luxpop-стиль)")
        elif name == "ag_dense":
            _cache[name] = _Dense("ag_nk_dense.csv", "Ag_ (плотная)")
        elif name == "si_dense":
            _cache[name] = _Dense("si_nk_dense.csv", "Si_ (плотная)")
        elif name == "zno_dense":
            _cache[name] = _Dense("zno_nk_dense.csv", "ZnO_ (плотная)")
        elif name == "cu_sparse":
            _cache[name] = _SparseEps("eps_met_Cu.csv", "Cu (EPS_MET.c, 517-3900 нм)")
        elif name == "ag_sparse":
            _cache[name] = _SparseEps("eps_met_Ag.csv", "Ag (EPS_MET.c, старая)")
        elif name == "au_sparse":
            _cache[name] = _SparseEps("eps_met_Au.csv", "Au (EPS_MET.c, старая)")
        elif name == "vo2_cold":
            _cache[name] = _VO2("vo2_cold_nk.csv", "VO2_cold")
        elif name == "vo2_hot":
            _cache[name] = _VO2("vo2_hot_nk.csv", "VO2_hot")
        elif name == "azo":
            _cache[name] = _AZO()
        else:
            raise KeyError(name)
    return _cache[name]


# ---------------------------------------------------------------------------
# Диспетчер libmat - порт закомментированного BIBL_mat.c.
# В присланном снимке ВЕСЬ BIBL_mat.c закомментирован (файл компилируется
# пустым), а libmat линкуется из готовой newlib+mat.lib - поэтому мапа
# восстановлена по тексту комментария и по доступным в снимке функциям.
# Случаи 9-19 (EPS_R, AZO_2, AZO_wt_*, GZO_, HfN_, ITO_, TaN_, TiN_, ZrN_)
# в снимке отсутствуют (нет исходников функций) - недоступны.
# Случаи 20/21 добавлены по тестовому коду в VIBR.c (сравнение VO2_cold/hot
# с libmat(20)/libmat(21)).
# ---------------------------------------------------------------------------

MATERIAL_NAMES = {
    1: "Cu (EPS_MET.c, 517-3900 нм, в УФ - экстраполяция!)",
    2: "Ag (Ag_.c, плотная 200-2000 нм)",
    3: "Au (Au_.c, плотная 200-2000 нм)",
    4: "ZnO (ZnO_.c, плотная 200-2000 нм)",
    5: "Si (Si_.c, плотная 200-2000 нм)",
    6: "SiO2 (Селлмейер)",
    7: "HfO2 (Селлмейер)",
    8: "AZO (таблица 1503-2695 нм)",
    20: "VO2 cold (1000-6000 нм)",
    21: "VO2 hot (1000-6000 нм)",
}


def libmat(is_: int, lam_nm: float, strict: bool = False) -> complex:
    """eps(lambda) материала номер is_ в конвенции Лерера (Im eps <= 0).

    strict=True: ошибка при выходе за диапазон таблицы (вместо
    молчаливой экстраполяции, как делает C-код для Cu в УФ).
    """
    if is_ == 1:
        return _get("cu_sparse").eps(lam_nm, strict=strict)
    if is_ == 2:
        return _get("ag_dense").eps(lam_nm)
    if is_ == 3:
        return _get("au_dense").eps(lam_nm)
    if is_ == 4:
        return _get("zno_dense").eps(lam_nm)
    if is_ == 5:
        return _get("si_dense").eps(lam_nm)
    if is_ == 6:
        return _sio2(lam_nm)
    if is_ == 7:
        return _hfo(lam_nm)
    if is_ == 8:
        return _get("azo").eps(lam_nm)
    if is_ == 20:
        return _get("vo2_cold").eps(lam_nm)
    if is_ == 21:
        return _get("vo2_hot").eps(lam_nm)
    raise KeyError(f"libmat: материал {is_} недоступен в этом снимке программы "
                   f"(есть: {sorted(MATERIAL_NAMES)})")


def cu_table_range() -> tuple[float, float]:
    """Диапазон реальных данных Cu (для демонстрации УФ-экстраполяции)."""
    t = _get("cu_sparse")
    return t.lo, t.hi


def au_sparse_eps(lam_nm: float, strict: bool = False) -> complex:
    """Старая разреженная таблица Au из EPS_MET.c (250-4000 нм) - для сверки."""
    return _get("au_sparse").eps(lam_nm, strict=strict)


def ag_sparse_eps(lam_nm: float, strict: bool = False) -> complex:
    """Старая разреженная таблица Ag из EPS_MET.c (400-2000 нм) - для сверки."""
    return _get("ag_sparse").eps(lam_nm, strict=strict)
