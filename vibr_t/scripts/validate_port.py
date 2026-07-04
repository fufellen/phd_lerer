"""
Валидация Python-порта материальной библиотеки и композитных формул VIBR_T.

Проверки (все должны проходить ДО использования порта в статье):
1. Узловые значения таблиц: eps в точности из n,k CSV-строк (конверсия
   eps'=n^2-k^2, eps''=-2nk) и из eps-строк EPS_MET.c.
2. Линейная интерполяция между узлами.
3. Экстраполяция Cu ниже 517 нм воспроизводит C-алгоритм (окно из двух
   первых точек) и демонстрирует нефизичную положительную Re eps в УФ.
4. Формульный паритет COMPOSITE с независимым портом composite_ema
   (composite_from_eps) на сетке аргументов.
5. composite3_fixed(C2=0) == mg_single (машинная точность).
6. composite3_as_written(C2=0) != mg_single - численная демонстрация бага.
7. Условие Фрёлиха Re eps_p = -2*eps_m для констант Лерера (Au, Cu).
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "composite_ema" / "scripts"))

import vibr_materials as vm  # noqa: E402
import vibr_composite as vc  # noqa: E402
import composite_ema as ema  # noqa: E402

FAIL = 0


def check(name: str, ok: bool, detail: str = "") -> None:
    global FAIL
    status = "OK " if ok else "FAIL"
    print(f"[{status}] {name}" + (f" - {detail}" if detail else ""))
    if not ok:
        FAIL += 1


def read_csv(name: str):
    with (HERE.parent / "data" / name).open("r", encoding="utf-8") as fh:
        rows = list(csv.reader(fh))
    return [[float(x) for x in r] for r in rows[1:]]


def main() -> None:
    # --- 1. Узловые значения: плотная таблица Au (материал 3) ------------
    au = read_csv("au_nk_dense.csv")
    row520 = next(r for r in au if r[0] == 520.0)
    n_, k_ = row520[1], row520[2]
    want = complex(n_ * n_ - k_ * k_, -2 * n_ * k_)
    got = vm.libmat(3, 520.0)
    check("Au_ узел 520 нм", abs(got - want) < 1e-12, f"eps={got:.6f}")

    # узловые значения: разреженная таблица Cu (материал 1)
    cu = read_csv("eps_met_Cu.csv")
    lam0, re0, im0 = cu[0]
    got = vm.libmat(1, lam0)
    check(f"Cu узел {lam0:.0f} нм (первая точка таблицы)",
          abs(got - complex(re0, im0)) < 1e-12, f"eps={got:.4f}")

    # --- 2. Линейная интерполяция между узлами ---------------------------
    lam_mid = 0.5 * (cu[0][0] + cu[1][0])
    want = complex(0.5 * (cu[0][1] + cu[1][1]), 0.5 * (cu[0][2] + cu[1][2]))
    got = vm.libmat(1, lam_mid)
    check(f"Cu середина отрезка {lam_mid:.1f} нм", abs(got - want) < 1e-12)

    # --- 3. Экстраполяция Cu в УФ (главный научный факт) ------------------
    lo, hi = vm.cu_table_range()
    check("Диапазон таблицы Cu начинается с 517 нм", lo == 517.0, f"lo={lo}, hi={hi}")
    # C-алгоритм: при x < xi[0] окно = точки 0,1 -> линейная экстраполяция
    x = 250.0
    slope_re = (cu[1][1] - cu[0][1]) / (cu[1][0] - cu[0][0])
    want_re = cu[0][1] + slope_re * (x - cu[0][0])
    got = vm.libmat(1, x)
    check("Cu 250 нм: экстраполяция по двум первым точкам (как в C)",
          abs(got.real - want_re) < 1e-9, f"Re eps={got.real:.3f}")
    check("Cu 250 нм: Re eps > 0 - в УФ медь у Лерера НЕ металл",
          got.real > 0, f"Re eps={got.real:.3f} (нефизично для Cu)")
    # где экстраполированная Re eps пересекает ноль
    lams = np.arange(250.0, 520.0, 0.5)
    re_vals = np.array([vm.libmat(1, float(l)).real for l in lams])
    crossings = lams[np.where(np.diff(np.sign(re_vals)) != 0)]
    if len(crossings):
        print(f"       Re eps_Cu(экстрап.) = 0 при ~{crossings[0]:.0f} нм; "
              f"ниже - 'диэлектрическая' медь")
    # строгий режим должен ловить выход за таблицу
    try:
        vm.libmat(1, 250.0, strict=True)
        check("Cu strict=True вне таблицы -> ошибка", False)
    except ValueError:
        check("Cu strict=True вне таблицы -> ошибка", True)

    # --- 4. Паритет COMPOSITE с composite_ema ----------------------------
    eps_h = complex(1.77 ** 2, 0.0)
    rng = np.random.default_rng(42)
    max_d = 0.0
    for _ in range(200):
        eps_p = complex(rng.uniform(-150, 5), rng.uniform(-30, 0))
        C = rng.uniform(0.01, 0.3)
        d = abs(vc.mg_single(eps_h, eps_p, C) - ema.composite_from_eps(eps_h, eps_p, C))
        max_d = max(max_d, d)
    check("COMPOSITE == composite_ema.composite_from_eps (200 случайных точек)",
          max_d < 1e-12, f"max|delta|={max_d:.2e}")

    # --- 5. composite3_fixed(C2=0) == mg_single ---------------------------
    max_d = 0.0
    for lam in np.arange(250.0, 900.0, 10.0):
        eps_p = vm.libmat(3, float(lam))
        d = abs(vc.composite3_fixed(eps_h, eps_h, eps_p, eps_p, 0.10, 0.0)
                - vc.mg_single(eps_h, eps_p, 0.10))
        max_d = max(max_d, d)
    check("composite3_fixed(C2=0) == COMPOSITE", max_d < 1e-12, f"max|delta|={max_d:.2e}")

    # --- 6. Демонстрация бага composite3_as_written ----------------------
    lam_demo = 554.0
    eps_p = vm.libmat(3, lam_demo)
    bugged = vc.composite3_as_written(eps_h, eps_h, eps_p, eps_p, 0.10, 0.0)
    correct = vc.mg_single(eps_h, eps_p, 0.10)
    rel = abs(bugged - correct) / abs(correct)
    check("composite3_as_written(C2=0) отличается от COMPOSITE (баг воспроизведён)",
          rel > 0.01,
          f"lam={lam_demo} нм: bug={bugged:.4f}, correct={correct:.4f}, delta={rel:.1%}")

    # --- 7. Условие Фрёлиха на константах Лерера --------------------------
    target = -2.0 * eps_h.real
    for name, mat, lam_grid in (("Au", 3, np.arange(400.0, 700.0, 0.25)),
                                ("Cu", 1, np.arange(517.0, 700.0, 0.25))):
        re_vals = np.array([vm.libmat(mat, float(l)).real for l in lam_grid])
        idx = np.where(np.diff(np.sign(re_vals - target)) != 0)[0]
        if len(idx):
            l0, l1 = lam_grid[idx[0]], lam_grid[idx[0] + 1]
            r0, r1 = re_vals[idx[0]], re_vals[idx[0] + 1]
            lam_f = l0 + (target - r0) * (l1 - l0) / (r1 - r0)
            print(f"       Фрёлих {name} (Re eps=-2*{eps_h.real:.4f}): {lam_f:.1f} нм "
                  f"(константы Лерера)")
            check(f"Фрёлих {name} найден в видимом диапазоне", 400 <= lam_f <= 700)
        else:
            check(f"Фрёлих {name} найден в видимом диапазоне", False, "нет пересечения")

    # --- 8. Ошибки диапазонов у остальных материалов ----------------------
    for mat, bad_lam in ((3, 150.0), (4, 2500.0), (20, 500.0)):
        try:
            vm.libmat(mat, bad_lam)
            check(f"материал {mat}: lambda={bad_lam} -> ошибка", False)
        except ValueError:
            check(f"материал {mat}: lambda={bad_lam} -> ошибка", True)

    print()
    if FAIL:
        print(f"ИТОГ: {FAIL} проверок ПРОВАЛЕНО")
        sys.exit(1)
    print("ИТОГ: все проверки пройдены - порт согласован с C-исходниками")


if __name__ == "__main__":
    main()
