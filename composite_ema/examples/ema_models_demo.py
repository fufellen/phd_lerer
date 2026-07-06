"""Пример использования трёх моделей эффективной среды (Python).

Мини-версия compare_effective_models.py: показывает, как из проницаемости
частиц eps_p и матрицы eps_m получить eps_eff композита каждой из трёх моделей,
и воспроизводит таблицу 2 статьи (пик |Im eps_eff| для Au и Cu).

Запуск:  python ema_models_demo.py
Требует numpy (matplotlib НЕ нужен). Соответствует C-примеру ema_models_demo.c.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

# читаемый вывод кириллицы в консоли Windows (cp866/cp1251 -> utf-8)
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

# модели и материальная база лежат в ../scripts
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import composite_ema as ce               # noqa: E402
import compare_effective_models as cm    # noqa: E402


def material_eps(name: str, lam_nm: float) -> complex:
    return ce.eps_au(lam_nm) if name == "Au" else ce.eps_cu(lam_nm)


def peak(name: str, model, lam_lo=300.0, lam_hi=900.0, step=2.0):
    """Возвращает (lambda пика, Im eps_eff в пике) для заданной модели."""
    lam = np.arange(lam_lo, lam_hi + 1e-9, step)
    im = np.array([model(material_eps(name, float(x)), float(x)).imag for x in lam])
    k = int(np.argmax(np.abs(im)))
    return lam[k], im[k]


def main() -> None:
    n_host = 1.77
    eps_m = complex(n_host ** 2, 0.0)   # 3.1329
    C = 0.10

    print(f"eps_eff композита (n_host={n_host}, C={C:.0%}), 3 модели. "
          "Материалы: Johnson & Christy 1972.\n")

    # --- точечный пример на 535 нм, Au ---
    ep = material_eps("Au", 535.0)
    print("Точечный пример на 535 нм, Au:")
    print(f"  Клаузиус-Моссотти -> {ce.composite_from_eps(eps_m, ep, C):+.4f}")
    print(f"  Бруггеман         -> {cm.bruggeman_from_eps(eps_m, ep, C):+.4f}")
    print(f"  MLWA (R_p=30 нм)  -> {cm.mlwa_mg_from_eps(eps_m, ep, C, 535.0, 30.0):+.4f}\n")

    # --- воспроизведение таблицы 2 ---
    models = [
        ("Клаузиус-Моссотти (1)", lambda ep, lam: ce.composite_from_eps(eps_m, ep, C)),
        ("Бруггеман (4)",         lambda ep, lam: cm.bruggeman_from_eps(eps_m, ep, C)),
        ("MLWA, R_p=10 нм",       lambda ep, lam: cm.mlwa_mg_from_eps(eps_m, ep, C, lam, 10.0)),
        ("MLWA, R_p=30 нм",       lambda ep, lam: cm.mlwa_mg_from_eps(eps_m, ep, C, lam, 30.0)),
        ("MLWA, R_p=50 нм",       lambda ep, lam: cm.mlwa_mg_from_eps(eps_m, ep, C, lam, 50.0)),
    ]
    for name in ("Au", "Cu"):
        print(f"{name}:")
        for label, model in models:
            lam_pk, im_pk = peak(name, model)
            print(f"  {label:<24} пик {lam_pk:.0f} нм,  Im eps_eff = {im_pk:+.3f}")
        print()

    print("Ожидание (таблица 2 статьи):")
    print("  Au: КМ 572/-6.05, Бруггеман 708/-2.66, MLWA 580/-6.69, 640/-6.49, 790/-3.59")
    print("  Cu: КМ 596/-3.32, Бруггеман 702/-2.56, MLWA 600/-3.80, 640/-5.60, 786/-3.37")


if __name__ == "__main__":
    main()
