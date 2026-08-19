"""Прогон всех четырёх сцен подряд с итоговой сводкой проверок.

Каждая сцена сама себя проверяет и возвращает ненулевой код, если хотя бы одна
проверка не прошла. Здесь они запускаются последовательно, а в конце выводится
общий итог: сколько сцен построено и какие проверки провалились.

Запуск:
    python plasmon_fields_3d/scenes/run_all.py
"""

from __future__ import annotations

import io
import re
import sys
import time
from contextlib import redirect_stdout
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT.parent / "lrspp_coupling"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

SCENES = [
    ("1. Обычный плазмон", "scene1_single_spp"),
    ("2. Связанный плазмон", "scene2_coupled_spp"),
    ("3. Переход диэлектрик -> плазмон", "scene3_dielectric_to_plasmonic"),
    ("4. Фазовращатель на PCM", "scene4_pcm_phase_shifter"),
]


def main() -> int:
    import importlib

    total_ok = total_fail = 0
    report: list[str] = []

    for title, module_name in SCENES:
        print(f"=== {title}")
        started = time.time()
        buf = io.StringIO()
        module = importlib.import_module(module_name)
        with redirect_stdout(buf):
            code = module.main()
        text = buf.getvalue()

        ok = len(re.findall(r"\[OK\]", text))
        fail = len(re.findall(r"\[СБОЙ\]", text))
        total_ok += ok
        total_fail += fail
        elapsed = time.time() - started
        print(f"    проверок пройдено {ok}, провалено {fail}, время {elapsed:.0f} с")
        report.append(f"{title}: {ok} пройдено, {fail} провалено, {elapsed:.0f} с")
        if fail:
            for line in text.splitlines():
                if "[СБОЙ]" in line:
                    print(f"    {line.strip()}")
        if code != 0 and fail == 0:
            print(f"    сцена вернула код {code}")

    print()
    print("Итог")
    for line in report:
        print(f"  {line}")
    print(f"  всего: {total_ok} пройдено, {total_fail} провалено")
    print(f"  результаты в {ROOT / 'results'}")
    return 1 if total_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
