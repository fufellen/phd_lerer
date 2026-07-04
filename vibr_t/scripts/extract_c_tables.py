"""
Автоматическая экстракция табличных оптических констант из C-исходников
материальной библиотеки VIBR_T (программа А.М. Лерера, получена 2026-07-04
через Льва Яковлева, папка Obsidian "Лев Яковлев/Программа попытка отследить
баг", скопирована в vibr_t/src_c/).

Зачем скрипт, а не ручная перепечатка: таблицы содержат тысячи чисел
(Au_.c/Ag_.c/Si_.c/ZnO_.c - шаг 1 нм, 200-2000(2001) нм), ручной перенос
не проверяем. Скрипт парсит присваивания массивов прямо из C и пишет CSV
в vibr_t/data/. Повторный запуск полностью воспроизводит data/ из src_c/.

Форматы, которые разбираются:
1. "Плотные" n,k-файлы (Au_.c, Ag_.c, Si_.c, ZnO_.c):
       Lam[m] = <нм>; n[m] = <n>; k[m] = <k>;
   индекс m совпадает с длиной волны в нм (200..2000, у Ag_/Si_/ZnO_ есть
   ещё точка 2001). Конвертация в eps в самих C-функциях: eps' = n^2-k^2,
   eps'' (хранимая) = -2nk (конвенция e^{+i w t}, Im eps <= 0).
   Здесь в CSV сохраняются исходные n,k без конвертации.
2. Разреженные eps-таблицы EPS_MET.c: Ag(xi/fs), Au(xi/fs), Cu(xii/fss);
   fs[m][0] = Re eps, fs[m][1] = Im eps (уже отрицательная).
   Файл разрезается по границам функций "double Ag(", "double Au(",
   "double Cu(".
3. AZO_.c: две независимые сетки (мкм) для Re eps (Lam/n) и для |Im eps|
   (Lam/k) - сохраняются двумя CSV.
4. VO2_cold.c / VO2_hot.c: xi (нм) + fs[m][0]=n', fs[m][1]=n'' (это n,k, не
   eps - конвертация в eps в C-функции).

SiO2_/HfO_ (SiO2_HfO.c) - аналитические формулы Селлмейера, таблиц нет,
портированы напрямую в vibr_materials.py.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

SRC_C = Path(__file__).resolve().parent.parent / "src_c"
DATA = Path(__file__).resolve().parent.parent / "data"

_ASSIGN_RE = re.compile(
    r"\b([A-Za-z_][A-Za-z0-9_]*)\s*\[\s*(\d+)\s*\]\s*(?:\[\s*(\d+)\s*\])?\s*=\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\s*;"
)


def parse_assignments(text: str) -> dict[str, dict]:
    """Собирает все присваивания вида name[i]=v; и name[i][j]=v; -> словари."""
    arrays: dict[str, dict] = {}
    for m in _ASSIGN_RE.finditer(text):
        name, i, j, val = m.group(1), int(m.group(2)), m.group(3), float(m.group(4))
        key = (i, int(j)) if j is not None else i
        arrays.setdefault(name, {})[key] = val
    return arrays


def strip_block_comments(text: str) -> str:
    """Убирает /* ... */ (данные в комментариях не должны попадать в CSV)."""
    return re.sub(r"/\*.*?\*/", "", text, flags=re.S)


def write_csv(path: Path, header: list[str], rows: list[list[float]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    print(f"[written] {path.name}: {len(rows)} точек")


def extract_dense_nk(fname: str, out_name: str) -> None:
    """Au_.c / Ag_.c / Si_.c / ZnO_.c -> lambda_nm,n,k"""
    text = strip_block_comments((SRC_C / fname).read_text(encoding="cp1251", errors="replace"))
    arr = parse_assignments(text)
    lam, n, k = arr["Lam"], arr["n"], arr["k"]
    idx = sorted(set(lam) & set(n) & set(k))
    if len(idx) != len(lam) or len(idx) != len(n) or len(idx) != len(k):
        raise ValueError(f"{fname}: несогласованные индексы Lam/n/k "
                         f"({len(lam)}/{len(n)}/{len(k)}, пересечение {len(idx)})")
    rows = [[lam[i], n[i], k[i]] for i in idx]
    write_csv(DATA / out_name, ["lambda_nm", "n", "k"], rows)


def split_functions(text: str, names: list[str]) -> dict[str, str]:
    """Разрезает C-файл на тела функций по 'double <name>('."""
    marks = []
    for nm in names:
        m = re.search(rf"double\s+{nm}\s*\(", text)
        if not m:
            raise ValueError(f"функция {nm} не найдена")
        marks.append((m.start(), nm))
    marks.sort()
    out = {}
    for (start, nm), nxt in zip(marks, marks[1:] + [(len(text), None)]):
        out[nm] = text[start:nxt[0]]
    return out


def extract_eps_met() -> None:
    """EPS_MET.c: Ag/Au (xi,fs), Cu (xii,fss) -> lambda_nm,eps_re,eps_im"""
    text = strip_block_comments((SRC_C / "EPS_MET.c").read_text(encoding="cp1251", errors="replace"))
    bodies = split_functions(text, ["Ag", "Au", "Cu"])
    for nm, (xi_name, fs_name) in {"Ag": ("xi", "fs"), "Au": ("xi", "fs"), "Cu": ("xii", "fss")}.items():
        arr = parse_assignments(bodies[nm])
        xi = arr[xi_name]
        fs = arr[fs_name]
        idx = sorted(i for i in xi if (i, 0) in fs and (i, 1) in fs)
        dropped = sorted(set(xi) - set(idx))
        if dropped:
            print(f"[note] EPS_MET.c {nm}: индексы {dropped} имеют xi, но не имеют fs - "
                  f"незаполненные точки исходного массива (в C это было бы чтение мусора)")
        rows = [[xi[i], fs[(i, 0)], fs[(i, 1)]] for i in idx]
        write_csv(DATA / f"eps_met_{nm}.csv", ["lambda_nm", "eps_re", "eps_im"], rows)


def extract_azo() -> None:
    """AZO_.c: сетки в мкм; отдельно Re eps (n) и |Im eps| (k)."""
    text = (SRC_C / "AZO_.c").read_text(encoding="cp1251", errors="replace")
    text = strip_block_comments(text)
    # Re-таблица идёт до 'if (lam < Lam[0] || lam > Lam[18])', k-таблица после
    parts = text.split("k[0]")
    arr_re = parse_assignments(parts[0])
    arr_im = parse_assignments("k[0]" + parts[1])
    rows_re = [[arr_re["Lam"][i], arr_re["n"][i]] for i in sorted(set(arr_re["Lam"]) & set(arr_re["n"]))]
    rows_im = [[arr_im["Lam"][i], arr_im["k"][i]] for i in sorted(set(arr_im["Lam"]) & set(arr_im["k"]))]
    write_csv(DATA / "azo_eps_re.csv", ["lambda_um", "eps_re"], rows_re)
    write_csv(DATA / "azo_eps_im_abs.csv", ["lambda_um", "eps_im_abs"], rows_im)


def extract_vo2(fname: str, out_name: str) -> None:
    text = strip_block_comments((SRC_C / fname).read_text(encoding="cp1251", errors="replace"))
    arr = parse_assignments(text)
    xi, fs = arr["xi"], arr["fs"]
    idx = sorted(i for i in xi if (i, 0) in fs and (i, 1) in fs)
    rows = [[xi[i], fs[(i, 0)], fs[(i, 1)]] for i in idx]
    write_csv(DATA / out_name, ["lambda_nm", "n", "k"], rows)


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    extract_dense_nk("Au_.c", "au_nk_dense.csv")
    extract_dense_nk("Ag_.c", "ag_nk_dense.csv")
    extract_dense_nk("Si_.c", "si_nk_dense.csv")
    extract_dense_nk("ZnO_.c", "zno_nk_dense.csv")
    extract_eps_met()
    extract_azo()
    extract_vo2("VO2_cold.c", "vo2_cold_nk.csv")
    extract_vo2("VO2_hot.c", "vo2_hot_nk.csv")
    print("\n[done] все таблицы извлечены в vibr_t/data/")


if __name__ == "__main__":
    main()
