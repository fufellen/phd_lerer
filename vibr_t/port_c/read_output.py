# -*- coding: utf-8 -*-
"""Разбор вывода VIBR_T: блоки по числу базисных функций, для каждой длины волны R, T, delta P нулевого порядка по поляризациям.

  python read_output.py <вывод> [<эталон>]   -> сводка по блокам; с эталоном - наибольшие отличия
Сумма R и T по всем порядкам здесь не считается: при периоде меньше длины волны порядок один.
"""
import re
import sys


def read(path):
    text = open(path, "rb").read().decode("cp1251", "replace")
    blocks, cur, lam, pol = [], None, None, None
    for line in text.splitlines():
        m = re.search(r"N_g=\s*(\d+)\s+M_z=\s*(\d+)\s+M_fi=\s*(\d+)\s+M_r=\s*(\d+)", line)
        if m:
            cur = {"N_g": int(m[1]), "M_z": int(m[2]), "M_fi": int(m[3]), "M_r": int(m[4]), "rows": {}}
            blocks.append(cur)
            continue
        m = re.search(r"lamda\s*=\s*([\d.]+)", line)
        if m:
            lam = float(m[1])
            continue
        if re.search(r"\bs - POLARIZATION", line):
            pol = "s"
            continue
        if re.search(r"\bp - POLARIZATION", line):
            pol = "p"
            continue
        m = re.match(r"\s*0\s+0\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s*$", line)
        if m and cur is not None and lam is not None and pol:
            cur["rows"][(lam, pol)] = (float(m[1]), float(m[2]), float(m[3]))
    return blocks


def main():
    a = read(sys.argv[1])
    for i, b in enumerate(a):
        lams = sorted({k[0] for k in b["rows"]})
        print("блок %d: N_g=%d M_z=%d M_fi=%d M_r=%d, длин волн %d, %s" % (
            i, b["N_g"], b["M_z"], b["M_fi"], b["M_r"], len(lams), "%g..%g" % (lams[0], lams[-1]) if lams else "-"))
    if len(sys.argv) > 2:
        r = read(sys.argv[2])
        print("эталон:")
        for i, b in enumerate(r):
            lams = sorted({k[0] for k in b["rows"]})
            print("блок %d: N_g=%d M_z=%d M_fi=%d M_r=%d, длин волн %d, %s" % (
                i, b["N_g"], b["M_z"], b["M_fi"], b["M_r"], len(lams), "%g..%g" % (lams[0], lams[-1]) if lams else "-"))
        for i, x in enumerate(a):
            for j, y in enumerate(r):
                keys = set(x["rows"]) & set(y["rows"])
                if not keys:
                    continue
                d = [max(abs(p - q) for p, q in zip(x["rows"][k], y["rows"][k])) for k in keys]
                k = max(keys, key=lambda k: max(abs(p - q) for p, q in zip(x["rows"][k], y["rows"][k])))
                print("наш блок %d - эталон %d: общих строк %d, наибольшее отличие %.5f при %s: %s против %s" % (
                    i, j, len(keys), max(d), k, x["rows"][k], y["rows"][k]))


if __name__ == "__main__":
    main()
