# -*- coding: utf-8 -*-
"""Сборка VIBR_T обычным gcc, без Visual Studio 2013 и Intel MKL.

src_c/ и newLIb/ не меняются: они копируются в build/, правки PATCHES применяются к копии.
  python build_vibr.py            -> port_c/build/vibr_t.exe
  python build_vibr.py --trap     -> port_c/build_trap/vibr_t.exe: отладочная сборка, недопустимая операция останавливает программу (для gdb)
Нужен gcc в PATH. Ключи совпадают с Release исходного проекта по смыслу: оптимизация, без быстрой арифметики.
"""
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
VIBR = HERE.parent
LIB = VIBR.parent / "newLIb"
TRAP = "--trap" in sys.argv
WORK = HERE / ("build_trap" if TRAP else "build")

# файлы из VIBR_T.vcxproj
PROJECT = ["Ancillary.c", "BASIS.c", "BIBL_mat.c", "CONST.c", "CONVERG.c", "DIAGON_N.c", "D_N_PER.c", "DIAGON.c",
           "DROITE.c", "E_extern.c", "EXTERN.c", "F_df.c", "Fun.c", "MATR_.c", "MEMORY.c", "MLAP.C", "NOMER.c",
           "RECUR.c", "RECUR_ZERO.c", "VIBR.c", "Vvod.c", "ZERO_new.c"]
# таблицы материалов: в исходном проекте они приходили из библиотеки newlib+mat
MATERIALS = ["EPS_MET.c", "Ag_.c", "Au_.c", "Si_.c", "ZnO_.c", "SiO2_HfO.c", "AZO_.c", "VO2_cold.c", "VO2_hot.c"]
# из newLIb: только то, что нужно компоновщику
LIBRARY = ["MAT1D.C", "MAT2D.C", "Mat3d.c", "MAT4D.c", "FSKCOM.C", "TIM.C", "HORKUZ.C", "MEMCHK.C"]

# (файл, что заменить, на что). Только расхождения компиляторов, расчёт не затрагивается.
PATCHES = [
    ("MLAP.C", b"(char *)p_blas += acol*size", b"p_blas = (char *)p_blas + acol*size"),
    ("VIBRr.h", b"__inline void", b"static inline void"),
]
HEADERS = {
    "mkl.h": '#include "mkl_stub.h"\n', "mkl_blas.h": '#include "mkl_stub.h"\n', "mkl_cblas.h": '#include "mkl_stub.h"\n',
    "mkl_types.h": '#include "mkl_stub.h"\n', "mkl_vml.h": '#include "mkl_stub.h"\n',
    "newlib+mat.h": '#pragma once\n#include "newLib.h"\nvoid libmat(int is, double lam, double *Re_eps, double *Im_eps);\n',
}


def main():
    if WORK.exists():
        shutil.rmtree(WORK)
    (WORK / "lib").mkdir(parents=True)
    for p in (VIBR / "src_c").iterdir():
        if p.is_file():
            shutil.copy2(p, WORK / p.name)
    for p in LIB.iterdir():
        if p.is_file():
            shutil.copy2(p, WORK / "lib" / p.name)
    for p in (HERE / "compat").iterdir():
        shutil.copy2(p, WORK / p.name)
    for name, text in HEADERS.items():
        (WORK / name).write_text(text, encoding="ascii")
    for name, old, new in PATCHES:
        f = WORK / name
        data = f.read_bytes()
        if old not in data:
            sys.exit("правка не применилась: %s" % name)
        f.write_bytes(data.replace(old, new))

    flags = (["-g", "-O0", "-DFP_TRAP"] if TRAP else ["-O2"]) + ["-w", "-fcommon", "-I.", "-Ilib"]
    objs = []
    jobs = [(n, ["-std=gnu89"]) for n in PROJECT + MATERIALS]
    jobs += [("lib/" + n, ["-std=gnu89"]) for n in LIBRARY]
    jobs += [("vibr_compat.c", ["-std=gnu99"])]
    bad = 0
    for name, std in jobs:
        obj = name.replace("/", "_").rsplit(".", 1)[0] + ".o"
        r = subprocess.run(["gcc"] + std + flags + ["-c", "-x", "c", name, "-o", obj], cwd=WORK, capture_output=True,
                           text=True, encoding="utf-8", errors="replace")
        if r.returncode:
            bad += 1
            print("ОШИБКА", name, r.stderr[:500])
        objs.append(obj)
    if bad:
        sys.exit("не скомпилировано файлов: %d" % bad)
    r = subprocess.run(["gcc"] + objs + ["-o", "vibr_t.exe", "-lm", "-Wl,--stack,268435456"], cwd=WORK, capture_output=True,
                       text=True, encoding="utf-8", errors="replace")
    if r.returncode:
        print(r.stderr[-3000:])
        sys.exit("компоновка не прошла")
    print("собрано:", WORK / "vibr_t.exe")


main()
