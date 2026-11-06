# -*- coding: utf-8 -*-
"""Запуск программы А. М. Лерера для прямоугольных столбиков (exe из письма 06.11.2026) из Python: одна длина волны - один процесс.

Исходников этой версии нет, только VIBR_T.__e и Vvod.c. Exe лежит вне git, путь - переменная VIBR_RECT_EXE или
Obsidian-служебное/.../от Лерера/Письма 06.11.2026/статья1/VIBR_T.__e на Диске.

Вопрос «Will you build the schedules» exe читает клавишей с консоли (getch), а не из стандартного ввода, поэтому без консоли
он висит без расхода процессора. Здесь процесс получает скрытую консоль (CREATE_NO_WINDOW), а клавиша пишется в неё.

Проницаемость столбика и слоя подаётся показателем преломления n + ik (n_eps = 1, k > 0), модель среды выбирает вызывающий код.
"""
import cmath
import ctypes
import os
import re
import shutil
import subprocess
import tempfile
import threading
import time
from ctypes import wintypes as W
from dataclasses import dataclass
from pathlib import Path

BUILD = Path(__file__).resolve().parent / "build"
_LOCK = threading.Lock()
_K = ctypes.WinDLL("kernel32", use_last_error=True)
_K.CreateFileW.restype = W.HANDLE
_K.CreateFileW.argtypes = [W.LPCWSTR, W.DWORD, W.DWORD, W.LPVOID, W.DWORD, W.DWORD, W.HANDLE]
_K.CloseHandle.argtypes = [W.HANDLE]


class _Key(ctypes.Structure):
    _fields_ = [("down", W.BOOL), ("rep", W.WORD), ("vk", W.WORD), ("sc", W.WORD), ("ch", W.WCHAR), ("ctrl", W.DWORD)]


class _Rec(ctypes.Structure):
    _fields_ = [("type", W.WORD), ("key", _Key)]


def exe_path():
    """Копия exe с обычным расширением в build/, источник - см. описание модуля."""
    dst = BUILD / "vibr_rect.exe"
    if dst.exists():
        return dst
    env = os.environ.get("VIBR_RECT_EXE")
    src = Path(env) if env else next(iter(Path.home().glob(
        "Мой диск*/Obsidian-служебное/PhD/Scientific_study/Scientific_study/paper_analysis/from_Lerer/от Лерера/Письма 06.11.2026/статья1/VIBR_T.__e")))
    BUILD.mkdir(exist_ok=True)
    shutil.copy2(src, dst)
    return dst


def _press(pid, ch="n"):
    """Нажатие клавиши в консоли процесса pid. Вызов закрыт замком: AttachConsole действует на весь процесс Python."""
    with _LOCK:
        _K.FreeConsole()
        if not _K.AttachConsole(pid):
            return False
        h = _K.CreateFileW("CONIN$", 0xC0000000, 3, None, 3, 0, None)
        recs = (_Rec * 2)(_Rec(1, _Key(True, 1, ord(ch.upper()), 0, ch, 0)), _Rec(1, _Key(False, 1, ord(ch.upper()), 0, ch, 0)))
        n = W.DWORD()
        ok = _K.WriteConsoleInputW(h, recs, 2, ctypes.byref(n))
        _K.CloseHandle(h)
        _K.FreeConsole()
        return bool(ok)


@dataclass
class Grating:
    period: float = 400.0
    ax: float = 250.0
    ay: float = 250.0
    h_pillar: float = 100.0
    h_layer: float = 2000.0
    h_buffer: float = 325.0
    n_buffer: float = 1.45
    n_sub: float = 1.45


@dataclass
class Basis:
    n_g: int = 20
    m_z: int = 5
    m_xy: int = 3


def nk(eps):
    """n + ik из eps с Im eps >= 0."""
    s = cmath.sqrt(eps)
    if s.imag < 0:
        s = -s
    return s.real, s.imag


def input_text(g, b, lam, eps_pillar, eps_layer):
    """Текст входного файла. Vvod.c умножает начальную длину волны на 1.00001, здесь это снято."""
    np_, kp = nk(eps_pillar)
    nl, kl = nk(eps_layer)
    x = lam / 1.00001
    return "\n".join([
        "/ GR / 0",
        "/ n_eps / 1",
        "/ if (GR) T /",
        "/ fl AF / 0 %.7f 100000 %.7f" % (x, lam + 0.5),
        "/ dx dy / %.6f %.6f" % (g.period, g.period),
        "/ N_fig / 0",
        "/ Nsl N_hol Nreg / 4 1 1",
        "/ is_l A_t B_t A_b B_b Xc Yc n k / 0 %.6f %.6f %.6f %.6f 0.0 0.0 %.8f %.8f" % (g.ax, g.ay, g.ax, g.ay, np_, kp),
        "/ n_[1] k_[1] / 1.0 0.0",
        "/ pillar layer / 0 %.6f 1.0 0.0" % g.h_pillar,
        "/ composite layer / 0 %.6f %.8f %.8f" % (g.h_layer, nl, kl),
        "/ buffer / 0 %.6f %.6f 0.0" % (g.h_buffer, g.n_buffer),
        "/ substrate / 0 %.6f 0.0" % g.n_sub,
        "/ ATET FI / 0.0 100 0.0 0.0",
        "/ AM_z / %d 100000 %d" % (b.m_z, b.m_z),
        "/ AM_y / %d 100000 %d" % (b.m_xy, b.m_xy),
        "/ ravn AM_x / 1 %d 100000 %d" % (b.m_xy, b.m_xy),
        "/ AN_g / %d 100000 %d" % (b.n_g, b.n_g),
        "/ Near / 0",
        "/ if Near /", "/ if Near /", "/ if Near /", "/ if Near /",
        "/ dB / 0",
        "", ""])


def read(path):
    """Строки нулевого порядка: {поляризация: (R, T, delta P)} и найденные числа базиса; строка '1.#QNAN' в разбор не попадёт."""
    text = open(path, "rb").read().decode("cp1251", "replace")
    head = re.search(r"N_g=\s*(\d+)\s+M_z=\s*(\d+)\s+M_x=\s*(\d+)\s+M_y=\s*(\d+)", text)
    calc = re.search(r"M_z_calc=\s*(\d+)\s+MS_=\s*(\d+)", text)
    rows, pol = {}, None
    for line in text.splitlines():
        m = re.match(r"\s*([sp]) - POLARIZATION", line)
        if m:
            pol = m[1]
            continue
        m = re.match(r"\s*0\s+0\s+([\d.]+)\s+([\d.]+)\s+([\d.eE+-]+)\s*$", line)
        if m and pol and pol not in rows:
            rows[pol] = (float(m[1]), float(m[2]), float(m[3]))
    return {"rows": rows, "basis": tuple(int(x) for x in head.groups()) if head else None,
            "M_z_calc": int(calc[1]) if calc else None, "MS": int(calc[2]) if calc else None, "text": text}


def run_text(text, timeout=1800, threads=1):
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), KMP_ALL_THREADS=str(threads))
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        (d / "in.txt").write_text(text, encoding="ascii")
        with open(d / "console.log", "wb") as con:
            p = subprocess.Popen([str(exe_path())], cwd=d, stdin=subprocess.PIPE, stdout=con, stderr=subprocess.STDOUT, env=env,
                                 creationflags=subprocess.CREATE_NO_WINDOW)
            p.stdin.write(b"1\nin.txt\nout.txt\n")
            p.stdin.close()
            t0 = time.time()
            # answer() сначала вычищает буфер клавиш, поэтому раннее нажатие может потеряться: повторять, лишние не мешают
            delay = 0.2
            while p.poll() is None:
                if time.time() - t0 > timeout:
                    p.kill()
                    p.wait()
                    raise TimeoutError("VIBR_T прямоугольные: нет ответа за %g с" % timeout)
                if time.time() - t0 < 60:
                    _press(p.pid)
                try:
                    p.wait(delay)
                except subprocess.TimeoutExpired:
                    pass
                delay = min(delay * 1.5, 2.0)
        r = read(d / "out.txt")
        r["seconds"] = time.time() - t0
        r["rc"] = p.returncode
        return r


def solve(g, lam, eps_pillar, eps_layer, basis=None, pol="s", timeout=1800, threads=1):
    """R, T, A на длине волны lam; A = delta P из вывода программы. В словаре - ещё время счёта и M_z_calc."""
    r = run_text(input_text(g, basis or Basis(), lam, eps_pillar, eps_layer), timeout, threads)
    if pol not in r["rows"]:
        raise RuntimeError("в выводе программы нет чисел (NaN?), lam = %g, код %s" % (lam, r["rc"]))
    return dict(zip("RTA", r["rows"][pol]), seconds=r["seconds"], M_z_calc=r["M_z_calc"], MS=r["MS"])
