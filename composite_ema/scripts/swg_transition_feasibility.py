"""Пригоден ли композит поглотителя Лерера для волноводного перехода.

Вопрос. Поглотитель из прямоугольных столбиков согласует свободное пространство
с поглощающей стопкой на 400-750 нм. Можно ли ту же структуру - субволновую
решётку из композитных столбиков - применить как переход из диэлектрического
волновода в плазмонный?

Здесь считаются два числа, от которых ответ зависит количественно.

1. Погонные потери самого композита Au (10 %) в матрице n = 1.77 как функция
   длины волны. У поглотителя они и есть рабочий механизм, у перехода - потеря.
   Резонанс Фрёлиха наночастиц золота лежит внутри рабочей полосы поглотителя,
   поэтому там композит поглощает предельно сильно; вопрос в том, что с ним
   происходит на 1550 нм.

2. Граница субволнового режима. Решётка работает как эффективная среда, пока
   период меньше брэгговского: Lambda < lambda / (2 n_eff). Выше этой границы
   структура не согласует, а отражает. Период поглотителя 400 нм проверяется
   на этот критерий для 1550 нм.

Запуск:
    python composite_ema/scripts/swg_transition_feasibility.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from composite_ema import composite_from_eps, eps_au  # noqa: E402

OUT = HERE.parent / "results"
OUT.mkdir(parents=True, exist_ok=True)

N_MATRIX = 1.77          # матрица композита у Лерера
C_AU = 0.10              # объёмная доля золота, 10 %
PERIOD_NM = 400.0        # период решётки поглотителя
ABSORBER_BAND = (400.0, 750.0)
LAMBDA_TELECOM = 1550.0


def composite_index(lam_nm: float) -> complex:
    """Показатель композита в соглашении n = n' + i n'', n'' > 0 у поглощающей среды.

    В composite_ema принято соглашение с Im(eps) < 0, поэтому знак мнимой части
    приводится здесь явно, а не подразумевается.
    """
    eps_h = complex(N_MATRIX**2, 0.0)
    eps = composite_from_eps(eps_h, eps_au(lam_nm), C_AU)
    n = np.sqrt(complex(eps.real, -abs(eps.imag)))
    return complex(abs(n.real), abs(n.imag))


def loss_db_per_um(lam_nm: float, n: complex) -> float:
    k0 = 2.0 * np.pi / (lam_nm * 1e-3)
    return float(8.686 * k0 * n.imag)


def main() -> int:
    lam = np.arange(400.0, 1601.0, 5.0)
    n_re, n_im, alpha = [], [], []
    for lm in lam:
        n = composite_index(float(lm))
        n_re.append(n.real)
        n_im.append(n.imag)
        alpha.append(loss_db_per_um(float(lm), n))
    n_re = np.array(n_re)
    alpha = np.array(alpha)

    i_peak = int(np.argmax(alpha))
    n_tel = composite_index(LAMBDA_TELECOM)
    a_tel = loss_db_per_um(LAMBDA_TELECOM, n_tel)
    a_peak = float(alpha[i_peak])

    # граница субволнового режима: период должен быть меньше брэгговского
    period_max = lam / (2.0 * n_re)

    fig, axes = plt.subplots(1, 3, figsize=(16.0, 4.6))

    ax = axes[0]
    ax.semilogy(lam, alpha, color="tab:red")
    ax.axvspan(*ABSORBER_BAND, color="tab:orange", alpha=0.18,
               label="рабочая полоса поглотителя")
    ax.axvline(LAMBDA_TELECOM, color="tab:blue", ls="--", lw=1.2, label="1550 нм")
    ax.plot([lam[i_peak]], [a_peak], "o", color="k", ms=5)
    ax.annotate(f"{a_peak:.0f} дБ/мкм\nрезонанс Фрёлиха\n{lam[i_peak]:.0f} нм",
                (lam[i_peak], a_peak), textcoords="offset points", xytext=(24, -34), fontsize=8)
    ax.plot([LAMBDA_TELECOM], [a_tel], "o", color="tab:blue", ms=5)
    ax.annotate(f"{a_tel:.3f} дБ/мкм", (LAMBDA_TELECOM, a_tel),
                textcoords="offset points", xytext=(-70, 14), fontsize=8, color="tab:blue")
    ax.set_xlabel("длина волны, нм")
    ax.set_ylabel("погонные потери композита, дБ/мкм")
    ax.set_title(f"Композит Au {C_AU * 100:.0f} % в матрице n = {N_MATRIX}\n"
                 "то, что делает поглотитель поглотителем")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    ax = axes[1]
    ax.plot(lam, n_re, color="tab:green")
    ax.axvspan(*ABSORBER_BAND, color="tab:orange", alpha=0.18)
    ax.axvline(LAMBDA_TELECOM, color="tab:blue", ls="--", lw=1.2)
    ax.axhline(N_MATRIX, color="0.5", ls=":", lw=1, label=f"матрица без золота, n = {N_MATRIX}")
    ax.set_xlabel("длина волны, нм")
    ax.set_ylabel(r"$\mathrm{Re}\,n$ композита")
    ax.set_title(f"Показатель композита\nна 1550 нм n = {n_tel.real:.3f}")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[2]
    ax.plot(lam, period_max, color="tab:purple", label=r"граница $\Lambda=\lambda/(2n)$")
    ax.fill_between(lam, 0, period_max, color="tab:purple", alpha=0.12,
                    label="субволновый режим: работает как среда")
    ax.axhline(PERIOD_NM, color="k", ls="--", lw=1.4,
               label=f"период поглотителя {PERIOD_NM:.0f} нм")
    ax.axvline(LAMBDA_TELECOM, color="tab:blue", ls="--", lw=1.2)
    p_tel = LAMBDA_TELECOM / (2.0 * n_tel.real)
    ax.plot([LAMBDA_TELECOM], [p_tel], "o", color="tab:blue", ms=5)
    ax.annotate(f"предел {p_tel:.0f} нм", (LAMBDA_TELECOM, p_tel),
                textcoords="offset points", xytext=(-105, -26), fontsize=8, color="tab:blue")
    ax.set_xlabel("длина волны, нм")
    ax.set_ylabel("допустимый период, нм")
    ax.set_ylim(0, 600)
    ax.set_title("Где решётка ещё среда, а где уже зеркало")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.3)

    fig.suptitle(
        "Композит поглотителя как материал волноводного перехода: "
        "что мешает на 400-750 нм и что меняется на 1550 нм",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(OUT / "swg_transition_feasibility.png", dpi=150)
    plt.close(fig)

    lines = []
    add = lines.append
    add("Пригодность композита поглотителя для волноводного перехода")
    add(f"Композит: Au {C_AU * 100:.0f} % в матрице n = {N_MATRIX}, формула Клаузиуса-Моссотти")
    add("")
    add(f"{'lambda, нм':>11} {'Re n':>9} {'Im n':>11} {'потери, дБ/мкм':>16} {'предел периода, нм':>20}")
    for lm in (450.0, 520.0, 550.0, 633.0, 750.0, 1310.0, 1550.0):
        n = composite_index(lm)
        add(f"{lm:11.0f} {n.real:9.4f} {n.imag:11.5f} {loss_db_per_um(lm, n):16.4f} "
            f"{lm / (2.0 * n.real):20.0f}")
    add("")
    add(f"Максимум потерь {a_peak:.1f} дБ/мкм при {lam[i_peak]:.0f} нм - резонанс Фрёлиха")
    add(f"наночастиц золота в матрице n = {N_MATRIX}; он лежит внутри рабочей полосы")
    add(f"поглотителя {ABSORBER_BAND[0]:.0f}-{ABSORBER_BAND[1]:.0f} нм, и в этом весь его смысл.")
    add("")
    add(f"На {LAMBDA_TELECOM:.0f} нм тот же композит теряет {a_tel:.4f} дБ/мкм, то есть")
    add(f"в {a_peak / a_tel:.0f} раз меньше, при показателе {n_tel.real:.3f}.")
    add("")
    add(f"Предел периода на {LAMBDA_TELECOM:.0f} нм: {p_tel:.0f} нм.")
    add(f"Период поглотителя {PERIOD_NM:.0f} нм " +
        ("НЕ проходит" if PERIOD_NM > p_tel else "проходит") +
        " этот критерий, поэтому переносить его без пересчёта нельзя.")

    (OUT / "swg_transition_feasibility.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))

    failures = 0
    ok = a_peak / a_tel > 100
    print(f"\n  [{'OK' if ok else 'СБОЙ'}] потери на резонансе и на 1550 нм различаются "
          f"более чем на два порядка: в {a_peak / a_tel:.0f} раз")
    failures += 0 if ok else 1

    ok = ABSORBER_BAND[0] <= lam[i_peak] <= ABSORBER_BAND[1]
    print(f"  [{'OK' if ok else 'СБОЙ'}] максимум потерь лежит внутри рабочей полосы "
          f"поглотителя: {lam[i_peak]:.0f} нм")
    failures += 0 if ok else 1

    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
