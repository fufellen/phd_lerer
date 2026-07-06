# Примеры: три модели эффективной среды

Как вызывать модели эффективной диэлектрической проницаемости композита из
статьи «Влияние модели эффективной среды и источника оптических констант…» —
в C (для проекта Лерера) и в Python. Все примеры воспроизводят **таблицу 2**
статьи (положение и амплитуда пика `|Im eps_eff|` для Au и Cu).

## Быстрый запуск

```powershell
# Windows PowerShell — собирает C-пример (нужен gcc/MinGW) и запускает Python
.\composite_ema\examples\run_examples.ps1
```

```bat
REM Windows cmd
composite_ema\examples\run_examples.bat
```

```bash
# только Python (numpy; matplotlib не нужен)
python composite_ema/examples/ema_models_demo.py
```

Ожидаемый пик `|Im eps_eff|` (n_host=1.77, C=10%, константы Джонсона-Кристи):

| Модель | Au: пик, нм | Au: Im | Cu: пик, нм | Cu: Im |
|---|---:|---:|---:|---:|
| Клаузиус-Моссотти (1) | 572 | −6.05 | 596 | −3.32 |
| Бруггеман (4) | 708 | −2.66 | 702 | −2.56 |
| MLWA, R_p=10 нм | 580 | −6.69 | 600 | −3.80 |
| MLWA, R_p=30 нм | 640 | −6.49 | 640 | −5.60 |
| MLWA, R_p=50 нм | 790 | −3.59 | 786 | −3.37 |

## Использование в C (в проекте Лерера)

Функции живут в `../source_c/COMPOSIT.c` (Клаузиус-Моссотти) и
`../source_c/COMPOSITE_MODELS.c` (Бруггеман + MLWA). Они используют глобальную
длину волны `f` (нм), комплексные помощники `divide_`/`mult_` и выбор материала
`libmat(n_mat, f, &Re, &Im)` — то же окружение, что у существующей `COMPOSITE()`.

```c
double eps_m_re = 1.77 * 1.77, eps_m_im = 0.0;  /* матрица, eps = eps' - i eps'' */
double C = 0.10;                                /* объёмная доля частиц */
double re, im;
int n_mat = 4;                                  /* 4 = Au (см. комментарий в COMPOSIT.c) */

f = 535.0;                                      /* длина волны, нм (глобальная) */

COMPOSITE     (eps_m_re, eps_m_im, C, n_mat,        &re, &im);  /* формула (1) */
BRUGGEMAN     (eps_m_re, eps_m_im, C, n_mat,        &re, &im);  /* формула (4) */
COMPOSITE_MLWA(eps_m_re, eps_m_im, C, n_mat, 30.0,  &re, &im);  /* MLWA, R_p=30 нм */
```

Смесь двух сортов частиц (формула (5), исправленная реализация):

```c
COMPOSITE_3(eps_m_re, eps_m_im, C1, C2, n_mat_1, n_mat_2, &re, &im);
```

Если нужна модель без обращения к базе (eps частиц задаётся напрямую), есть
`*_eps`-ядра: `BRUGGEMAN_eps(...)`, `COMPOSITE_MLWA_eps(...)` — см. сигнатуры в
`COMPOSITE_MODELS.c`.

### Собрать автономный C-пример

`ema_models_demo.c` включает оба реальных `.c` и прогоняет их со встроенной
базой Au/Cu (Джонсон-Кристи). `VIBRr.h` в этой папке — минимальная заглушка
только для автономной сборки (в проекте Лерера используется его заголовок).

```powershell
cd composite_ema\examples
gcc -O2 -I. -o ema_models_demo.exe ema_models_demo.c -lm
.\ema_models_demo.exe
```

## Использование в Python

Модели: `composite_ema.composite_from_eps(eps_m, eps_p, C)` (Клаузиус-Моссотти),
`compare_effective_models.bruggeman_from_eps(eps_m, eps_p, C)` (Бруггеман),
`compare_effective_models.mlwa_mg_from_eps(eps_m, eps_p, C, lam_nm, R_p_nm)` (MLWA).

```python
import composite_ema as ce
import compare_effective_models as cm

eps_m = complex(1.77**2, 0.0)
eps_p = ce.eps_au(535.0)   # Au по Джонсону-Кристи, eps = eps' - i eps''
C = 0.10

print(ce.composite_from_eps(eps_m, eps_p, C))              # формула (1)
print(cm.bruggeman_from_eps(eps_m, eps_p, C))              # формула (4)
print(cm.mlwa_mg_from_eps(eps_m, eps_p, C, 535.0, 30.0))   # MLWA, R_p=30 нм
```

Полное диагностическое сравнение с графиком и CSV — `../scripts/compare_effective_models.py`.

## Соответствие C и Python

`ema_models_demo.c` и `ema_models_demo.py` печатают одинаковые числа. C-модели
`BRUGGEMAN`/`COMPOSITE_MLWA` сверены с Python-портом до 6 знаков (точечный
пример на 535 нм и все пики таблицы 2 совпадают).
