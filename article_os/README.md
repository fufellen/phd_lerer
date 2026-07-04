# article_os - фигуры и числа для статьи в «Оптика и спектроскопия»

Генерация журнальных рисунков и числовых фактов для рукописи
«Влияние модели эффективной среды и источника оптических констант на
расчётные характеристики поглощающих метаповерхностей с плазмонными
наночастицами» (рукопись живёт в Obsidian-вольте:
`.../от Лерера/Статья ОиС/Метаповерхности с ПНЧ_рукопись.md`).

## Скрипты

```powershell
python article_os/scripts/make_figures.py      # fig1..fig7 (PNG 600 dpi + EPS) в article_os/figures/
python article_os/scripts/metrics_for_text.py  # числа, цитируемые в тексте рукописи
```

Оба скрипта только читают уже провалидированные CSV из `composite_ema/`,
`grating_2d_rcwa/` и `vibr_t/` - ничего не пересчитывают. Поэтому числа в
тексте статьи всегда согласованы с кривыми на рисунках.

## Требования журнала, реализованные в make_figures.py

- подписи осей и весь текст в поле рисунка - только английский;
- кривые пронумерованы (расшифровка в подрисуночных подписях рукописи);
- без сетки; стили линий различимы в градациях серого;
- каждый рисунок отдельным файлом, PNG (600 dpi) + EPS (вектор).

## Соответствие рисунков данным

| Рисунок | Источник данных |
|---|---|
| fig1 | схема (рисуется кодом, данных нет) |
| fig2 | `composite_ema/results/composite_eps_Au_vs_Cu.csv` |
| fig3 | `composite_ema/results/effective_model_comparison_AuCu.csv` |
| fig4 | `grating_2d_rcwa/results/reproduce_fig2_au_vs_cu.csv` + `..._2d.csv` |
| fig5 | `grating_2d_rcwa/results/reproduce_fig5_hfn_vs_tin.csv` |
| fig6 | `vibr_t/results/compare_constants_au_cu.csv` |
| fig7 | `vibr_t/results/rerun_fig2_lerer_constants.csv` |

Статус: `validated` - исходные CSV получены провалидированными скриптами
соответствующих подпроектов (см. их README), фигуры и метрики - чистое
отображение этих данных.
