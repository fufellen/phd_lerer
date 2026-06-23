# ЭДП/FEM-проверка Ag-полоски W=800 нм

Эта папка хранит воспроизводимые скрипты метода эффективной диэлектрической проницаемости / effective-index method (ЭДП/EIM) и постобработку COMSOL-экспортов для проверки Ag-полоски на подложке.

## Структура

- `../Ag_.c` - исходные оптические константы серебра.
- `scripts/repeat_lerer_metal_strip.py` - планарная TM-задача для тонкой Ag-пленки и ЭДП-подобная горизонтальная оценка для `W=800 nm`.
- `scripts/width_sweep_ag_strip.py` - ЭДП-подобный sweep по ширине.
- `scripts/compare_comsol_width_sweep_500nm.py` - выбор COMSOL-мод для sweep по ширине при `lambda=500 nm`.
- `scripts/compare_comsol_lambda_sweep_w800.py` - выбор спектральной ветви через полевые диагностические признаки.
- `scripts/CreateAgStrip*.java` - Java-скрипты генерации COMSOL-моделей.
- `comsol_exports/` - маленькие CSV-таблицы, экспортированные из COMSOL Mode Analysis.
- `results/` - итоговые CSV/PNG-сводки для сравнения.

Большие COMSOL-модели (`*.mph`), class-файлы, логи, recovery- и status-файлы намеренно не хранятся в Git.

## Python-окружение

Из этой папки:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

SciPy удобен, но не обязателен для первого ЭДП-скрипта: `repeat_lerer_metal_strip.py` содержит fallback-решатель Ньютона, если SciPy не установлен. `numpy` и `matplotlib` нужны обязательно.

## Повторить ЭДП-расчеты

Из корня репозитория:

```powershell
python .\edp\metal_strip_w800\scripts\repeat_lerer_metal_strip.py
python .\edp\metal_strip_w800\scripts\width_sweep_ag_strip.py
```

Результаты пишутся в `edp/metal_strip_w800/results/`.

Ключевые файлы:

- `analytic_metal_strip_w800.csv`
- `analytic_metal_strip_w800.png`
- `analytic_width_sweep_ag_strip.csv`
- `analytic_width_sweep_vs_lambda.png`
- `analytic_width_dependence_500_600_700nm.png`

## Постобработка COMSOL CSV

Постобработчики предполагают, что CSV-экспорты COMSOL уже лежат в `comsol_exports/`.

```powershell
python .\edp\metal_strip_w800\scripts\compare_comsol_width_sweep_500nm.py
python .\edp\metal_strip_w800\scripts\compare_comsol_lambda_sweep_w800.py
```

Результаты пишутся в `results/`.

Спектральный постобработчик не выбирает фиксированный номер моды. Сначала он отбрасывает слабопотерянные light-line-подобные моды, затем использует диагностические полевые доли из `intopAg(ewfd.normE^2)`, `intopSub(...)` и интегралов по воздуху, а уже внутри такого пула выбирает ближайший комплексный `n_eff` к ЭДП-оценке со слабой проверкой непрерывности. Эти доли - полевые прокси локализации, а не нормированная мощность моды.

## COMSOL Java automation

Java-скрипты компилируются через COMSOL, например:

```powershell
& 'C:\Program Files\COMSOL\COMSOL62\Multiphysics\bin\win64\comsolcompile.exe' .\edp\metal_strip_w800\scripts\CreateAgStripW800LambdaSweep.java
```

Запускать compiled class нужно через `comsolbatch.exe`. После запуска проверяй batch log и наличие ожидаемых CSV, а не только код возврата процесса.

Java-скрипты хранятся как исходная автоматизация. Они могут создавать `.mph`, `.class`, `.log`, `.status` и `.recovery`; эти файлы игнорируются Git. Маленькие CSV-экспорты, которые нужны Python-постобработке, можно держать в `comsol_exports/`.

## Текущий диагностический результат

Для спектрального sweep `W=800 nm` коротковолновая область `450...550 nm` хорошо согласуется с ЭДП-подобной оценкой: отличие по `n'` меньше `1%`. В области `700...800 nm` ЭДП все заметнее завышает `n'` и занижает относительные потери `n''/n'`. Полевой прокси выбранной ветви меняется плавно: доля Ag падает, а доля подложки растет по мере приближения моды к световой линии подложки.
