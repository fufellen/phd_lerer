@echo off
REM Build and run the effective-medium model examples (Windows cmd).
REM Needs gcc (MinGW-w64) in PATH for the C example, and python + numpy for Python.
chcp 65001 >nul
cd /d "%~dp0"

echo === C example (ema_models_demo.c) ===
where gcc >nul 2>nul
if %errorlevel%==0 (
    gcc -O2 -I. -o ema_models_demo.exe ema_models_demo.c -lm
    if %errorlevel%==0 ( ema_models_demo.exe ) else ( echo build failed )
) else (
    echo gcc not in PATH - C example skipped ^(Python example below prints the same numbers^).
)

echo.
echo === Python example (ema_models_demo.py) ===
python "%~dp0ema_models_demo.py"
