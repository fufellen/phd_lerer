# Build and run the effective-medium model examples (Windows PowerShell).
# Needs gcc (MinGW-w64) in PATH for the C example, and python + numpy for Python.
$ErrorActionPreference = "Stop"
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $here
$OutputEncoding = [System.Text.Encoding]::UTF8
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8

Write-Host "=== C example (ema_models_demo.c) ===" -ForegroundColor Cyan
$gcc = Get-Command gcc -ErrorAction SilentlyContinue
if ($gcc) {
    & gcc -O2 -I. -o ema_models_demo.exe ema_models_demo.c -lm
    if ($LASTEXITCODE -eq 0) { & .\ema_models_demo.exe } else { Write-Host "build failed" -ForegroundColor Red }
} else {
    Write-Host "gcc not in PATH - C example skipped (Python example below prints the same numbers)." -ForegroundColor Yellow
}

Write-Host ""
Write-Host "=== Python example (ema_models_demo.py) ===" -ForegroundColor Cyan
python .\ema_models_demo.py
