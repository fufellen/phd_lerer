@echo off
setlocal

set "CSC=%SystemRoot%\Microsoft.NET\Framework64\v4.0.30319\csc.exe"
if not exist "%CSC%" set "CSC=%SystemRoot%\Microsoft.NET\Framework\v4.0.30319\csc.exe"

if not exist "%CSC%" (
    echo Error: C# compiler was not found.
    exit /b 1
)

"%CSC%" /nologo /target:exe /out:FileInputConsole.exe Program.cs
if errorlevel 1 exit /b %errorlevel%

echo Built: %CD%\FileInputConsole.exe
