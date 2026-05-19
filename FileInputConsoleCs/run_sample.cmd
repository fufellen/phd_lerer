@echo off
setlocal

call "%~dp0build.cmd"
if errorlevel 1 exit /b %errorlevel%

"%~dp0FileInputConsole.exe" "%~dp0sample_input.txt"
