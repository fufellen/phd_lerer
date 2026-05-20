@echo off
setlocal
chcp 65001 >nul

set "ROOT=%~dp0"
set "CONFIG=%~1"
if "%CONFIG%"=="" set "CONFIG=Debug"

set "CMAKE_EXE="
for /f "delims=" %%I in ('where cmake 2^>nul') do (
    if not defined CMAKE_EXE set "CMAKE_EXE=%%I"
)

if not defined CMAKE_EXE (
    set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
    if exist "%VSWHERE%" (
        for /f "usebackq delims=" %%I in (`"%VSWHERE%" -latest -products * -property installationPath`) do (
            if exist "%%I\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe" (
                set "CMAKE_EXE=%%I\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe"
            )
        )
    )
)

if not defined CMAKE_EXE (
    echo Error: CMake was not found in PATH.
    echo Install CMake, or install Visual Studio/Build Tools with the CMake tools component.
    exit /b 1
)

"%CMAKE_EXE%" -S "%ROOT%" -B "%ROOT%build" -DCMAKE_BUILD_TYPE=%CONFIG%
if errorlevel 1 exit /b %errorlevel%

"%CMAKE_EXE%" --build "%ROOT%build" --config %CONFIG%
if errorlevel 1 exit /b %errorlevel%

echo.
echo Built executable:
echo   %ROOT%build\bin\%CONFIG%\phd_lerer.exe
