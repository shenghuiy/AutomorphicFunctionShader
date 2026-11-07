@echo off
rem Builds poincare.exe with the MSVC x64 toolchain (Visual Studio 2022 or Build Tools).
setlocal
set "OUTPUT=poincare.exe"
if not "%~1"=="" set "OUTPUT=%~1"
cd /d "%~dp0"
for /f "usebackq tokens=*" %%i in (`"%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set "VSPATH=%%i"
if not defined VSPATH (
  echo Visual Studio with the C++ x64 tools was not found.
  exit /b 1
)
call "%VSPATH%\VC\Auxiliary\Build\vcvars64.bat" >nul 2>&1 || exit /b 1
cl /nologo /O2 /EHsc /std:c++17 /utf-8 /W3 /D_CRT_SECURE_NO_WARNINGS poincare.cpp /Fe:"%OUTPUT%" /link user32.lib d3d12.lib dxgi.lib d3dcompiler.lib
