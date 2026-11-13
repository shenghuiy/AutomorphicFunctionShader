@echo off
rem Builds the native renderers with the MSVC x64 toolchain (Visual Studio 2022 or Build Tools).
rem
rem   build.bat                          poincare.exe, icosahedral.exe and pearcey.exe in native\build
rem   build.bat poincare                 native\build\poincare.exe only
rem   build.bat icosahedral              native\build\icosahedral.exe only
rem   build.bat pearcey                  native\build\pearcey.exe only
rem   build.bat poincare other.exe      native\build\other.exe (e.g. while poincare.exe is running)
rem
rem All outputs (.exe, .lib, .exp, .obj) go to native\build, the parent folder's build folder.
setlocal
cd /d "%~dp0"
set "TARGET=%~1"
if "%TARGET%"=="" set "TARGET=all"
if /i not "%TARGET%"=="all" if /i not "%TARGET%"=="poincare" if /i not "%TARGET%"=="icosahedral" if /i not "%TARGET%"=="pearcey" if /i not "%TARGET%"=="kantorovich" (
  echo Unknown target "%TARGET%". Use poincare, icosahedral, pearcey, kantorovich, or no argument for all.
  exit /b 1
)
for /f "usebackq tokens=*" %%i in (`"%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set "VSPATH=%%i"
if not defined VSPATH (
  echo Visual Studio with the C++ x64 tools was not found.
  exit /b 1
)
call "%VSPATH%\VC\Auxiliary\Build\vcvars64.bat" >nul 2>&1 || exit /b 1
set "OUT=..\build"
if not exist "%OUT%" mkdir "%OUT%" || exit /b 1
if /i "%TARGET%"=="all" (
  call :build poincare poincare.exe || exit /b 1
  call :build icosahedral icosahedral.exe || exit /b 1
  call :build pearcey pearcey.exe || exit /b 1
  call :build kantorovich kantorovich.exe || exit /b 1
  exit /b 0
)
set "OUTPUT=%~2"
if "%OUTPUT%"=="" set "OUTPUT=%TARGET%.exe"
call :build %TARGET% "%OUTPUT%"
exit /b %errorlevel%

:build
rem %1 = program folder (also the source name), %2 = output executable name in %OUT%
cl /nologo /O2 /EHsc /std:c++17 /utf-8 /W3 /D_CRT_SECURE_NO_WARNINGS /Icommon %1\%1.cpp ^
   /Fo:"%OUT%\\" /Fe:"%OUT%\%~2" /link user32.lib d3d12.lib dxgi.lib d3dcompiler.lib
exit /b %errorlevel%
