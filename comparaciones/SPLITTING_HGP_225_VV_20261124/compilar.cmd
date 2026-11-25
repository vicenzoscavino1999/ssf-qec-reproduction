@echo off
setlocal
for /f "usebackq tokens=*" %%i in (`"%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe" -latest -products * -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath`) do set "QEC_VS=%%i"
if not defined QEC_VS exit /b 1
call "%QEC_VS%\VC\Auxiliary\Build\vcvars64.bat" >nul
if errorlevel 1 exit /b 1
cl /nologo /O2 /EHsc /std:c++17 /arch:AVX2 /W4 /Fe:motor_225.exe motor_225.cpp
exit /b %errorlevel%
