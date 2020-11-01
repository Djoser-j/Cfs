@echo off
:: ---------------------------------------------------------------
:: Compile all cf demo programs

:: To be run in the Cfs base directory
:: ---------------------------------------------------------------

:: library name
set name=cfr_math

:: Set path to the freeBasic compiler:
set Fbas=".."
:: without closing backslash.

:: library directory
set ldir=.\library
:: source files directory
set tool=.\modules

if not exist %ldir%\lib%name%.dll.a call 1_make_cfr_dll.bat
echo.
echo  compiling all Cf demos...

set Cmdlin=%Fbas%\fbc -p %ldir% -w pedantic "%%~fa"

for %%a in (%tool%\*.bas) do %Cmdlin% >>compile.log

set Movlin=move "%%~fa" %tool%\bin

for %%a in (%tool%\*.exe) do %Movlin% >nul

echo.
pause
