@echo off
:: ---------------------------------------------------------------
:: Compile a single cf demo

:: To be run in the Cfs base directory
:: ---------------------------------------------------------------

:: module to compile
set file=Cf_Gosper

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
echo  compiling %file%.exe

set opts= -p %ldir% -w pedantic

%Fbas%\fbc %opts% %tool%\%file%.bas >compile.log

move %tool%\%file%.exe %tool%\bin >nul

echo.
rem pause
