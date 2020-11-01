@echo off
:: ---------------------------------------------------------------
:: Delete all compiles

:: To be run in the Cfs base directory
:: ---------------------------------------------------------------

:: Library name
set name=cfr_math

:: working directory
set work=.\workdir
:: library directory
set ldir=.\library
:: source files directory
set tool=.\modules\bin

echo.

del %ldir%\lib%name%.dll.a
del %tool%\%name%.dll

for %%a in (%tool%\*.exe) do del "%%~fa" >nul

for %%a in (%work%\*.log) do del "%%~fa" >nul

exit
