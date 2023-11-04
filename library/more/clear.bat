@echo off
:: Delete all compiles and logfiles
:: run this script in the Cfs base directory

:: library name
set name=cfr_math

:: executables directory
set tool=.\modules\bin
:: working directory
set work=.\workdir

del %tool%\%name%.dll >nul

del %tool%\*.exe >nul

del %work%\*.log >nul

del .\compare.txt >nul
del .\compile.log >nul
