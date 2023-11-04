@echo off
:: ---------------------------------------------------------------
:: Command line: run <demo name> without .bas

:: run this script in the Cfs base directory
:: ---------------------------------------------------------------

:: executables directory
set tool=.\modules\bin
:: working directory
set work=.\workdir

echo Run demo: 
set demo=%1

echo.
if not exist %tool%\%demo%.exe (
  echo not found: %tool%\%demo%.exe
  echo.
  pause
  exit
)

%tool%\%demo%.exe > %work%\%demo%.log

more /e < %work%\%demo%.log
