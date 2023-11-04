@echo off
:: ---------------------------------------------------------------
:: Run all cf demo programs

:: run this script in the Cfs base directory
:: ---------------------------------------------------------------

:: executables directory
set tool=.\modules\bin
:: working directory
set work=.\workdir

echo.
if not exist %tool%\*.exe (
  echo  demo.exes not found
  echo.
  pause
  exit
)

echo  running all cf demos...

del %work%\*.log >nul 2>&1

for %%a in (%tool%\*.exe) do "%%~fa" > %work%\"%%~na".log
