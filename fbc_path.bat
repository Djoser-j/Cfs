@echo off
:: ---------------------------------------------------------------
:: Relative path to the FreeBASIC compiler

:: run this script in the base directory
:: ---------------------------------------------------------------

set Fbas="..\fbc.exe"

if not exist %Fbas% (
  echo.
  echo  FB compiler not in path
  echo.
  pause
)
