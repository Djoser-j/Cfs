@echo off
:: ---------------------------------------------------------------
:: Command line: make <demo name> without .bas

:: run this script in the Cfs base directory
:: ---------------------------------------------------------------

:: library name
set name=cfr_math

call fbc_path.bat
if not exist %Fbas% exit

:: source files directory
set tool=.\modules

set bins=%tool%\bin

if not exist %bins%\%name%.dll call make_cfr_dll.bat

echo Make demo: 
set demo=%1

echo.
if not exist %tool%\%demo%.bas (
  echo not found: %tool%\%demo%.bas
  echo.
  pause
  exit
)

echo  compiling demo %demo%.bas

set opts=-O 2 -p %bins% -s console -w pedantic -x

%Fbas% %tool%\%demo%.bas %opts% %bins%\%demo%.exe >compile.log
