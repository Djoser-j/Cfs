@echo off
:: ---------------------------------------------------------------
:: Compile all cf demo programs

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
echo.
echo  compiling all cf demos...

set cmnd=%Fbas% -O 2 -p %bins% -s console -w pedantic "%%~fa"

for %%a in (%tool%\*.bas) do %cmnd% >>compile.log

move %tool%\*.exe %bins% >nul

echo.
pause
