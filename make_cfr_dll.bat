@echo off
:: ---------------------------------------------------------------
:: Compile dynamic cfr library

:: run this script in the Cfs base directory
:: ---------------------------------------------------------------

:: library names
set lib1=cfr_arith
set lib2=cfr_math

call fbc_path.bat
if not exist %Fbas% exit

:: library directory
set ldir=.\library
:: path to the include files
set tool=.\modules

set bins=%tool%\bin

echo.
echo  compiling %lib2% DLL

set opts=-dll -O 2 -i %tool% -p %bins% -w pedantic -x

%Fbas% %ldir%\%lib1%.bas %ldir%\%lib2%.bas %opts% %bins%\%lib2%.dll >compile.log

del %bins%\*.a >nul
