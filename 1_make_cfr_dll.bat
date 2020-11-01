@echo off
:: ---------------------------------------------------------------
:: Compile dynamic cfr libraries
:: To be run in the Cfs base directory

:: Note: the Cfs archive you downloaded is assumed to be
:: unpacked to the base directory of your freeBasic installation.
:: ---------------------------------------------------------------

:: Library names
set lib1=cfr_arith
set lib2=cfr_math

:: Set path to the library directory:
set ldir=.\library
:: without closing backslash.

echo.
echo  compiling %lib2% DLL

pushd %ldir%

:: Path to the freeBasic compiler
set Fbas="..\.."
:: Source files directory
set tool=..\modules

set Cmdlin=%Fbas%\fbc -c -i %tool% -w pedantic "%%~fa"

for %%a in (*.bas) do %Cmdlin% >>compile.log

%Fbas%\fbc -dll -a %lib2%.o %lib1%.o >>compile.log

move %lib2%.dll %tool%\bin

for %%a in (*.o) do del "%%~fa" >nul

echo.
rem pause
