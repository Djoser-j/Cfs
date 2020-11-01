@echo off
:: ---------------------------------------------------------------
:: Run a single cf demo

:: To be run in the Cfs base directory
:: ---------------------------------------------------------------

:: demo to run
set file=Cf_Gosper

:: Set path to the source files directory:
set tool=.\modules\bin
:: Path to the working directory:
set work=.\workdir
:: without closing backslashes.

echo.
if not exist %tool%\%file%.exe goto :fail

echo  running %file%.exe

%tool%\%file% > %work%\%file%.log

echo.
exit

:fail
echo  %file%.exe not found
echo.
pause
