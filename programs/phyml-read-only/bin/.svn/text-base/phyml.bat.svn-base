@echo off
cls
color 1B
mode con:cols=98
mode con:lines=50
if not exist phyml.exe goto error
if exist phyml.exe then goto launch
echo on


:launch
phyml
goto end

:error
echo Error - can't find `phyml.exe'
pause
goto end


:end
echo Execution finished
pause
