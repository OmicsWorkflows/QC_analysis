@echo off

setlocal enabledelayedexpansion
set "file=D:\Desktop\QC analysis\QC_runs_to_process_Exp_250602.txt"
set "count=0"
for /f "tokens=* delims=" %%a in ('type "%file%"') do (
    set /a count+=1
    if !count! equ 2 set "line=%%a"
)
for /f "tokens=2 delims=:" %%b in ("!line!") do set "title=%%b"
for /f "tokens=* delims= " %%c in ("!title!") do set "title=%%c"
title !title!
endlocal

"D:\Desktop\QC analysis\src\R-4.3.2\bin\Rscript.exe" "D:\Desktop\QC analysis\src\QC scripts\Version 3.3.9\Scripts\QC_analysis_3.3.9.R" "D:\Desktop\QC analysis\\" "D:\Desktop\QC analysis\QC_runs_to_process_Exp_250602.txt" "D:\Desktop\QC analysis\src\QC scripts\Version 3.3.9\Functions\\" "" "C:\Program Files\Mozilla Firefox\firefox.exe"

if exist "Rplots.pdf" del Rplots.pdf

pause
