@echo off

setlocal enabledelayedexpansion
set "file=E:\Desktop\QC_analysis\QC_runs_to_process_TIMS.txt"
set "count=0"
for /f "tokens=* delims=" %%a in ('type "%file%"') do (
    set /a count+=1
    if !count! equ 2 set "line=%%a"
)
for /f "tokens=2 delims=:" %%b in ("!line!") do set "title=%%b"
for /f "tokens=* delims= " %%c in ("!title!") do set "title=%%c"
title !title!
endlocal

"E:\Desktop\QC_analysis\src\R-4.3.2\bin\Rscript.exe" "E:\Desktop\QC_analysis\src\QC scripts\Version 3.3.6\Scripts\QC_analysis_3.3.6.R" "E:\Desktop\QC_analysis\log.txt" "E:\Desktop\QC_analysis\QC_runs_to_process_TIMS.txt" "E:\Desktop\QC_analysis\src\QC scripts\Version 3.3.6\Functions\\" "E:\Desktop\QC_analysis\src\Python-3.11\python.exe"

if exist "Rplots.pdf" del Rplots.pdf

pause