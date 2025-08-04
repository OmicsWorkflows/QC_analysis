@echo off

setlocal enabledelayedexpansion
set "file=D:\Desktop\QC_analysis\QC_runs_to_process_HT_QC_HeLa.txt"
set "count=0"
for /f "tokens=* delims=" %%a in ('type "%file%"') do (
    set /a count+=1
    if !count! equ 2 set "line=%%a"
)
for /f "tokens=2 delims=:" %%b in ("!line!") do set "title=%%b"
for /f "tokens=* delims= " %%c in ("!title!") do set "title=%%c"
title !title!
endlocal

"D:\Desktop\QC_analysis\src\R-4.3.2\bin\Rscript.exe" "D:\Desktop\QC_analysis\src\QC scripts\Version 3.3.11\Scripts\QC_analysis_3.3.11.R" "D:\Desktop\QC_analysis\\" "D:\Desktop\QC_analysis\QC_runs_to_process_HT_QC_HeLa.txt" "D:\Desktop\QC_analysis\src\QC scripts\Version 3.3.11\Functions\\" "D:\Desktop\QC_analysis\src\Python-3.11\python.exe" "C:\Program Files\Mozilla Firefox\firefox.exe"

if exist "Rplots.pdf" del Rplots.pdf

pause