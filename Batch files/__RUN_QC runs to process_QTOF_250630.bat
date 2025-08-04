@echo off

setlocal enabledelayedexpansion
set "file=C:\Users\Admin\Desktop\Calibrations and QC analyses\QC analyses\QC_runs_to_process_QTOF_250630.txt"
set "count=0"
for /f "tokens=* delims=" %%a in ('type "%file%"') do (
    set /a count+=1
    if !count! equ 2 set "line=%%a"
)
for /f "tokens=2 delims=:" %%b in ("!line!") do set "title=%%b"
for /f "tokens=* delims= " %%c in ("!title!") do set "title=%%c"
title !title!
endlocal

"C:\Users\Admin\Desktop\Calibrations and QC analyses\src\R-4.3.2\bin\Rscript.exe" "C:\Users\Admin\Desktop\Calibrations and QC analyses\src\QC scripts\Version 3.3.11\Scripts\QC_analysis_3.3.11.R" "C:\Users\Admin\Desktop\Calibrations and QC analyses\QC analyses\\" "C:\Users\Admin\Desktop\Calibrations and QC analyses\QC analyses\QC_runs_to_process_QTOF_250630.txt" "C:\Users\Admin\Desktop\Calibrations and QC analyses\src\QC scripts\Version 3.3.11\Functions\\" "" "C:\Program Files (x86)\Mozilla Firefox\firefox.exe"

if exist "Rplots.pdf" del Rplots.pdf

pause