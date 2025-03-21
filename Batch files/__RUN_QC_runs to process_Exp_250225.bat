@echo off

setlocal enabledelayedexpansion
set "file=D:\Desktop\QC analysis\QC_runs_to_process_Exp.txt"
set "count=0"
for /f "tokens=* delims=" %%a in ('type "%file%"') do (
    set /a count+=1
    if !count! equ 2 set "line=%%a"
)
for /f "tokens=2 delims=:" %%b in ("!line!") do set "title=%%b"
for /f "tokens=* delims= " %%c in ("!title!") do set "title=%%c"
title !title!
endlocal

"D:\Desktop\QC analysis\src\R-4.3.2\bin\Rscript.exe" "D:\Desktop\QC analysis\src\QC scripts\Version 3.3.6\Scripts\QC_analysis_3.3.6.R" "D:\Desktop\QC analysis\" "D:\Desktop\QC analysis\QC_runs_to_process_Exp.txt" "D:\Desktop\QC analysis\src\QC scripts\Version 3.3.6\Functions\\" ""

if exist "Rplots.pdf" del Rplots.pdf

pause
