@echo off
REM Bit-Pop Windows build script for conda
REM Note: Bioconda focuses on Linux. Windows support via conda-forge.

REM Extract binary (Windows build)
tar -xf bit-pop-x86_64-windows.tar.gz

REM Install to conda bin directory
if not exist %PREFIX%\bin mkdir %PREFIX%\bin
copy /Y bit-pop.exe %PREFIX%\bin\bit-pop.exe
