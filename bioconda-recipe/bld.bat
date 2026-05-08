@echo off
copy /Y bit-pop-x86_64-windows.zip bit-pop.zip
tar -xf bit-pop.zip bit-pop.exe
move /Y bit-pop.exe %PREFIX%\bin\bit-pop.exe
