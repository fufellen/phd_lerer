# QtFileInputConsole

Small Qt console application that reads an input file path from an argument or from stdin, opens the file, and prints its contents with line numbers.

Build with CMake:

```powershell
cmake -S . -B build
cmake --build build
.\build\Debug\QtFileInputConsole.exe .\sample_input.txt
```

Build with qmake:

```powershell
qmake QtFileInputConsole.pro
nmake
.\QtFileInputConsole.exe .\sample_input.txt
```
