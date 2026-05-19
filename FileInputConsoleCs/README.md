# FileInputConsoleCs

Console file-input application built with the C# compiler included in .NET Framework.
It does not require Qt, CMake, Visual Studio, NuGet packages, or downloaded runtime files.

Build with the included script:

```powershell
.\build.cmd
```

Run the sample:

```powershell
.\run_sample.cmd
```

Build:

```powershell
& C:\Windows\Microsoft.NET\Framework64\v4.0.30319\csc.exe /nologo /target:exe /out:FileInputConsole.exe Program.cs
```

Run:

```powershell
.\FileInputConsole.exe .\sample_input.txt
```
