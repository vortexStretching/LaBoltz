param(
  [string]$Configuration = "Debug",
  [switch]$Clean
)

$ErrorActionPreference = "Stop"

$Root = Resolve-Path "$PSScriptRoot\.."
$VcVars = "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat"

if (-not (Test-Path $VcVars)) {
  throw "Could not find Visual Studio vcvars64.bat at '$VcVars'. Install the C++ Build Tools workload first."
}

if ($Clean -and (Test-Path "$Root\build")) {
  Remove-Item -LiteralPath "$Root\build" -Recurse -Force
}

$Command = @(
  "call `"$VcVars`"",
  "cmake -S `"$Root`" -B `"$Root\build`"",
  "cmake --build `"$Root\build`" --config $Configuration",
  "ctest --test-dir `"$Root\build`" -C $Configuration --output-on-failure"
) -join " && "

cmd.exe /c $Command

