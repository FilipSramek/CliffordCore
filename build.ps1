$ErrorActionPreference = 'Stop'

$root = Split-Path -Parent $MyInvocation.MyCommand.Path
$driveLetter = 'Z:'

$compiler = Get-Command g++ -ErrorAction SilentlyContinue
if (-not $compiler) {
    throw 'g++ was not found. Please install a C++ compiler and make sure it is available on PATH.'
}

try {
    subst $driveLetter /D | Out-Null
} catch {}

subst $driveLetter $root | Out-Null

$buildDir = Join-Path $driveLetter 'build'
$source = Join-Path $driveLetter 'tests/test_core.cpp'
$output = Join-Path $buildDir 'test_core.exe'

New-Item -ItemType Directory -Force -Path $buildDir | Out-Null
Set-Location $driveLetter

& g++ -std=c++17 -I (Join-Path $driveLetter 'include') $source -o $output
& $output

subst $driveLetter /D | Out-Null
