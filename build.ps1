# Build and run the CliffordCore test suites.
#
#   .\build.ps1                       build and run every suite
#   .\build.ps1 pga                   only suites whose name contains "pga"
#
# Every directory under tests\ except common\ and compile_fail\ is one suite --
# one per algebra -- linked into one binary, build\test_<dir>.exe. Each .cpp in
# it is its own translation unit (one per kind of test, plus main.cpp); two
# algebras never meet in one translation unit or one binary. Binaries go to
# build\, which is gitignored.
#
# This runs the suites only. The standalone-header and compile-fail checks are
# bash scripts; ./build.sh runs them first.
#
# The repository path contains non-ASCII characters, so this maps the repo root
# to a spare drive letter before invoking g++ (see build_examples.ps1 for why
# the mapping is probed and released in a finally block).

[CmdletBinding()]
param(
    [string]$Filter = ''
)

$ErrorActionPreference = 'Stop'

$root = Split-Path -Parent $MyInvocation.MyCommand.Path

if (-not (Get-Command g++ -ErrorAction SilentlyContinue)) {
    throw 'g++ was not found. Install a C++ compiler and make sure it is on PATH.'
}

# Pick a drive letter that is not already in use.
$driveLetter = $null
foreach ($candidate in 'Z:', 'Y:', 'X:', 'W:') {
    if (-not (Test-Path "$candidate\")) { $driveLetter = $candidate; break }
}
if (-not $driveLetter) {
    throw 'No free drive letter available for the subst workaround (tried Z:, Y:, X:, W:).'
}

$previousLocation = Get-Location
$total = 0
$failed = 0

subst $driveLetter $root | Out-Null
try {
    Set-Location $driveLetter

    # Windows PowerShell 5.1 turns a native command's stderr into an ErrorRecord,
    # so with 'Stop' any g++ warning would terminate the whole script. We want
    # to report a broken suite and carry on to the next one.
    $ErrorActionPreference = 'Continue'

    $buildDir = Join-Path $driveLetter 'build'
    New-Item -ItemType Directory -Force -Path $buildDir | Out-Null

    $suites = @(Get-ChildItem -Directory (Join-Path $driveLetter 'tests') |
        Where-Object { $_.Name -notin @('common', 'compile_fail') } |
        Where-Object { $_.Name -like "*$Filter*" } |
        Sort-Object Name)

    foreach ($suite in $suites) {
        $name = "test_$($suite.Name)"
        $exe = Join-Path $buildDir "$name.exe"
        $total++

        $sources = @(Get-ChildItem (Join-Path $suite.FullName '*.cpp') |
            Sort-Object Name | ForEach-Object { $_.FullName })

        Write-Host "--- building $name ($($sources.Count) sources)"
        if ($sources.Count -eq 0) {
            Write-Host "BUILD FAILED: $name has no sources"
            $failed++
            continue
        }

        # Never leave an old binary behind for a suite that no longer builds.
        Remove-Item -Force -ErrorAction SilentlyContinue $exe

        # One g++ for the whole suite: it compiles each source in turn, then links.
        & g++ -std=c++17 -Wall -Wextra -Wpedantic -I (Join-Path $driveLetter 'include') $sources -o $exe
        if ($LASTEXITCODE -ne 0) {
            Write-Host "BUILD FAILED: $name"
            $failed++
            continue
        }

        Write-Host "=== $name ==="
        & $exe
        if ($LASTEXITCODE -ne 0) {
            Write-Host "RUN FAILED: $name"
            $failed++
        }
    }
}
finally {
    # Leave the drive before unmapping it, or the shell is stranded.
    Set-Location $previousLocation
    subst $driveLetter /D | Out-Null
}

Write-Host ''
if ($total -eq 0) {
    if ($Filter) {
        Write-Host "No test suites matched '$Filter'."
    } else {
        Write-Host "No test suites found in $root\tests\."
    }
    exit 1
}

Write-Host "$total test suites, $failed failed."
if ($failed -ne 0) { exit 1 }
