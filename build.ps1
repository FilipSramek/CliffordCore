# Build and run the CliffordCore test suites.
#
#   .\build.ps1                       build and run every tests\*.cpp
#
# Every tests\*.cpp is its own translation unit and its own binary -- one per
# algebra, so the algebras never meet in a single TU. Binaries go to build\,
# which is gitignored.
#
# The repository path contains non-ASCII characters, so this maps the repo root
# to a spare drive letter before invoking g++ (see build_examples.ps1 for why
# the mapping is probed and released in a finally block).

[CmdletBinding()]
param()

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

    $sources = @(Get-ChildItem (Join-Path $driveLetter 'tests\*.cpp') | Sort-Object Name)

    foreach ($src in $sources) {
        $name = $src.BaseName
        $exe = Join-Path $buildDir "$name.exe"
        $total++

        Write-Host "--- building $name"
        & g++ -std=c++17 -I (Join-Path $driveLetter 'include') $src.FullName -o $exe
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
    Write-Host "No test sources found in $root\tests\."
    exit 1
}

Write-Host "$total test suites, $failed failed."
if ($failed -ne 0) { exit 1 }
