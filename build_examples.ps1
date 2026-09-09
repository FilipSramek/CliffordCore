# Build and run the CliffordCore examples.
#
#   .\build_examples.ps1              build and run everything
#   .\build_examples.ps1 03           only examples whose name contains "03"
#   .\build_examples.ps1 -NoRun       build only, do not execute
#
# Binaries go to build\examples\, which is gitignored.
#
# The repository path contains non-ASCII characters, so this maps the repo root
# to a spare drive letter before invoking g++, matching build.ps1. Three things
# build.ps1 gets wrong and this does not:
#   1. a compile failure there leaves the drive still mapped   -> try/finally
#   2. it leaves the shell sitting on the unmapped drive       -> restore location
#   3. it runs `subst <letter> /D` before mapping, destroying a
#      mapping the user may already have                       -> probe for a free one

[CmdletBinding()]
param(
    [string]$Filter = '',
    [switch]$NoRun
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
    # so with 'Stop' any g++ warning or error would terminate the whole script.
    # We want to report a broken example and carry on to the next one, so relax
    # this for the build loop. Cleanup is handled by the finally block regardless.
    $ErrorActionPreference = 'Continue'

    $outDir = Join-Path $driveLetter 'build\examples'
    New-Item -ItemType Directory -Force -Path $outDir | Out-Null

    $sources = @(Get-ChildItem (Join-Path $driveLetter 'examples\*.cpp') | Sort-Object Name)
    if ($Filter) {
        $sources = @($sources | Where-Object { $_.BaseName -like "*$Filter*" })
    }

    foreach ($src in $sources) {
        $name = $src.BaseName
        $exe = Join-Path $outDir "$name.exe"
        $total++

        Write-Host "--- building $name"
        & g++ -std=c++17 -I (Join-Path $driveLetter 'include') $src.FullName -o $exe
        if ($LASTEXITCODE -ne 0) {
            Write-Host "BUILD FAILED: $name"
            $failed++
            continue
        }

        if (-not $NoRun) {
            Write-Host "=== $name ==="
            & $exe
            if ($LASTEXITCODE -ne 0) {
                Write-Host "RUN FAILED: $name"
                $failed++
            }
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
    if ($Filter) { Write-Host "No examples matched '$Filter'." }
    else { Write-Host "No examples found in $root\examples\." }
    exit 1
}

Write-Host "$total examples, $failed failed."
if ($failed -ne 0) { exit 1 }
