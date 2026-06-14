param(
    [string]$Distro = $env:GEMELLI_WSL_DISTRO
)

if ([string]::IsNullOrWhiteSpace($Distro)) {
    $Distro = "Ubuntu"
}

$ErrorActionPreference = "Stop"

Write-Host "CorrOmics Joint-RPCA WSL setup"
Write-Host "WSL distribution: $Distro"

function Test-Command($Command) {
    $null = Get-Command $Command -ErrorAction SilentlyContinue
    return $?
}

if (-not (Test-Command "wsl.exe")) {
    Write-Host "WSL is not available. Requesting WSL installation..."
    Write-Host "Windows may ask for administrator approval and may require a reboot."
    Start-Process powershell -Verb RunAs -ArgumentList "wsl --install -d $Distro"
    Write-Host "After WSL finishes installing and Ubuntu is initialized, run this setup again."
    exit 1
}

$status = & wsl.exe --status 2>$null
if ($LASTEXITCODE -ne 0) {
    Write-Host "WSL is installed but not initialized. Requesting WSL installation/update..."
    Start-Process powershell -Verb RunAs -ArgumentList "wsl --install -d $Distro"
    Write-Host "After WSL finishes installing and Ubuntu is initialized, run this setup again."
    exit 1
}

& wsl.exe -d $Distro -- bash -lc "printf 'WSL ready: '; uname -a"
if ($LASTEXITCODE -ne 0) {
    Write-Host "The $Distro distribution is not ready."
    Write-Host "Installing it now. Windows may ask for approval or require first-launch setup."
    wsl.exe --install -d $Distro
    Write-Host "Open $Distro once, create the Linux username/password, then run this setup again."
    exit 1
}

$ScriptPath = Join-Path $PSScriptRoot "setup_gemelli_wsl.sh"
$WslScriptPath = (& wsl.exe -d $Distro -- wslpath -a "$ScriptPath").Trim()
if ($LASTEXITCODE -ne 0 -or [string]::IsNullOrWhiteSpace($WslScriptPath)) {
    throw "Could not translate setup script path for WSL: $ScriptPath"
}

Write-Host "Installing/updating Gemelli worker inside WSL..."
& wsl.exe -d $Distro -- bash "$WslScriptPath"
if ($LASTEXITCODE -ne 0) {
    throw "Gemelli WSL setup failed."
}

Write-Host ""
Write-Host "Joint-RPCA WSL setup complete."
Write-Host "Worker Python:"
& wsl.exe -d $Distro -- bash -lc 'printf "%s/.corromics/miniforge/envs/gemelli-standalone/bin/python\n" "$HOME"'
