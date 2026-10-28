$ErrorActionPreference = 'Stop'
$script:LauncherRoot = $PSScriptRoot
$script:Url = 'http://localhost:8503'

function Test-DockerReady {
    # Windows PowerShell can treat native stderr as an error even when redirected.
    $ErrorActionPreference = 'Continue'
    & docker info *> $null
    return ($LASTEXITCODE -eq 0)
}

function Invoke-DockerCompose {
    param([string[]]$Arguments)

    $ErrorActionPreference = 'Continue'
    & docker compose @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw ('Docker Compose failed (exit code ' + $LASTEXITCODE + '). Check the Docker error above.')
    }
}

function Wait-AppReady {
    param([int]$TimeoutSeconds = 60)

    $deadline = (Get-Date).AddSeconds($TimeoutSeconds)
    do {
        $response = $null
        try {
            $response = Invoke-WebRequest -Uri 'http://127.0.0.1:8503/api/health' -UseBasicParsing -TimeoutSec 3
        } catch {
            # The container may still be starting.
        }
        if ($null -ne $response -and $response.StatusCode -eq 200) {
            $health = $null
            try {
                $health = $response.Content | ConvertFrom-Json
            } catch {
                # Do not open a browser for an unrelated HTTP server.
            }
            if ($null -eq $health -or $health.app -cne 'ambs-mapper') {
                throw 'Port 8503 responded with an unexpected app. Check which program is using this port.'
            }
            return
        }
        Start-Sleep -Seconds 2
    } while ((Get-Date) -lt $deadline)

    throw 'amBS Mapper did not respond within 60 seconds. Run: docker compose logs --tail 100 ambs-mapper'
}

function Invoke-Launcher {
    Set-Location -LiteralPath $script:LauncherRoot
    New-Item -ItemType Directory -Force -Path (Join-Path $script:LauncherRoot 'data\input') | Out-Null

    if (-not (Get-Command docker -ErrorAction SilentlyContinue)) {
        throw 'Docker was not found. Install Docker Desktop, then run start.bat again.'
    }

    if (-not (Test-DockerReady)) {
        $desktopPath = Join-Path $env:ProgramFiles 'Docker\Docker\Docker Desktop.exe'
        if (-not (Test-Path -LiteralPath $desktopPath)) {
            throw 'Docker is not running. Start Docker Desktop in Linux containers mode, then try again.'
        }
        Write-Host 'Starting Docker Desktop...'
        Start-Process -FilePath $desktopPath -WindowStyle Hidden
        $dockerDeadline = (Get-Date).AddSeconds(120)
        do {
            Start-Sleep -Seconds 3
            if (Test-DockerReady) { break }
        } while ((Get-Date) -lt $dockerDeadline)
        if (-not (Test-DockerReady)) {
            throw 'Docker did not become ready within 120 seconds. Check Docker Desktop, then try again.'
        }
    }

    $runningIds = @(Invoke-DockerCompose -Arguments @('ps', '--status', 'running', '--quiet', 'ambs-mapper') |
        Where-Object { -not [string]::IsNullOrWhiteSpace($_) })
    if ($runningIds.Count -gt 0) {
        Write-Host 'amBS Mapper is already running. Keeping the current container and analysis.'
    } else {
        Write-Host 'Starting amBS Mapper (the first build may take several minutes)...'
        Invoke-DockerCompose -Arguments @('up', '-d', '--build', 'ambs-mapper')
    }

    Write-Host 'Waiting for amBS Mapper...'
    Wait-AppReady
    Start-Process $script:Url
    Write-Host ('amBS Mapper is ready: ' + $script:Url)
}

if ($MyInvocation.InvocationName -ne '.') {
    try {
        Invoke-Launcher
        exit 0
    } catch {
        Write-Host ('Error: ' + $_.Exception.Message) -ForegroundColor Red
        exit 1
    }
}
