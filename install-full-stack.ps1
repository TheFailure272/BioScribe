# ============================================================================
# BioScribe AI - Full Stack Installation with Build Tools Check
# ============================================================================

Write-Host "🧬 BioScribe AI - Full Production Installation" -ForegroundColor Cyan
Write-Host "=" * 60 -ForegroundColor Cyan
Write-Host ""

# Check for Visual C++ Build Tools
Write-Host "Checking for Microsoft Visual C++ Build Tools..." -ForegroundColor Yellow

$vcVarsPath = "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
$vsVarsPath = "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat"

$buildToolsFound = $false

if (Test-Path $vcVarsPath) {
    Write-Host "✅ Found Visual Studio Build Tools 2022" -ForegroundColor Green
    $buildToolsFound = $true
} elseif (Test-Path $vsVarsPath) {
    Write-Host "✅ Found Visual Studio Community 2022" -ForegroundColor Green
    $buildToolsFound = $true
} else {
    Write-Host "❌ Microsoft Visual C++ Build Tools not found!" -ForegroundColor Red
    Write-Host ""
    Write-Host "You need to install Microsoft C++ Build Tools first:" -ForegroundColor Yellow
    Write-Host ""
    Write-Host "1. Download from: https://visualstudio.microsoft.com/visual-cpp-build-tools/" -ForegroundColor White
    Write-Host "2. Run the installer" -ForegroundColor White
    Write-Host "3. Select 'Desktop development with C++'" -ForegroundColor White
    Write-Host "4. Install (requires 5-10 GB, takes 15-30 minutes)" -ForegroundColor White
    Write-Host "5. RESTART YOUR COMPUTER" -ForegroundColor Red
    Write-Host "6. Run this script again" -ForegroundColor White
    Write-Host ""
    
    $response = Read-Host "Do you want to open the download page now? (Y/N)"
    if ($response -eq "Y" -or $response -eq "y") {
        Start-Process "https://visualstudio.microsoft.com/visual-cpp-build-tools/"
    }
    
    Write-Host ""
    Write-Host "Alternatively, you can install minimal dependencies without Build Tools:" -ForegroundColor Yellow
    Write-Host "  pip install -r requirements-windows.txt" -ForegroundColor White
    Write-Host ""
    
    exit 1
}

Write-Host ""
Write-Host "=" * 60 -ForegroundColor Cyan
Write-Host "Installing Backend Dependencies..." -ForegroundColor Cyan
Write-Host "=" * 60 -ForegroundColor Cyan

# Navigate to backend
Set-Location -Path "backend"

# Create virtual environment if it doesn't exist
if (-not (Test-Path "venv")) {
    Write-Host "Creating Python virtual environment..." -ForegroundColor Yellow
    python -m venv venv
    Write-Host "✅ Virtual environment created" -ForegroundColor Green
}

# Activate virtual environment
Write-Host "Activating virtual environment..." -ForegroundColor Yellow
& ".\venv\Scripts\Activate.ps1"

# Upgrade pip
Write-Host "Upgrading pip..." -ForegroundColor Yellow
python -m pip install --upgrade pip

# Ask which requirements to install
Write-Host ""
Write-Host "Choose installation type:" -ForegroundColor Cyan
Write-Host "1. Standard (requirements.txt) - ~2GB, recommended" -ForegroundColor Green
Write-Host "2. Production (requirements-production.txt) - ~5GB, all features" -ForegroundColor Yellow
Write-Host "3. Minimal (requirements-windows.txt) - ~500MB, basic features" -ForegroundColor White
Write-Host ""
$choice = Read-Host "Enter choice (1-3, default: 1)"

if ([string]::IsNullOrWhiteSpace($choice)) {
    $choice = "1"
}

switch ($choice) {
    "1" {
        $reqFile = "requirements.txt"
        Write-Host "Installing standard dependencies (this may take 10-20 minutes)..." -ForegroundColor Green
    }
    "2" {
        $reqFile = "requirements-production.txt"
        Write-Host "Installing production dependencies (this may take 20-40 minutes)..." -ForegroundColor Yellow
    }
    "3" {
        $reqFile = "requirements-windows.txt"
        Write-Host "Installing minimal dependencies..." -ForegroundColor White
    }
    default {
        $reqFile = "requirements.txt"
        Write-Host "Installing standard dependencies..." -ForegroundColor Green
    }
}

Write-Host ""
Write-Host "Installing from $reqFile..." -ForegroundColor Yellow
Write-Host "This will download and compile packages. Please be patient..." -ForegroundColor Yellow
Write-Host ""

# Install dependencies
pip install -r $reqFile

if ($LASTEXITCODE -eq 0) {
    Write-Host ""
    Write-Host "=" * 60 -ForegroundColor Green
    Write-Host "✅ Installation Complete!" -ForegroundColor Green
    Write-Host "=" * 60 -ForegroundColor Green
    Write-Host ""
    Write-Host "To start the backend:" -ForegroundColor Cyan
    Write-Host "  cd backend" -ForegroundColor White
    Write-Host "  .\venv\Scripts\Activate.ps1" -ForegroundColor White
    Write-Host "  python main_real.py" -ForegroundColor White
    Write-Host ""
    Write-Host "Or use the start script:" -ForegroundColor Cyan
    Write-Host "  .\start-backend.ps1" -ForegroundColor White
    Write-Host ""
} else {
    Write-Host ""
    Write-Host "=" * 60 -ForegroundColor Red
    Write-Host "❌ Installation Failed!" -ForegroundColor Red
    Write-Host "=" * 60 -ForegroundColor Red
    Write-Host ""
    Write-Host "Common solutions:" -ForegroundColor Yellow
    Write-Host "1. Make sure you restarted your computer after installing Build Tools" -ForegroundColor White
    Write-Host "2. Try running from 'Developer Command Prompt for VS 2022'" -ForegroundColor White
    Write-Host "3. Install Miniconda and use conda instead" -ForegroundColor White
    Write-Host ""
    Write-Host "See INSTALL_COMPILERS.md for detailed troubleshooting" -ForegroundColor Yellow
}

# Return to root
Set-Location -Path ".."
