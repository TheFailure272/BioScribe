# ============================================================================
# BioScribe AI - Automated Setup Script (Windows PowerShell)
# ============================================================================

Write-Host "🧬 BioScribe AI - Automated Setup" -ForegroundColor Cyan
Write-Host "=================================" -ForegroundColor Cyan
Write-Host ""

# Check Python installation
Write-Host "Checking Python installation..." -ForegroundColor Yellow
try {
    $pythonVersion = python --version 2>&1
    Write-Host "✅ Found: $pythonVersion" -ForegroundColor Green
} catch {
    Write-Host "❌ Python not found!" -ForegroundColor Red
    Write-Host "Please install Python 3.11+ from https://www.python.org/downloads/" -ForegroundColor Red
    Write-Host "Make sure to check 'Add Python to PATH' during installation" -ForegroundColor Yellow
    exit 1
}

# Check Node.js installation
Write-Host "Checking Node.js installation..." -ForegroundColor Yellow
try {
    $nodeVersion = node --version 2>&1
    Write-Host "✅ Found Node.js: $nodeVersion" -ForegroundColor Green
} catch {
    Write-Host "❌ Node.js not found!" -ForegroundColor Red
    Write-Host "Please install Node.js 18+ from https://nodejs.org/" -ForegroundColor Red
    exit 1
}

Write-Host ""
Write-Host "=================================" -ForegroundColor Cyan
Write-Host "Setting up Backend..." -ForegroundColor Cyan
Write-Host "=================================" -ForegroundColor Cyan

# Navigate to backend
Set-Location -Path "backend"

# Create virtual environment
Write-Host "Creating Python virtual environment..." -ForegroundColor Yellow
if (Test-Path "venv") {
    Write-Host "Virtual environment already exists, skipping..." -ForegroundColor Yellow
} else {
    python -m venv venv
    Write-Host "✅ Virtual environment created" -ForegroundColor Green
}

# Activate virtual environment
Write-Host "Activating virtual environment..." -ForegroundColor Yellow
& ".\venv\Scripts\Activate.ps1"

# Upgrade pip
Write-Host "Upgrading pip..." -ForegroundColor Yellow
python -m pip install --upgrade pip

# Ask user which requirements to install
Write-Host ""
Write-Host "Choose installation type:" -ForegroundColor Cyan
Write-Host "1. Minimal (Fast, ~500MB) - Recommended for quick start" -ForegroundColor Green
Write-Host "2. Standard (Moderate, ~2GB) - Full features" -ForegroundColor Yellow
Write-Host "3. Production (Heavy, ~5GB) - All features + AI models" -ForegroundColor Red
Write-Host ""
$choice = Read-Host "Enter choice (1-3, default: 1)"

if ([string]::IsNullOrWhiteSpace($choice)) {
    $choice = "1"
}

switch ($choice) {
    "1" {
        $reqFile = "requirements-minimal.txt"
        Write-Host "Installing minimal dependencies..." -ForegroundColor Green
    }
    "2" {
        $reqFile = "requirements.txt"
        Write-Host "Installing standard dependencies..." -ForegroundColor Yellow
    }
    "3" {
        $reqFile = "requirements-production.txt"
        Write-Host "Installing production dependencies (this may take a while)..." -ForegroundColor Red
    }
    default {
        $reqFile = "requirements-minimal.txt"
        Write-Host "Installing minimal dependencies..." -ForegroundColor Green
    }
}

# Install dependencies
if (Test-Path $reqFile) {
    pip install -r $reqFile
    Write-Host "✅ Dependencies installed" -ForegroundColor Green
} else {
    Write-Host "⚠️  $reqFile not found, using requirements.txt" -ForegroundColor Yellow
    pip install -r requirements.txt
}

# Create .env file if it doesn't exist
if (-not (Test-Path ".env")) {
    Write-Host "Creating .env file..." -ForegroundColor Yellow
    Copy-Item ".env.example" ".env" -ErrorAction SilentlyContinue
    Write-Host "✅ .env file created (please update with your settings)" -ForegroundColor Green
}

# Return to root
Set-Location -Path ".."

Write-Host ""
Write-Host "=================================" -ForegroundColor Cyan
Write-Host "Setting up Frontend..." -ForegroundColor Cyan
Write-Host "=================================" -ForegroundColor Cyan

# Navigate to frontend
Set-Location -Path "bioscribe-ai"

# Install npm dependencies
Write-Host "Installing npm dependencies..." -ForegroundColor Yellow
npm install

if ($LASTEXITCODE -eq 0) {
    Write-Host "✅ Frontend dependencies installed" -ForegroundColor Green
} else {
    Write-Host "⚠️  Frontend installation had warnings (this is usually okay)" -ForegroundColor Yellow
}

# Return to root
Set-Location -Path ".."

Write-Host ""
Write-Host "=================================" -ForegroundColor Green
Write-Host "✅ Setup Complete!" -ForegroundColor Green
Write-Host "=================================" -ForegroundColor Green
Write-Host ""
Write-Host "To start the application:" -ForegroundColor Cyan
Write-Host ""
Write-Host "Backend:" -ForegroundColor Yellow
Write-Host "  cd backend" -ForegroundColor White
Write-Host "  .\venv\Scripts\Activate.ps1" -ForegroundColor White
Write-Host "  python main_real.py" -ForegroundColor White
Write-Host "  → http://localhost:8000" -ForegroundColor Gray
Write-Host ""
Write-Host "Frontend:" -ForegroundColor Yellow
Write-Host "  cd bioscribe-ai" -ForegroundColor White
Write-Host "  npm run dev" -ForegroundColor White
Write-Host "  → http://localhost:3000" -ForegroundColor Gray
Write-Host ""
Write-Host "Or use the start scripts:" -ForegroundColor Cyan
Write-Host "  .\start-backend.ps1" -ForegroundColor White
Write-Host "  .\start-frontend.ps1" -ForegroundColor White
Write-Host ""
Write-Host "🧬 Happy drug discovering!" -ForegroundColor Magenta
