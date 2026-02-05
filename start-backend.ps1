# ============================================================================
# BioScribe AI - Start Backend Server
# ============================================================================

Write-Host "🧬 Starting BioScribe AI Backend..." -ForegroundColor Cyan

# Navigate to backend
Set-Location -Path "backend"

# Activate virtual environment
if (Test-Path "venv\Scripts\Activate.ps1") {
    Write-Host "Activating virtual environment..." -ForegroundColor Yellow
    & ".\venv\Scripts\Activate.ps1"
} else {
    Write-Host "⚠️  Virtual environment not found. Run setup.ps1 first!" -ForegroundColor Red
    exit 1
}

# Start the server
Write-Host "Starting FastAPI server..." -ForegroundColor Green
Write-Host "API will be available at: http://localhost:8000" -ForegroundColor Cyan
Write-Host "API Documentation: http://localhost:8000/docs" -ForegroundColor Cyan
Write-Host ""
Write-Host "Press Ctrl+C to stop the server" -ForegroundColor Yellow
Write-Host ""

python main_real.py
