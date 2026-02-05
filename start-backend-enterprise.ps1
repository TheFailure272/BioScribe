# BioScribe AI - Enterprise Edition Startup Script
# Starts the full-featured enterprise backend

Write-Host "============================================" -ForegroundColor Cyan
Write-Host "🧬 BioScribe AI - ENTERPRISE EDITION" -ForegroundColor Green
Write-Host "============================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "Starting Enterprise Backend..." -ForegroundColor Yellow
Write-Host "All Features: ENABLED" -ForegroundColor Green
Write-Host ""

# Change to backend directory
Set-Location -Path "backend"

# Activate virtual environment if it exists
if (Test-Path ".venv\Scripts\Activate.ps1") {
    Write-Host "✓ Activating virtual environment..." -ForegroundColor Green
    & .venv\Scripts\Activate.ps1
} else {
    Write-Host "⚠ Virtual environment not found. Using global Python..." -ForegroundColor Yellow
}

Write-Host ""
Write-Host "Enterprise Features:" -ForegroundColor Cyan
Write-Host "  ✓ AI-Powered Target Discovery" -ForegroundColor Green
Write-Host "  ✓ Novel Molecule Generation" -ForegroundColor Green
Write-Host "  ✓ Drug Combination Prediction" -ForegroundColor Green
Write-Host "  ✓ Patient Stratification" -ForegroundColor Green
Write-Host "  ✓ Clinical Trial Optimization" -ForegroundColor Green
Write-Host "  ✓ Molecular Dynamics Simulation" -ForegroundColor Green
Write-Host "  ✓ RNA Aptamer Design" -ForegroundColor Green
Write-Host "  ✓ CRISPR Guide Design" -ForegroundColor Green
Write-Host "  ✓ mRNA Therapeutic Design" -ForegroundColor Green
Write-Host "  ✓ Blockchain Recording" -ForegroundColor Green
Write-Host "  ✓ FAIR Data Principles" -ForegroundColor Green
Write-Host "  ✓ Causal AI Validation" -ForegroundColor Green
Write-Host ""
Write-Host "API Documentation: http://localhost:8000/docs" -ForegroundColor Cyan
Write-Host "Health Check: http://localhost:8000/api/health" -ForegroundColor Cyan
Write-Host ""
Write-Host "Press Ctrl+C to stop the server" -ForegroundColor Yellow
Write-Host "============================================" -ForegroundColor Cyan
Write-Host ""

# Start the enterprise backend
python main_enterprise.py
