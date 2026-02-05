# 🧬 BioScribe AI - Complete Setup Guide

## Prerequisites Installation

### 1. Install Python 3.11+
**Download:** https://www.python.org/downloads/

**Windows Installation:**
- ✅ Check "Add Python to PATH"
- ✅ Check "Install pip"
- ✅ Choose "Customize installation" → Enable all optional features

**Verify Installation:**
```powershell
python --version
# Should show: Python 3.11.x or higher
```

### 2. Install Node.js 18+
**Download:** https://nodejs.org/

**Verify Installation:**
```powershell
node --version
npm --version
```

---

## Quick Start (Recommended)

### Option 1: Automated Setup (Windows)
```powershell
# Run the automated setup script
.\setup.ps1
```

### Option 2: Manual Setup

#### Backend Setup
```powershell
# Navigate to backend
cd backend

# Create virtual environment
python -m venv venv

# Activate virtual environment
.\venv\Scripts\Activate.ps1

# Install dependencies (choose one)
pip install -r requirements-minimal.txt    # Quick start (recommended)
# OR
pip install -r requirements.txt            # Full features
# OR
pip install -r requirements-production.txt # Production-grade (heavy)

# Copy environment file
copy .env.example .env

# Start backend server
python main_real.py
# OR
uvicorn main_real:app --reload --host 0.0.0.0 --port 8000
```

Backend will be available at: **http://localhost:8000**

#### Frontend Setup
```powershell
# Navigate to frontend
cd bioscribe-ai

# Install dependencies
npm install

# Start development server
npm run dev
```

Frontend will be available at: **http://localhost:3000**

---

## Dependency Breakdown

### Minimal Setup (Fast, ~500MB)
**File:** `requirements-minimal.txt`
- Core FastAPI framework
- BioPython for protein analysis
- Basic scientific computing (NumPy, SciPy, Pandas)
- MongoDB support
- **Use for:** Quick testing, development

### Standard Setup (Moderate, ~2GB)
**File:** `requirements.txt`
- Everything in minimal
- RDKit for molecular calculations
- Basic AI models
- Enhanced docking capabilities
- **Use for:** Full feature development

### Production Setup (Heavy, ~5GB)
**File:** `requirements-production.txt`
- Everything in standard
- Advanced AI models (PyTorch, Transformers)
- Molecular dynamics (OpenMM)
- Full monitoring and logging
- Security features
- **Use for:** Production deployment

---

## Troubleshooting

### Python Not Found
```powershell
# Check if Python is in PATH
$env:Path

# Add Python to PATH manually (replace with your Python path)
$env:Path += ";C:\Python311;C:\Python311\Scripts"
```

### Virtual Environment Issues
```powershell
# If activation fails, enable script execution
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser

# Then activate again
.\venv\Scripts\Activate.ps1
```

### Dependency Installation Errors

**RDKit Installation (if needed):**
```powershell
# Use conda for easier RDKit installation
conda install -c conda-forge rdkit
```

**PyTorch Installation (if needed):**
```powershell
# CPU version (lighter)
pip install torch --index-url https://download.pytorch.org/whl/cpu

# GPU version (if you have CUDA)
pip install torch --index-url https://download.pytorch.org/whl/cu118
```

### MongoDB Not Available
The application works without MongoDB (uses in-memory storage as fallback).

To install MongoDB locally:
- Download: https://www.mongodb.com/try/download/community
- Or use MongoDB Atlas (cloud): https://www.mongodb.com/cloud/atlas

---

## Project Structure

```
windsurf-project/
├── backend/
│   ├── models/              # AI and analysis models
│   │   ├── protein.py       # Protein analysis
│   │   ├── drug_generator_simple.py
│   │   ├── docking_simple.py
│   │   ├── ai_molecular_engine.py
│   │   └── advanced_docking.py
│   ├── services/            # Real processing services
│   │   ├── real_ai_analysis.py
│   │   ├── real_drug_generation.py
│   │   ├── real_docking_service.py
│   │   └── real_interaction_analysis.py
│   ├── database/            # Database management
│   │   └── mongodb.py
│   ├── main.py              # Main API (full features)
│   ├── main_real.py         # Real processing API
│   ├── main_demo.py         # Demo version
│   ├── requirements.txt     # Dependencies
│   └── .env.example         # Environment template
├── bioscribe-ai/            # Next.js frontend
│   ├── src/
│   │   ├── app/
│   │   └── components/
│   └── package.json
└── README.md
```

---

## Available Backend Versions

### 1. `main_real.py` (Recommended)
- Real calculations without heavy dependencies
- Fast startup
- Production-ready
- **Use:** `python main_real.py`

### 2. `main.py` (Full Features)
- All AI models
- Advanced docking
- Complete feature set
- **Use:** `python main.py`

### 3. `main_demo.py` (Demo)
- Quick demo version
- Minimal dependencies
- **Use:** `python main_demo.py`

---

## Testing the Setup

### 1. Test Backend
```powershell
# Visit API documentation
http://localhost:8000/docs

# Test health endpoint
curl http://localhost:8000/api/health
```

### 2. Test Frontend
```powershell
# Visit application
http://localhost:3000
```

### 3. Run Tests (if available)
```powershell
cd backend
pytest
```

---

## Next Steps

1. ✅ Install Python and Node.js
2. ✅ Run setup scripts or manual installation
3. ✅ Start backend server
4. ✅ Start frontend application
5. ✅ Open http://localhost:3000
6. 🧬 Start discovering drugs!

---

## Production Deployment

### Docker Deployment (Recommended)
```powershell
# Build and run with Docker Compose
docker-compose up --build
```

### Manual Deployment
- **Frontend:** Deploy to Vercel, Netlify
- **Backend:** Deploy to Railway, Render, AWS
- **Database:** MongoDB Atlas

See `deploy.md` for detailed deployment instructions.

---

## Support

- **Issues:** Open a GitHub issue
- **Documentation:** Check README.md and API docs at /docs
- **Examples:** See example proteins in the UI

---

**Built with ❤️ for the biotech community**
