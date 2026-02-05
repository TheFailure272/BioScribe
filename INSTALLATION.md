# 🧬 BioScribe AI - Installation Guide

## Quick Installation (Recommended)

### 1. Install Prerequisites

**Python 3.11+**
- Download: https://www.python.org/downloads/
- ✅ During installation, check "Add Python to PATH"

**Node.js 18+**
- Download: https://nodejs.org/

### 2. Run Automated Setup

```powershell
# Clone or navigate to project directory
cd path\to\windsurf-project

# Run automated setup
.\setup.ps1
```

The script will:
- ✅ Verify Python and Node.js installation
- ✅ Create Python virtual environment
- ✅ Install backend dependencies
- ✅ Install frontend dependencies
- ✅ Create configuration files

### 3. Start the Application

**Option A: Use start scripts (Easiest)**
```powershell
# Terminal 1 - Backend
.\start-backend.ps1

# Terminal 2 - Frontend
.\start-frontend.ps1
```

**Option B: Manual start**
```powershell
# Terminal 1 - Backend
cd backend
.\venv\Scripts\Activate.ps1
python main_real.py

# Terminal 2 - Frontend
cd bioscribe-ai
npm run dev
```

### 4. Access the Application

- **Frontend:** http://localhost:3000
- **Backend API:** http://localhost:8000
- **API Docs:** http://localhost:8000/docs

---

## Manual Installation

### Backend Setup

```powershell
# Navigate to backend
cd backend

# Create virtual environment
python -m venv venv

# Activate virtual environment
.\venv\Scripts\Activate.ps1

# Upgrade pip
python -m pip install --upgrade pip

# Install dependencies (choose one)
pip install -r requirements-minimal.txt    # Quick start (recommended)
pip install -r requirements.txt            # Standard features
pip install -r requirements-production.txt # Full production

# Create environment file
copy .env.example .env

# Start server
python main_real.py
```

### Frontend Setup

```powershell
# Navigate to frontend
cd bioscribe-ai

# Install dependencies
npm install

# Start development server
npm run dev
```

---

## Dependency Options

### Minimal (Fast - ~500MB)
**File:** `requirements-minimal.txt`
- Core FastAPI framework
- BioPython for protein analysis
- Basic scientific computing
- **Install time:** ~5 minutes
- **Use for:** Quick testing, development

```powershell
pip install -r requirements-minimal.txt
```

### Standard (Moderate - ~2GB)
**File:** `requirements.txt`
- Everything in minimal
- RDKit for molecular calculations
- Enhanced features
- **Install time:** ~15 minutes
- **Use for:** Full development

```powershell
pip install -r requirements.txt
```

### Production (Heavy - ~5GB)
**File:** `requirements-production.txt`
- Everything in standard
- PyTorch and AI models
- Advanced molecular dynamics
- Full monitoring and security
- **Install time:** ~30 minutes
- **Use for:** Production deployment

```powershell
pip install -r requirements-production.txt
```

---

## Troubleshooting

### Python Not Found

**Check if Python is installed:**
```powershell
python --version
```

**If not found, add to PATH:**
```powershell
# Find Python installation path
where python

# Add to PATH (temporary)
$env:Path += ";C:\Python311;C:\Python311\Scripts"

# Or add permanently via System Environment Variables
```

### Virtual Environment Activation Fails

**Enable script execution:**
```powershell
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
```

**Then activate:**
```powershell
.\venv\Scripts\Activate.ps1
```

### Dependency Installation Errors

**RDKit Installation Issues:**
```powershell
# Option 1: Use conda (easier)
conda install -c conda-forge rdkit

# Option 2: Use pre-built wheel
pip install rdkit-pypi
```

**PyTorch Installation:**
```powershell
# CPU version (lighter)
pip install torch --index-url https://download.pytorch.org/whl/cpu

# GPU version (CUDA 11.8)
pip install torch --index-url https://download.pytorch.org/whl/cu118
```

**NumPy/SciPy Issues:**
```powershell
# Install Microsoft Visual C++ Build Tools
# Download: https://visualstudio.microsoft.com/visual-cpp-build-tools/

# Or use pre-built wheels
pip install --only-binary :all: numpy scipy
```

### Port Already in Use

**Backend (port 8000):**
```powershell
# Find process using port 8000
netstat -ano | findstr :8000

# Kill process (replace PID)
taskkill /PID <PID> /F

# Or use different port
uvicorn main_real:app --port 8001
```

**Frontend (port 3000):**
```powershell
# Use different port
npm run dev -- -p 3001
```

### MongoDB Connection Issues

The application works without MongoDB (uses in-memory storage).

**To install MongoDB locally:**
- Download: https://www.mongodb.com/try/download/community
- Or use MongoDB Atlas (cloud): https://www.mongodb.com/cloud/atlas

**Update .env file:**
```
MONGODB_URL=mongodb://localhost:27017
```

---

## Docker Installation (Alternative)

### Prerequisites
- Docker Desktop for Windows
- Download: https://www.docker.com/products/docker-desktop

### Quick Start with Docker

```powershell
# Build and start all services
docker-compose up --build

# Access application
# Frontend: http://localhost:3000
# Backend: http://localhost:8000
# MongoDB: localhost:27017
```

### Docker Commands

```powershell
# Start services
docker-compose up

# Start in background
docker-compose up -d

# Stop services
docker-compose down

# View logs
docker-compose logs -f

# Rebuild after changes
docker-compose up --build
```

---

## Verification

### Test Backend

```powershell
# Check health endpoint
curl http://localhost:8000/api/health

# Or visit in browser
http://localhost:8000/docs
```

### Test Frontend

```powershell
# Visit in browser
http://localhost:3000
```

### Run Tests

```powershell
cd backend
pytest tests/ -v
```

---

## Next Steps

1. ✅ Complete installation
2. ✅ Start both servers
3. ✅ Open http://localhost:3000
4. 🧬 Try example proteins
5. 🧪 Generate drug candidates
6. 📊 View 3D visualizations

---

## Production Deployment

### Docker Deployment (Recommended)
```powershell
docker-compose -f docker-compose.yml up -d
```

### Cloud Deployment
- **Frontend:** Vercel, Netlify
- **Backend:** Railway, Render, AWS
- **Database:** MongoDB Atlas

See `deploy.md` for detailed deployment instructions.

---

## Getting Help

- **Setup Issues:** Check SETUP_GUIDE.md
- **API Documentation:** http://localhost:8000/docs
- **GitHub Issues:** Open an issue for bugs
- **Examples:** See example proteins in the UI

---

**Built with ❤️ for the biotech community**

🧬 Transform protein sequences into life-saving drugs through AI
