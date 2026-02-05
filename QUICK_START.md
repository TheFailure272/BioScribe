# 🚀 BioScribe AI - Quick Start Guide

## ⚡ 5-Minute Setup

### Step 1: Install Python
1. Download Python 3.11+ from https://www.python.org/downloads/
2. ✅ **IMPORTANT:** Check "Add Python to PATH" during installation
3. Verify: Open PowerShell and run `python --version`

### Step 2: Run Automated Setup
```powershell
# Navigate to project
cd path\to\windsurf-project

# Run setup script
.\setup.ps1
```

**Choose option 1 (Minimal)** when prompted for fastest setup.

### Step 3: Start the Application

**Terminal 1 - Backend:**
```powershell
.\start-backend.ps1
```

**Terminal 2 - Frontend:**
```powershell
.\start-frontend.ps1
```

### Step 4: Open Application
🌐 **http://localhost:3000**

---

## 🎯 First Use

1. **Select Example Protein** (e.g., HIV-1 Protease)
2. Click **"Analyze Protein"**
3. Navigate to **"Drug Generation"** tab
4. Click **"Generate Molecules"**
5. View results in **"Docking"** and **"Visualization"** tabs

---

## 🆘 Troubleshooting

### Python Not Found?
```powershell
# Check if Python is installed
python --version

# If not found, add to PATH or reinstall Python
```

### Port Already in Use?
```powershell
# Backend (port 8000)
# Kill process or use different port:
uvicorn main_real:app --port 8001

# Frontend (port 3000)
npm run dev -- -p 3001
```

### Dependencies Failed?
```powershell
# Try minimal installation
cd backend
pip install -r requirements-minimal.txt
```

---

## 📚 Full Documentation

- **Installation:** See `INSTALLATION.md`
- **Setup Guide:** See `SETUP_GUIDE.md`
- **Architecture:** See `ARCHITECTURE.md`
- **Enhancements:** See `ENHANCEMENTS.md`

---

## 🧬 What You Get

✅ **Industry-Grade Backend**
- FastAPI with structured logging
- Global error handling
- Request/response tracking
- Health monitoring
- Comprehensive validation

✅ **Modern Frontend**
- Next.js 15 with React 19
- Beautiful UI with TailwindCSS
- Interactive 3D visualization
- Real-time updates

✅ **Production-Ready**
- Docker support
- Automated testing
- Multiple deployment options
- Comprehensive documentation

---

**Ready to discover drugs! 🧬💊**
