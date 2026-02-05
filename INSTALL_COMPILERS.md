# 🔧 Installing C++ Build Tools for BioScribe AI

## Why You Need This

Many Python scientific packages (NumPy, SciPy, BioPython, RDKit, PyTorch) require C++ compilation on Windows. You need Microsoft Visual C++ Build Tools.

## Step 1: Install Microsoft C++ Build Tools

### Option A: Visual Studio Build Tools (Recommended)

1. **Download:** https://visualstudio.microsoft.com/visual-cpp-build-tools/

2. **Run the installer**

3. **Select "Desktop development with C++"**
   - ✅ MSVC v143 - VS 2022 C++ x64/x86 build tools
   - ✅ Windows 10/11 SDK
   - ✅ C++ CMake tools for Windows

4. **Install** (This will take 5-10 GB of space and 15-30 minutes)

### Option B: Visual Studio Community (Full IDE)

1. **Download:** https://visualstudio.microsoft.com/downloads/

2. **Install with:**
   - ✅ Desktop development with C++
   - ✅ Python development (optional but helpful)

## Step 2: Restart Your Computer

After installation, **restart your computer** to ensure all environment variables are set.

## Step 3: Verify Installation

Open PowerShell and run:

```powershell
# Check if cl.exe (C++ compiler) is available
where cl

# If not found, you may need to run from "Developer Command Prompt for VS"
```

## Step 4: Install Python Dependencies

After restarting, run:

```powershell
cd E:\Bioscribe\CascadeProjects\windsurf-project\backend

# Activate virtual environment
.\venv\Scripts\Activate.ps1

# Install production dependencies
pip install -r requirements-production.txt

# Or standard dependencies
pip install -r requirements.txt
```

## Alternative: Use Pre-built Wheels

If you don't want to install Build Tools, you can use pre-built wheels:

```powershell
# Install from pre-built wheels
pip install --only-binary :all: numpy scipy pandas biopython

# For PyTorch (CPU version)
pip install torch --index-url https://download.pytorch.org/whl/cpu

# For other packages
pip install -r requirements.txt
```

## Troubleshooting

### "Microsoft Visual C++ 14.0 or greater is required"

This means Build Tools aren't installed or not found. Solutions:

1. Install Visual Studio Build Tools (see above)
2. Restart your computer
3. Try installing from Developer Command Prompt

### "error: command 'cl.exe' failed"

1. Make sure you restarted after installing Build Tools
2. Try running from "Developer Command Prompt for VS 2022"
3. Check if cl.exe is in PATH: `where cl`

### Still Having Issues?

Use Anaconda/Miniconda instead:

```powershell
# Download Miniconda: https://docs.conda.io/en/latest/miniconda.html

# Create environment
conda create -n bioscribe python=3.11

# Activate
conda activate bioscribe

# Install packages
conda install -c conda-forge numpy scipy pandas biopython rdkit
pip install fastapi uvicorn pydantic
```

---

## Quick Links

- **Build Tools:** https://visualstudio.microsoft.com/visual-cpp-build-tools/
- **Visual Studio:** https://visualstudio.microsoft.com/downloads/
- **Miniconda:** https://docs.conda.io/en/latest/miniconda.html
