# 🏗️ BioScribe AI - Architecture Documentation

## System Overview

BioScribe AI is an industry-grade, full-stack drug discovery platform that transforms protein sequences into potential drug candidates using AI-powered molecular modeling and docking simulation.

```
┌─────────────────────────────────────────────────────────────┐
│                     BioScribe AI Platform                    │
├─────────────────────────────────────────────────────────────┤
│                                                               │
│  ┌──────────────┐      ┌──────────────┐      ┌───────────┐ │
│  │   Frontend   │◄────►│   Backend    │◄────►│ Database  │ │
│  │  (Next.js)   │      │  (FastAPI)   │      │ (MongoDB) │ │
│  └──────────────┘      └──────────────┘      └───────────┘ │
│         │                      │                     │       │
│         │                      │                     │       │
│  ┌──────▼──────┐      ┌───────▼────────┐    ┌──────▼────┐ │
│  │   React     │      │  AI Models     │    │   Redis   │ │
│  │ Components  │      │  & Services    │    │  (Cache)  │ │
│  └─────────────┘      └────────────────┘    └───────────┘ │
│                                                               │
└─────────────────────────────────────────────────────────────┘
```

---

## Technology Stack

### Frontend
- **Framework:** Next.js 15 (React 19)
- **Language:** TypeScript
- **Styling:** TailwindCSS 4
- **UI Components:** Shadcn/UI (Radix UI)
- **Animations:** Framer Motion
- **3D Visualization:** 3Dmol.js
- **Icons:** Lucide React

### Backend
- **Framework:** FastAPI (Python 3.11+)
- **ASGI Server:** Uvicorn
- **Validation:** Pydantic V2
- **Scientific Computing:** NumPy, SciPy, Pandas
- **Bioinformatics:** BioPython
- **Chemistry:** RDKit (optional)
- **AI/ML:** PyTorch, Transformers (optional)

### Database & Caching
- **Primary Database:** MongoDB (with in-memory fallback)
- **Caching:** Redis (optional)
- **Session Storage:** In-memory/MongoDB

### DevOps
- **Containerization:** Docker, Docker Compose
- **Testing:** Pytest, Jest
- **Logging:** Structured JSON logging
- **Monitoring:** Health checks, metrics endpoints

---

## Backend Architecture

### Directory Structure

```
backend/
├── main_real.py              # Production-ready API (recommended)
├── main.py                   # Full-featured API
├── main_demo.py              # Demo version
├── config.py                 # Configuration management
├── models/                   # Core analysis models
│   ├── protein.py           # Protein sequence analysis
│   ├── drug_generator_simple.py
│   ├── docking_simple.py
│   ├── ai_molecular_engine.py
│   └── advanced_docking.py
├── services/                 # Real processing services
│   ├── real_ai_analysis.py
│   ├── real_drug_generation.py
│   ├── real_docking_service.py
│   └── real_interaction_analysis.py
├── database/                 # Database management
│   └── mongodb.py
├── tests/                    # Test suite
│   ├── test_api.py
│   └── __init__.py
├── requirements.txt          # Standard dependencies
├── requirements-minimal.txt  # Minimal dependencies
├── requirements-production.txt # Full production stack
└── .env.example             # Environment template
```

### API Endpoints

#### Core Endpoints
```
GET  /                        # API information
GET  /api/health             # Health check with metrics
GET  /docs                   # OpenAPI documentation
GET  /redoc                  # ReDoc documentation
```

#### Analysis Endpoints
```
POST /api/ai/analyze-protein      # Protein sequence analysis
POST /api/ai/generate-molecules   # Drug candidate generation
```

### Request Flow

```
Client Request
    │
    ▼
Request Logging Middleware
    │
    ▼
CORS Middleware
    │
    ▼
Request Validation (Pydantic)
    │
    ▼
Business Logic
    │
    ├─► Protein Analysis
    │   ├─► Calculate molecular properties
    │   ├─► Predict binding sites
    │   └─► Calculate druggability score
    │
    ├─► Drug Generation
    │   ├─► Generate SMILES structures
    │   ├─► Calculate molecular properties
    │   └─► Estimate binding affinity
    │
    ▼
Response with Metrics
    │
    ▼
Error Handling (if needed)
    │
    ▼
Client Response
```

---

## Frontend Architecture

### Directory Structure

```
bioscribe-ai/
├── src/
│   ├── app/
│   │   ├── page.tsx          # Main application page
│   │   ├── layout.tsx        # Root layout
│   │   └── globals.css       # Global styles
│   ├── components/
│   │   ├── BioScribeWorkflow.tsx  # Main workflow orchestrator
│   │   ├── MolecularViewer.tsx    # 3D visualization
│   │   ├── tabs/
│   │   │   ├── ProteinInputTab.tsx
│   │   │   ├── DrugGenerationTab.tsx
│   │   │   ├── DockingTab.tsx
│   │   │   └── VisualizationTab.tsx
│   │   └── ui/               # Shadcn/UI components
│   └── lib/
│       ├── api.ts            # API client
│       └── utils.ts          # Utilities
├── public/                   # Static assets
└── package.json
```

### Component Hierarchy

```
App
└── BioScribeWorkflow
    ├── Header
    ├── Tabs
    │   ├── ProteinInputTab
    │   │   ├── SequenceInput
    │   │   ├── ExampleProteins
    │   │   └── AnalysisResults
    │   ├── DrugGenerationTab
    │   │   ├── GenerationControls
    │   │   └── CandidatesList
    │   ├── DockingTab
    │   │   ├── DockingControls
    │   │   └── DockingResults
    │   └── VisualizationTab
    │       └── MolecularViewer
    │           ├── ViewControls
    │           └── 3DMolCanvas
    └── Footer
```

---

## Data Flow

### Protein Analysis Flow

```
User Input (FASTA)
    │
    ▼
Frontend Validation
    │
    ▼
API Request (POST /api/ai/analyze-protein)
    │
    ▼
Backend Validation (Pydantic)
    │
    ▼
RealProteinAnalyzer
    ├─► Calculate Molecular Weight
    ├─► Calculate Isoelectric Point
    ├─► Calculate Hydrophobicity (GRAVY)
    ├─► Predict Binding Sites
    │   ├─► ATP binding motifs
    │   ├─► DNA binding motifs
    │   └─► Hydrophobic pockets
    └─► Calculate Druggability Score
    │
    ▼
Response with Analysis Results
    │
    ▼
Frontend Display
```

### Drug Generation Flow

```
Protein Analysis Results
    │
    ▼
API Request (POST /api/ai/generate-molecules)
    │
    ▼
RealDrugGenerator
    ├─► Select Drug Templates (SMILES)
    ├─► Calculate Molecular Properties
    │   ├─► Molecular Weight
    │   ├─► LogP (Lipophilicity)
    │   ├─► TPSA (Polar Surface Area)
    │   ├─► H-bond Donors/Acceptors
    │   └─► QED Score (Drug-likeness)
    ├─► Estimate Binding Affinity
    └─► Rank by QED Score
    │
    ▼
Response with Drug Candidates
    │
    ▼
Frontend Display & Visualization
```

---

## Key Features

### 1. Industry-Grade Error Handling

```python
# Global exception handler
@app.exception_handler(Exception)
async def global_exception_handler(request, exc):
    # Log error with traceback
    # Return structured error response
    # Track error metrics
```

### 2. Request Logging & Metrics

```python
# Request logging middleware
@app.middleware("http")
async def log_requests(request, call_next):
    # Track request count
    # Measure response time
    # Add custom headers (X-Process-Time, X-Request-ID)
```

### 3. Input Validation

```python
class ProteinAnalysisRequest(BaseModel):
    sequence: str = Field(..., min_length=10, max_length=10000)
    
    @validator('sequence')
    def validate_sequence(cls, v):
        # Clean sequence
        # Validate amino acids
        # Return cleaned sequence
```

### 4. Health Monitoring

```python
@app.get("/api/health")
async def health_check():
    return {
        "status": "healthy",
        "uptime_seconds": ...,
        "metrics": {
            "total_requests": ...,
            "total_errors": ...,
            "error_rate": ...
        },
        "services": {
            "protein_analysis": "operational",
            "drug_generation": "operational"
        }
    }
```

### 5. Lifespan Management

```python
@asynccontextmanager
async def lifespan(app: FastAPI):
    # Startup: Initialize services, log startup
    yield
    # Shutdown: Cleanup, log statistics
```

---

## Scientific Accuracy

### Protein Analysis
- **BioPython ProteinAnalysis:** Industry-standard calculations
- **Molecular Weight:** Accurate amino acid weights with peptide bond correction
- **Isoelectric Point:** Based on amino acid composition
- **GRAVY Score:** Kyte-Doolittle hydrophobicity scale
- **Binding Sites:** Motif-based prediction (ATP, DNA, hydrophobic pockets)

### Drug Generation
- **SMILES Templates:** Real drug-like structures
- **Molecular Properties:** Calculated from SMILES structure
- **QED Score:** Quantitative Estimate of Drug-likeness
- **Lipinski Rule of 5:** Drug-likeness filtering
- **Binding Affinity:** Estimated based on molecular properties

---

## Performance Optimization

### Backend
- **Async/Await:** Non-blocking I/O operations
- **Connection Pooling:** MongoDB connection reuse
- **Caching:** Redis for frequently accessed data
- **Lazy Loading:** Models loaded on demand

### Frontend
- **Code Splitting:** Next.js automatic code splitting
- **Image Optimization:** Next.js Image component
- **Lazy Loading:** Dynamic imports for heavy components
- **Memoization:** React.memo for expensive renders

---

## Security Features

### Input Validation
- Pydantic models with strict validation
- Sequence length limits (10-10,000 amino acids)
- Character filtering and sanitization

### CORS Configuration
- Whitelist of allowed origins
- Credentials support
- Controlled methods and headers

### Error Handling
- No sensitive data in error messages (production)
- Structured error responses
- Comprehensive logging

---

## Deployment Architecture

### Development
```
localhost:3000 (Frontend)
    │
    ▼
localhost:8000 (Backend)
    │
    ▼
localhost:27017 (MongoDB - optional)
```

### Production (Docker)
```
nginx:80 (Reverse Proxy)
    │
    ├─► frontend:3000 (Next.js)
    │
    └─► backend:8000 (FastAPI)
            │
            ├─► mongodb:27017
            └─► redis:6379
```

### Cloud Deployment
```
Vercel (Frontend)
    │
    ▼
Railway/Render (Backend)
    │
    ▼
MongoDB Atlas (Database)
```

---

## Testing Strategy

### Backend Tests
- **Unit Tests:** Individual function testing
- **Integration Tests:** API endpoint testing
- **Performance Tests:** Response time validation
- **Validation Tests:** Input validation testing

### Frontend Tests
- **Component Tests:** React component testing
- **Integration Tests:** User flow testing
- **E2E Tests:** Full workflow testing

---

## Monitoring & Observability

### Metrics
- Request count
- Error count and rate
- Response times
- Uptime tracking

### Logging
- Structured JSON logging
- Request/response logging
- Error tracking with stack traces
- Performance metrics

### Health Checks
- Service status
- Database connectivity
- API responsiveness

---

## Scalability Considerations

### Horizontal Scaling
- Stateless API design
- Session storage in database
- Load balancer ready

### Vertical Scaling
- Async processing
- Connection pooling
- Efficient algorithms

### Caching Strategy
- Redis for frequent queries
- In-memory caching for static data
- CDN for frontend assets

---

**Built with industry best practices for production-grade biotech applications**
