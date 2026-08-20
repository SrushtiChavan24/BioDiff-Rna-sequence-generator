# RNA Genomics Intelligence Platform

A full-stack RNA analytics system for sequence profiling, secondary structure prediction, molecular similarity analysis, and ranked therapeutic candidate generation.

This repository contains the local application code for `RNA_Project_local`, including the backend API server and the frontend dashboard. The remaining backend analysis experiments and Colab-based development are captured in the notebook `project_1 (5).ipynb`.

---

## Project Overview

The platform is designed to support RNA research workflows by converting RNA sequence inputs into structured scientific outputs:

- Sequence validation and normalization
- Secondary structure prediction using a pure Python Nussinov-style folding engine
- Thermodynamic stability scoring and structured candidate ranking
- Molecular similarity analysis across sequence/structure pairs
- Structural family detection and conformational projection
- Interactive dashboard visualization and exportable reports

The system is split into two main layers:

1. **Backend**: FastAPI server, analysis services, SQLite persistence, PDF export, and scientific interpretation.
2. **Frontend**: Next.js dashboard with genome analytics pages, interactive plots, and export workflow.

---

## Repository Structure

```text
RNA_Project_local/
├── backend/                # FastAPI backend server and analysis engine
│   ├── run.py              # Application entrypoint
│   ├── requirements.txt    # Python dependencies
│   └── app/
│       ├── main.py         # FastAPI app, router registration, CORS, exceptions
│       ├── config.py       # Environment-driven settings and directories
│       ├── routes/         # API endpoints for analysis, structure, candidates, export, health
│       ├── schemas/        # Pydantic payloads and response models
│       ├── services/       # Core RNA analysis, structure prediction, similarity, clustering
│       ├── interpretations/# Domain-specific narrative explanation engine
│       ├── visualization/  # SVG and plot/chart generators
│       ├── database/       # Async SQLite connection and ORM models
│       └── utils/          # Sequence validators and helpers
├── frontend/               # React / Next.js application for the analytics dashboard
│   ├── package.json        # Frontend dependencies and scripts
│   └── src/
│       ├── app/            # Pages for analysis, structure, candidates, similarity, export
│       ├── components/     # Reusable UI and visualization components
│       ├── lib/            # API client, app-wide utilities, context, constants
│       └── types/          # TypeScript schema definitions matching backend models
├── outputs/                # Generated artifacts and persisted data
├── project_1 (5).ipynb     # Colab/Notebook backend experiments and development notes
└── structure.txt           # Project structure reference
```

---

## System Architecture

### High level

- **Frontend**: React + Next.js dashboard connects to the backend API using a shared REST contract.
- **Backend**: FastAPI exposes analysis endpoints, handles validation, runs RNA intelligence services, and persists results.
- **Database**: SQLite via SQLAlchemy and `aiosqlite` stores analysis records and candidate metadata.
- **Exports**: Reports can be generated as PDF, JSON, or FASTA artifacts.

### Working architecture

1. User submits an RNA sequence from the browser.
2. Frontend sends a POST request to the backend analysis endpoint.
3. Backend validates the sequence, predicts secondary structure, computes stability and similarity metrics, and stores the result.
4. Frontend displays interpreted results, visualizations, and candidate insights.
5. Users can compare multiple candidates, detect structural families, and export results.

---

## Frontend Details

The frontend is built with:

- **Next.js 16**
- **React 19**
- **TypeScript**
- **Tailwind CSS** / ShadCN UI primitives
- **Plotly** for interactive charts and landscape visualization

### Core frontend features

- Sequence analysis page with input form and result summary
- Candidate registry and ranked candidate view
- Structure visualization page with dot-bracket and circular structure display
- Molecular similarity heatmap and structural projection scatter
- Export pages supporting PDF, JSON, and FASTA generation
- Shared API client layer in `frontend/src/lib/api.ts`

### Frontend / backend integration

The frontend uses `NEXT_PUBLIC_API_URL` to target the backend API. It supports both modern endpoints such as `/api/analyze` and legacy compatibility paths like `/api/analysis/sequence`.

---

## Backend Details

The backend layer is implemented in Python using FastAPI and includes:

- **FastAPI** for HTTP routing, request validation, and error handling.
- **Pydantic** models for strict schema validation of request/response payloads.
- **FastAPI lifespan hooks** for initializing and closing the async SQLite database.
- **CORS middleware** allowing secure front-end access.

### Key backend services

- **Structure prediction**: Pure Python Nussinov-Jacobson dynamic programming with thermodynamic scoring and dot-bracket output.
- **Stability scoring**: Computes a Therapeutic Stability Score using GC content, folding energy, paired fraction, and ensemble confidence.
- **Similarity analysis**: Builds a molecular similarity landscape from both sequence and structure data.
- **Clustering**: Detects structural families and projects candidates into a conformational relationship space.
- **Interpretation engine**: Creates biotech-grade narrative descriptions from numeric analysis results.
- **Visualization support**: SVG generation and plot payload assembly for dashboard rendering.
- **Export generation**: PDF dossiers, JSON archives, and FASTA export for publication-ready artifacts.

### Backend endpoints

Major backend routes include:

- `GET /api/health`
- `POST /api/analyze`
- `POST /api/analysis/sequence`
- `POST /api/analysis/batch`
- `GET /api/analysis/history`
- `POST /api/analysis/similarity`
- `POST /api/analysis/families`
- `POST /api/structure/predict`
- `GET /api/structure/{id}/svg`
- `GET /api/candidates`
- `POST /api/export/pdf`
- `POST /api/export/json`
- `POST /api/export/fasta`

---

## Methodology

This project applies domain-aware sequence analysis and RNA biophysics principles:

- **Sequence validation**: Enforces minimum and maximum sequence lengths and valid nucleotide letters.
- **Secondary structure folding**: Uses the Nussinov algorithm to maximize base pair counts and derive dot-bracket notation.
- **Thermodynamic scoring**: Estimates free energy using simplified nearest-neighbor parameters and penalizes hairpin loops.
- **Confidence estimation**: Samples suboptimal structures to quantify structural confidence and ensemble diversity.
- **Similarity modeling**: Compares sequences and predicted structures to create a distance-based similarity matrix.
- **Candidate ranking**: Sorts results by Therapeutic Stability Score and structural quality.
- **Structural families**: Groups sequences into fold-based families and projects them into a relationship space for exploratory analysis.

---

## Notebook / Colab Support

The file `project_1 (5).ipynb` contains the remaining Colab/back-end experiments, prototype analysis workflows, and supporting notebook-based research. It is intended as a companion resource for developing, validating, and extending the RNA analysis pipeline.

---

## Running the project

### Backend

```powershell
cd backend
pip install -r requirements.txt
python run.py
```

### Frontend

```powershell
cd frontend
npm install
npm run dev
```

### Default URLs

- Backend API: `http://localhost:8000`
- Frontend app: `http://localhost:3000`
- API docs: `http://localhost:8000/api/docs`

---

## Notes for developers

- `backend/app/config.py` defines host/port, API prefix, and runtime settings.
- `backend/app/main.py` registers routes and handles CORS / error middleware.
- `frontend/src/lib/api.ts` is the centralized REST client for backend integration.
- `backend/app/services/structure_prediction.py` is the main RNA folding engine.
- Persistence is handled through async SQLite with models in `backend/app/database/models.py`.

---

## Recommended next steps

1. Review `project_1 (5).ipynb` for Colab-based backend experiments.
2. Extend the RNA interpretation rules in `backend/app/interpretations/genomic_interpreter.py`.
3. Add additional export formats or candidate comparison metrics.
4. Validate the model with real RNA datasets or synthetic test cases.

---

## License

Include your preferred license here or add one when the project is ready for sharing.
