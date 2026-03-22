# GEM Media Optimizer

**Genome-to-media optimization for microbial strains**

A web application that helps bioprocess engineers and non-specialists optimize microbial growth media using genome-scale metabolic models (GEMs). Upload a pre-built SBML model, explore and tune media bounds in real time, and let an AI-powered optimizer find the best composition under a budget constraint.

## Quick Start

### Prerequisites
- Python 3.9+

### Install dependencies

```bash
pip install -r requirements.txt
```

### Run the app

```bash
streamlit run app.py
```

The app will open at `http://localhost:8501`.

## Workflow

The app has two main tabs:

### Manual Simulator
1. **Upload** — Provide a protein FASTA file or upload a pre-built GEM model (SBML/XML)
2. **GEM Generation** — Build a draft genome-scale metabolic model, or skip if uploading SBML
3. **Media Design** — Explore extracted medium components with interactive bound sliders
4. **Optimize** — Run Flux Balance Analysis (FBA) with your configured bounds
5. **Results** — View predicted growth rate, key drivers, recommended ranges, and download outputs

### AI Media Optimizer
1. **Upload a GEM model** in the Manual Simulator tab first
2. **Configure budget** — Set a total flux limit or cost-weighted budget so the optimizer can't cheat
3. **Assign costs** (optional) — In cost mode, set per-metabolite prices (cheap glucose vs expensive amino acids)
4. **Run optimization** — Bayesian Optimization (Optuna TPE sampler) searches for the best media allocation
5. **Review results** — Convergence plot, optimal bounds, side-by-side comparison, JSON/CSV export

## How the AI Optimizer Works

This is a **bilevel optimization** problem:

- **Inner optimization**: COBRApy's FBA finds optimal fluxes for a given set of media bounds
- **Outer optimization**: Optuna's Tree-structured Parzen Estimator (TPE) searches for the best media bounds

The key insight is the **budget constraint** — without it, the optimizer would simply set every nutrient to infinity. Two budget modes are supported:

| Mode | Description |
|------|-------------|
| **Total Flux Limit** | Sum of all uptake bounds ≤ budget |
| **Cost-Weighted** | Each metabolite has a $/unit cost; total cost ≤ budget |

The optimizer uses **proportional allocation**: it samples "shares" for each component and distributes the budget proportionally, guaranteeing 100% trial feasibility.

## Features

- Upload pre-built GEM models (SBML/XML) — skip lengthy reconstruction
- COBRApy-powered FBA with real metabolic models
- Interactive media bound sliders extracted from the model
- AI-powered Bayesian Optimization with budget constraints
- Real-time optimization progress tracking
- Convergence visualization and default vs. optimized comparison charts
- Card-based media builder (no spreadsheet tables)
- Sensitivity-ranked media driver chart
- JSON and CSV export for both manual and AI results
- Premium scientific SaaS design

## Tech Stack

| Component | Technology |
|-----------|-----------|
| Frontend | Streamlit |
| Charts | Plotly |
| Data | pandas |
| Metabolic Modeling | COBRApy |
| AI Optimization | Optuna (TPE Bayesian Optimization) |
| Styling | Custom CSS |

## File Structure

```
app.py           — Main Streamlit application (Manual Simulator + AI Optimizer tabs)
gem_builder.py   — COBRApy integration (model loading, medium extraction, FBA)
optimizer.py     — Bayesian optimization engine (Optuna + FBA bilevel optimization)
requirements.txt — Python dependencies
CLAUDE.md        — AI assistant instructions
README.md        — This file
```

## License

Hackathon project — all rights reserved.
