# GEM Media Optimizer

**Genome-to-media optimization for microbial strains**

A web application that helps bioprocess engineers and non-specialists upload protein sequences, generate genome-scale metabolic models, and discover optimal media compositions for maximizing microbial growth.

## Quick Start

### Prerequisites
- Python 3.9+

### Install dependencies

```bash
pip install streamlit plotly pandas
```

### Run the app

```bash
streamlit run app.py
```

The app will open at `http://localhost:8501`.

## Workflow

1. **Upload** — Provide a protein FASTA file (`.faa`) or load the example dataset
2. **Generate GEM** — Build a draft genome-scale metabolic model from the input
3. **Build Media** — Select and configure media components using the card-based builder
4. **Optimize** — Run flux balance analysis to find the best media composition for growth
5. **Results** — View predicted growth, key drivers, recommended ranges, and download outputs

## Features

- Drag-and-drop file upload with example dataset
- Card-based media builder (no spreadsheet tables)
- Recommended starting set for quick setup
- Sensitivity-ranked media driver chart
- Recommended concentration ranges
- Optimized recipe comparison (original vs. optimized)
- Performance proxies (growth rate, yield proxy, product rate proxy)
- JSON and CSV export
- Premium scientific SaaS design

## Tech Stack

| Component | Technology |
|-----------|-----------|
| Frontend | Streamlit |
| Charts | Plotly |
| Data | pandas |
| Styling | Custom CSS |

## Notes

This is a hackathon MVP. All biological outputs are simulated with realistic mock data. No backend GEM solver is connected yet — the app demonstrates the complete user experience and is ready for integration with tools like COBRApy.

## License

Hackathon project — all rights reserved.
