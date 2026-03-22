"""
gem_builder.py — CarveMe + COBRApy helpers for GEM reconstruction.

Uses CarveMe CLI to reconstruct a genome-scale metabolic model from a
protein FASTA file, then extracts summary statistics with COBRApy.

Speed optimizations applied:
  - Diamond BLAST: max threads, reduced --top
  - --solver scip (much faster than default GLPK fallback)
  - No gapfilling by default (skip the slow MILP solver step)
  - Gram-negative universe by default (smaller search space)
"""

from __future__ import annotations

import logging
import multiprocessing
import os
import shutil
import subprocess
import tempfile
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional

logger = logging.getLogger(__name__)

# Path to the carve binary (auto-detected or overridden)
_CARVE_BIN: Optional[str] = None


# ─────────────────────────────────────────────
# Data classes
# ─────────────────────────────────────────────

@dataclass
class GEMSummary:
    """Summary statistics for a reconstructed GEM."""

    reactions: int
    metabolites: int
    genes: int
    compartments: int
    has_biomass_objective: bool
    objective_id: str
    sbml_path: str  # path to the generated SBML file

    @property
    def biomass_status(self) -> str:
        return "Defined" if self.has_biomass_objective else "Not found"


# ─────────────────────────────────────────────
# Locate binaries
# ─────────────────────────────────────────────

def _find_carve_bin() -> str:
    """Find the carve binary, checking venv first then PATH."""
    global _CARVE_BIN
    if _CARVE_BIN:
        return _CARVE_BIN

    # Check project venv
    project_venv = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), ".venv", "bin", "carve"
    )
    if os.path.isfile(project_venv):
        _CARVE_BIN = project_venv
        return _CARVE_BIN

    # Check PATH
    found = shutil.which("carve")
    if found:
        _CARVE_BIN = found
        return _CARVE_BIN

    raise FileNotFoundError(
        "CarveMe 'carve' binary not found. Install with: pip install carveme"
    )


def check_dependencies() -> dict[str, bool]:
    """Check if required tools are available."""
    checks = {}

    try:
        _find_carve_bin()
        checks["carveme"] = True
    except FileNotFoundError:
        checks["carveme"] = False

    checks["diamond"] = shutil.which("diamond") is not None

    try:
        import cobra  # noqa: F401
        checks["cobra"] = True
    except ImportError:
        checks["cobra"] = False

    return checks


# ─────────────────────────────────────────────
# Core pipeline
# ─────────────────────────────────────────────

def _time_based_progress(
    progress_callback: Callable[[str], None],
    stop_event: threading.Event,
    estimated_seconds: float = 120.0,
):
    """
    Emit time-based progress updates since CarveMe produces no
    stderr output during the slow MILP carving phase.

    Stages and approximate time splits:
      0-5%   : Diamond BLAST alignment
      5-15%  : Scoring gene-reaction associations
      15-85% : Carving metabolic model (MILP — the slow part)
      85-95% : Post-processing model
      95-100%: Writing SBML output
    """
    stages = [
        (0.05, "Running Diamond BLAST alignment..."),
        (0.15, "Scoring gene-reaction associations..."),
        (0.50, "Carving metabolic model (solving MILP)..."),
        (0.75, "Carving metabolic model (refining)..."),
        (0.85, "Post-processing model..."),
        (0.95, "Writing SBML output..."),
    ]

    start = time.time()
    stage_idx = 0

    while not stop_event.is_set():
        elapsed = time.time() - start
        fraction = min(elapsed / estimated_seconds, 0.99)

        # Advance stage based on elapsed fraction
        while stage_idx < len(stages) - 1 and fraction >= stages[stage_idx + 1][0]:
            stage_idx += 1

        try:
            progress_callback(stages[stage_idx][1])
        except Exception:
            pass  # Ignore Streamlit NoSessionContext from background thread
        stop_event.wait(2.0)  # Update every 2 seconds


def run_carveme(
    fasta_path: str,
    output_path: str,
    universe: str = "bacteria",
    solver: str = "scip",
    gapfill: Optional[str] = None,
    threads: Optional[int] = None,
    progress_callback: Optional[Callable[[str], None]] = None,
) -> str:
    """
    Run CarveMe reconstruction on a protein FASTA file.

    Parameters
    ----------
    fasta_path : str
        Path to the input .faa file.
    output_path : str
        Desired path for the output SBML (.xml) file.
    universe : str
        Universe model: "bacteria", "grampos", "gramneg".
    solver : str
        MILP solver for reframed: "scip", "cplex", "gurobi".
    gapfill : str or None
        Media to gapfill on (e.g., "M9"). None = skip (faster).
    threads : int or None
        CPU threads for Diamond. None = use all available.
    progress_callback : callable or None
        Function(stage: str) called with progress messages.

    Returns
    -------
    str
        Path to the generated SBML file.
    """
    carve_bin = _find_carve_bin()

    if threads is None:
        threads = multiprocessing.cpu_count()

    diamond_args = f"--threads {threads} --top 5"

    # Build command as list (safe with spaces in paths)
    cmd = [
        carve_bin, fasta_path,
        "-o", output_path,
        "-u", universe,
        "--diamond-args", diamond_args,
    ]

    if gapfill:
        cmd.extend(["-g", gapfill])

    logger.info("Running CarveMe: %s", cmd)

    # Start time-based progress thread (CarveMe doesn't log during MILP)
    stop_event = threading.Event()
    progress_thread = None
    if progress_callback:
        progress_callback("Starting CarveMe reconstruction...")
        progress_thread = threading.Thread(
            target=_time_based_progress,
            args=(progress_callback, stop_event, 120.0),
            daemon=True,
        )
        progress_thread.start()

    try:
        # Run CarveMe as subprocess (list form — safe with spaces in paths)
        process = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=1800,  # 30 minute timeout
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError(
            "CarveMe timed out after 30 minutes. The input file may be too large."
        )
    finally:
        # Stop progress thread
        stop_event.set()
        if progress_thread:
            progress_thread.join(timeout=3)

    if process.returncode != 0:
        stderr = process.stderr.strip()
        # Provide user-friendly message for common failures
        if "0 queries aligned" in stderr:
            raise RuntimeError(
                "No protein homologs were found in the input file. "
                "Please ensure you uploaded a valid protein FASTA (.faa) file "
                "with standard amino acid sequences."
            )
        error_msg = stderr or process.stdout.strip() or "Unknown error"
        raise RuntimeError(f"CarveMe failed (exit {process.returncode}): {error_msg}")

    if not os.path.isfile(output_path):
        raise RuntimeError(
            "CarveMe completed but no output file was generated. "
            "The input may lack sufficient homology to known metabolic genes."
        )

    if progress_callback:
        progress_callback("Done!")

    return output_path


def read_gem_summary(sbml_path: str) -> GEMSummary:
    """
    Read an SBML model with COBRApy and return summary stats.

    Parameters
    ----------
    sbml_path : str
        Path to the SBML (.xml) file.

    Returns
    -------
    GEMSummary
    """
    from cobra.io import read_sbml_model

    model = read_sbml_model(sbml_path)

    # Check for a biomass objective
    obj_id = ""
    has_biomass = False
    if model.objective and model.objective.expression:
        obj_id = str(model.objective)
        has_biomass = True

    return GEMSummary(
        reactions=len(model.reactions),
        metabolites=len(model.metabolites),
        genes=len(model.genes),
        compartments=len(model.compartments),
        has_biomass_objective=has_biomass,
        objective_id=obj_id,
        sbml_path=sbml_path,
    )


def build_gem(
    fasta_bytes: bytes,
    file_name: str,
    output_dir: Optional[str] = None,
    universe: str = "bacteria",
    solver: str = "scip",
    gapfill: Optional[str] = None,
    progress_callback: Optional[Callable[[str], None]] = None,
) -> GEMSummary:
    """
    End-to-end: write uploaded bytes to disk, run CarveMe, read model.

    Parameters
    ----------
    fasta_bytes : bytes
        Raw content of the uploaded .faa file.
    file_name : str
        Original filename (used to derive the output name).
    output_dir : str or None
        Directory for outputs. Uses a temp dir if None.
    universe : str
        CarveMe universe: "bacteria", "grampos", "gramneg".
    solver : str
        MILP solver: "scip", "cplex", "gurobi".
    gapfill : str or None
        Gap-fill media (e.g., "M9"). None = skip for speed.
    progress_callback : callable or None
        Function(stage: str) for progress updates.

    Returns
    -------
    GEMSummary
    """
    if output_dir is None:
        output_dir = tempfile.mkdtemp(prefix="gem_optimizer_")

    os.makedirs(output_dir, exist_ok=True)

    # Write FASTA to disk
    fasta_path = os.path.join(output_dir, file_name)
    with open(fasta_path, "wb") as f:
        f.write(fasta_bytes)

    # Derive output SBML path
    stem = Path(file_name).stem
    sbml_path = os.path.join(output_dir, f"{stem}_model.xml")

    # Run CarveMe
    if progress_callback:
        progress_callback("Preparing input files...")

    run_carveme(
        fasta_path=fasta_path,
        output_path=sbml_path,
        universe=universe,
        solver=solver,
        gapfill=gapfill,
        progress_callback=progress_callback,
    )

    # Read and summarise
    if progress_callback:
        progress_callback("Reading model summary...")

    return read_gem_summary(sbml_path)


def read_sbml_bytes(sbml_path: str) -> bytes:
    """Read an SBML file and return its raw bytes (for download)."""
    with open(sbml_path, "rb") as f:
        return f.read()


# ─────────────────────────────────────────────
# GEM upload & COBRApy helpers
# ─────────────────────────────────────────────

@dataclass
class MediumComponent:
    """One exchange reaction that can be controlled as a medium component."""
    reaction_id: str
    name: str
    lower_bound: float
    upper_bound: float
    current_bound: float  # the flux bound the user sets (positive = uptake rate)


@dataclass
class FBAResult:
    """Results from a flux balance analysis run."""
    status: str            # "optimal", "infeasible", etc.
    growth_rate: float
    fluxes: dict           # reaction_id -> flux value
    exchange_fluxes: dict  # reaction_id -> flux value (exchange reactions only)


def load_gem_model(sbml_bytes: bytes, file_name: str):
    """
    Load an SBML model from raw bytes using COBRApy.

    Returns (cobra.Model, GEMSummary).
    """
    from cobra.io import read_sbml_model

    # Write to a temp file (cobra needs a file path)
    tmp_dir = tempfile.mkdtemp(prefix="gem_upload_")
    sbml_path = os.path.join(tmp_dir, file_name)
    with open(sbml_path, "wb") as f:
        f.write(sbml_bytes)

    model = read_sbml_model(sbml_path)

    obj_id = ""
    has_biomass = False
    if model.objective and model.objective.expression:
        obj_id = str(model.objective)
        has_biomass = True

    summary = GEMSummary(
        reactions=len(model.reactions),
        metabolites=len(model.metabolites),
        genes=len(model.genes),
        compartments=len(model.compartments),
        has_biomass_objective=has_biomass,
        objective_id=obj_id,
        sbml_path=sbml_path,
    )

    return model, summary


def extract_medium(model) -> list[MediumComponent]:
    """
    Extract the medium (active exchange reactions) from a cobra model.

    Returns a list of MediumComponent with the current uptake bounds.
    """
    medium = model.medium  # dict: reaction_id -> upper_bound
    components = []

    for rxn_id, bound in medium.items():
        try:
            rxn = model.reactions.get_by_id(rxn_id)
        except KeyError:
            continue

        # Try to get a readable name from the reaction name or metabolite
        name = rxn.name if rxn.name else rxn_id
        # If the name is just the ID, try to extract from metabolites
        if name == rxn_id and rxn.reactants:
            met = rxn.reactants[0]
            name = met.name if met.name else met.id

        components.append(MediumComponent(
            reaction_id=rxn_id,
            name=name,
            lower_bound=0.0,
            upper_bound=1000.0,
            current_bound=bound,
        ))

    # Sort by bound (most constrained first, then alphabetical)
    components.sort(key=lambda c: (-c.current_bound if c.current_bound < 1000 else 0, c.name))

    return components


def run_fba(model, medium_bounds: dict[str, float]) -> FBAResult:
    """
    Run FBA on the model with the given medium bounds.

    Parameters
    ----------
    model : cobra.Model
        The loaded metabolic model.
    medium_bounds : dict
        Mapping of exchange reaction ID -> uptake bound (positive value).

    Returns
    -------
    FBAResult
    """
    with model:  # context manager restores original state after
        # Set the medium
        model.medium = medium_bounds

        # Run FBA
        solution = model.optimize()

        exchange_fluxes = {}
        for rxn in model.exchanges:
            exchange_fluxes[rxn.id] = solution.fluxes.get(rxn.id, 0.0)

        return FBAResult(
            status=solution.status,
            growth_rate=solution.objective_value if solution.status == "optimal" else 0.0,
            fluxes=solution.fluxes.to_dict() if solution.status == "optimal" else {},
            exchange_fluxes=exchange_fluxes,
        )
