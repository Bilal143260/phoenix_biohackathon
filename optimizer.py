"""
optimizer.py — Bayesian Optimization of media bounds using Optuna + COBRApy FBA.

Bilevel optimization:
  - Inner: FBA finds optimal fluxes for a given media (COBRApy)
  - Outer: Optuna searches for best media bounds under a budget constraint

Two budget modes:
  1. Flux budget: sum of all uptake bounds <= budget
  2. Cost budget: sum(bound_i * cost_i) <= budget
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Callable, Optional

import optuna

# Silence Optuna's verbose logging
optuna.logging.set_verbosity(optuna.logging.WARNING)

logger = logging.getLogger(__name__)


@dataclass
class OptimizationConfig:
    """Configuration for the outer (Bayesian) optimization loop."""
    budget: float = 100.0
    budget_mode: str = "flux"  # "flux" or "cost"
    costs: dict[str, float] = field(default_factory=dict)  # reaction_id -> $/unit
    n_trials: int = 50
    max_bound: float = 1000.0  # upper limit for any single exchange


@dataclass
class TrialRecord:
    """One Optuna trial's inputs and outputs."""
    number: int
    bounds: dict[str, float]
    growth_rate: float
    total_cost: float
    feasible: bool


@dataclass
class OptimizationResult:
    """Full result of the Bayesian optimization run."""
    best_growth_rate: float
    best_bounds: dict[str, float]
    best_cost: float
    trials: list[TrialRecord]
    config: OptimizationConfig


def optimize_media(
    model,
    medium_components: list,  # list of MediumComponent
    config: OptimizationConfig,
    progress_callback: Optional[Callable[[int, int, float], None]] = None,
) -> OptimizationResult:
    """
    Run Bayesian Optimization over media bounds to maximize growth rate.

    Parameters
    ----------
    model : cobra.Model
        The loaded metabolic model.
    medium_components : list[MediumComponent]
        Exchange reactions to optimize (from extract_medium).
    config : OptimizationConfig
        Budget, cost map, number of trials, etc.
    progress_callback : callable or None
        Called with (trial_number, total_trials, best_so_far) after each trial.

    Returns
    -------
    OptimizationResult
    """
    from gem_builder import run_fba

    rxn_ids = [c.reaction_id for c in medium_components]
    default_bounds = {c.reaction_id: c.current_bound for c in medium_components}

    # Cost per unit for each reaction (default 1.0 for flux mode)
    cost_map = {}
    for rid in rxn_ids:
        if config.budget_mode == "cost":
            cost_map[rid] = config.costs.get(rid, 1.0)
        else:
            cost_map[rid] = 1.0

    trials: list[TrialRecord] = []

    def objective(trial: optuna.Trial) -> float:
        # Strategy: sample a fraction of budget for each component.
        # Each component gets a "share" (0 to 1), then we allocate budget
        # proportionally. This guarantees every trial is feasible.
        shares = {}
        for comp in medium_components:
            rid = comp.reaction_id
            shares[rid] = trial.suggest_float(f"{rid}_share", 0.0, 1.0)

        total_share = sum(s / cost_map[rid] for rid, s in shares.items())
        if total_share < 1e-9:
            total_share = 1.0

        # Allocate budget proportionally, capped at max_bound
        bounds = {}
        total_cost = 0.0
        for comp in medium_components:
            rid = comp.reaction_id
            upper = min(comp.upper_bound, config.max_bound)
            raw_alloc = (shares[rid] / cost_map[rid]) / total_share * config.budget / max(cost_map[rid], 1e-9)
            bounds[rid] = min(raw_alloc, upper)
            total_cost += bounds[rid] * cost_map[rid]

        # Run FBA with these bounds
        try:
            result = run_fba(model, bounds)
            growth = result.growth_rate if result.status == "optimal" else 0.0
        except Exception:
            growth = 0.0

        record = TrialRecord(
            number=trial.number,
            bounds=bounds,
            growth_rate=growth,
            total_cost=total_cost,
            feasible=True,
        )
        trials.append(record)

        if progress_callback:
            best = max(t.growth_rate for t in trials)
            progress_callback(trial.number + 1, config.n_trials, best)

        return growth

    # Create and run the study
    study = optuna.create_study(
        direction="maximize",
        sampler=optuna.samplers.TPESampler(seed=42),
    )

    # Seed with an equal-shares starting point
    equal_share = {f"{rid}_share": 1.0 / len(rxn_ids) for rid in rxn_ids}
    study.enqueue_trial(equal_share)

    study.optimize(objective, n_trials=config.n_trials, show_progress_bar=False)

    # Extract best feasible trial
    feasible_trials = [t for t in trials if t.feasible and t.growth_rate > 0]
    if feasible_trials:
        best_trial = max(feasible_trials, key=lambda t: t.growth_rate)
    else:
        default_total = sum(default_bounds[rid] * cost_map[rid] for rid in rxn_ids)
        best_trial = TrialRecord(
            number=-1,
            bounds=default_bounds,
            growth_rate=0.0,
            total_cost=default_total,
            feasible=False,
        )

    return OptimizationResult(
        best_growth_rate=best_trial.growth_rate,
        best_bounds=best_trial.bounds,
        best_cost=best_trial.total_cost,
        trials=trials,
        config=config,
    )
