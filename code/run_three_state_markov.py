"""Driver: produce the three-state Markov Mpemba figure and numerical record.

Runs the model defined in `mpemba_three_state.py` at the parameter point found
by the grid search documented below, and writes

    code/out/<timestamp>_three_state_markov/
        metadata.json   -- all parameters and derived quantities
        results.npz     -- raw arrays behind the figure
        markov_a2.pdf   -- the figure (also copied into the report tree)

Chosen parameter point
----------------------
    E = (0, 1, 5),  Gamma = (0.3, 2.0, 0.001),  beta_bath = 0.5

This point was selected from a grid scan over (E_3, Gamma_12, Gamma_13,
Gamma_23, beta_bath) by requiring, in this order:

  1. a non-degenerate slowest mode, lambda_2 > lambda_3 (Lu-Raz condition);
  2. a *second* root beta* of a_2(beta_0) inside (0, beta_bath), i.e. a strong
     Mpemba effect, placed well away from both endpoints;
  3. a physically sensible relaxation rate, 0.05 <= |lambda_2| <= 5, which
     rejects the near-metastable solutions (|lambda_2| ~ 1e-3) that the raw
     scan otherwise favours because their gap ratio diverges;
  4. a clean but not pathological spectral gap, 3 <= lambda_3 / lambda_2 <= 25.

The resulting point has lambda_3 / lambda_2 = 12.6 and beta* = 0.1246.

Demonstration
-------------
Two initial inverse temperatures straddling beta* are propagated:

    beta_0 = 0.100  (hotter, |a_2| = 0.0110)
    beta_0 = 0.294  (colder, |a_2| = 0.0330)

The hotter state starts *farther* from stationarity in L1 distance (0.388 vs
0.172) yet overtakes the colder one at t = 0.305 and reaches ||p - pi||_1 <
1e-3 at t = 3.93 instead of t = 5.76.  A third trajectory started exactly at
beta* has a_2 = 0 and decays at the rate lambda_3, reaching the same threshold
at t = 0.78 -- the strong Mpemba effect.

The distance plotted is the unnormalised L1 norm ||p - pi||_1 = sum_i |p_i -
pi_i|, i.e. twice the total variation distance.

Usage
-----
    python3 run_three_state_markov.py
"""

from __future__ import annotations

import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mpemba_three_state import (
    ThreeStateModel,
    a2_of_beta,
    build_rate_matrix,
    evolve,
    gibbs,
    l1_distance,
    overlap_coefficients,
    spectral_decomposition,
    strong_mpemba_root,
)

# --- Parameter point (see module docstring for how it was selected) --------
MODEL = ThreeStateModel(energies=(0.0, 1.0, 5.0), gamma=(0.3, 2.0, 0.001), beta_bath=0.5)
BETA_HOT = 0.100
BETA_COLD = 0.294
T_MAX = 14.0
N_TIME = 2801
CONVERGENCE_THRESHOLD = 1e-3

HERE = Path(__file__).resolve().parent
OUT_ROOT = HERE / "out"
REPORT_FIGS = HERE.parent / "notes_tex" / "report_quantum_mpemba" / "figs"


def _style() -> None:
    """Match the figure style to the 11pt serif body text of the report."""
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 8,
            "axes.labelsize": 8,
            "axes.titlesize": 8,
            "legend.fontsize": 7,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 0.6,
            "lines.linewidth": 1.2,
            "xtick.direction": "in",
            "ytick.direction": "in",
        }
    )


def main() -> None:
    _style()

    stationary = gibbs(MODEL.energies, MODEL.beta_bath)
    rate = build_rate_matrix(MODEL)
    spectrum = spectral_decomposition(rate, stationary)
    if not spectrum.gap_condition:
        raise RuntimeError("lambda_2 > lambda_3 is violated at this parameter point")

    beta_star = strong_mpemba_root(MODEL, spectrum)
    if beta_star is None:
        raise RuntimeError("no strong Mpemba root at this parameter point")

    # --- panel (a): a_2 as a function of the initial inverse temperature ---
    betas = np.linspace(1e-3, MODEL.beta_bath, 1200)
    a2_curve = a2_of_beta(MODEL, spectrum, betas)

    # --- panel (b): relaxation from three initial temperatures -------------
    times = np.linspace(0.0, T_MAX, N_TIME)
    runs: dict[str, dict] = {}
    for label, beta0 in (
        ("hot", BETA_HOT),
        ("cold", BETA_COLD),
        ("star", float(beta_star)),
    ):
        p0 = gibbs(MODEL.energies, beta0)
        coeffs = overlap_coefficients(spectrum, p0)
        distance = l1_distance(evolve(spectrum, stationary, p0, times), stationary)
        below = np.nonzero(distance < CONVERGENCE_THRESHOLD)[0]
        runs[label] = {
            "beta_0": beta0,
            "a_2": float(coeffs[1]),
            "a_3": float(coeffs[2]),
            "initial_distance": float(distance[0]),
            "t_converged": float(times[below[0]]) if below.size else None,
            "distance": distance,
        }

    d_hot, d_cold = runs["hot"]["distance"], runs["cold"]["distance"]
    crossings = times[np.nonzero(np.diff(np.sign(d_hot - d_cold)))[0]]

    # --- figure ------------------------------------------------------------
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(5.4, 2.0))

    ax_a.axhline(0.0, color="0.7", lw=0.6)
    ax_a.plot(betas, a2_curve, color="C0")
    ax_a.plot([beta_star], [0.0], "o", ms=4, color="C3", zorder=5)
    ax_a.annotate(
        r"$\beta_0^{*}$",
        xy=(beta_star, 0.0),
        xytext=(beta_star + 0.02, -0.022),
        color="C3",
    )
    for beta0, color in ((BETA_HOT, "C1"), (BETA_COLD, "C2")):
        ax_a.axvline(beta0, color=color, ls=":", lw=0.9)
    ax_a.set_xlabel(r"initial inverse temperature $\beta_0$")
    ax_a.set_ylabel(r"overlap $a_2(\beta_0)$")
    ax_a.set_xlim(0.0, MODEL.beta_bath)
    ax_a.text(0.03, 0.92, "(a)", transform=ax_a.transAxes, va="top")
    ax_a.text(
        0.97,
        0.06,
        r"$\beta_{\rm bath}$",
        transform=ax_a.transAxes,
        ha="right",
        color="0.4",
    )

    ax_b.semilogy(
        times,
        d_hot,
        color="C1",
        label=rf"$\beta_0={BETA_HOT:.3f}$ (hotter)",
    )
    ax_b.semilogy(
        times,
        d_cold,
        color="C2",
        label=rf"$\beta_0={BETA_COLD:.3f}$ (colder)",
    )
    ax_b.semilogy(
        times,
        runs["star"]["distance"],
        color="C3",
        ls="--",
        label=rf"$\beta_0=\beta_0^{{*}}={beta_star:.3f}$",
    )
    if crossings.size:
        ax_b.axvline(crossings[0], color="0.7", lw=0.6, ls="-.")
    ax_b.set_xlabel(r"time $t$")
    ax_b.set_ylabel(r"$\|\varrho(t)-\varrho_{\rm ss}\|_1$")
    ax_b.set_xlim(0.0, 8.0)
    ax_b.set_ylim(1e-6, 1.0)
    ax_b.legend(frameon=False, loc="upper right")
    ax_b.text(0.03, 0.06, "(b)", transform=ax_b.transAxes)

    fig.tight_layout(pad=0.4)

    # --- output ------------------------------------------------------------
    stamp = datetime.now(timezone.utc).astimezone().strftime("%Y-%m-%d_%H%M%S")
    out_dir = OUT_ROOT / f"{stamp}_three_state_markov"
    out_dir.mkdir(parents=True, exist_ok=True)

    figure_path = out_dir / "markov_a2.pdf"
    fig.savefig(figure_path)
    fig.savefig(out_dir / "markov_a2.png", dpi=200)
    plt.close(fig)

    REPORT_FIGS.mkdir(parents=True, exist_ok=True)
    shutil.copy(figure_path, REPORT_FIGS / "markov_a2.pdf")

    metadata = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "script": Path(__file__).name,
        "model": MODEL.as_dict(),
        "stationary_distribution": stationary.tolist(),
        "eigenvalues": spectrum.eigenvalues.tolist(),
        "gap_ratio": float(spectrum.eigenvalues[2] / spectrum.eigenvalues[1]),
        "left_eigenvector_2": spectrum.left[1].tolist(),
        "beta_star": float(beta_star),
        "convergence_threshold": CONVERGENCE_THRESHOLD,
        "first_crossing_time": float(crossings[0]) if crossings.size else None,
        "runs": {
            k: {kk: vv for kk, vv in v.items() if kk != "distance"}
            for k, v in runs.items()
        },
    }
    (out_dir / "metadata.json").write_text(json.dumps(metadata, indent=2))

    np.savez_compressed(
        out_dir / "results.npz",
        betas=betas,
        a2_curve=a2_curve,
        times=times,
        distance_hot=d_hot,
        distance_cold=d_cold,
        distance_star=runs["star"]["distance"],
        eigenvalues=spectrum.eigenvalues,
        stationary=stationary,
    )

    print(f"eigenvalues        : {spectrum.eigenvalues}")
    print(f"gap ratio l3/l2    : {metadata['gap_ratio']:.2f}")
    print(f"beta*              : {beta_star:.5f}  (beta_bath = {MODEL.beta_bath})")
    for label, run in runs.items():
        print(
            f"{label:<5}: beta_0={run['beta_0']:.5f}  a_2={run['a_2']:+.5f}  "
            f"d(0)={run['initial_distance']:.4f}  t_conv={run['t_converged']}"
        )
    print(f"first crossing     : t = {crossings[0]:.3f}" if crossings.size else "no crossing")
    print(f"written to         : {out_dir}")


if __name__ == "__main__":
    main()
