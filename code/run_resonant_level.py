"""Driver: resonant-level test of the dynamical mechanism.

Solves the quench of `resonant_level.py` in two regimes that differ only in the
bath bandwidth, and writes

    code/out/<timestamp>_resonant_level/
        metadata.json      -- parameters and all extracted numbers
        results.npz        -- raw arrays behind the figure
        rlm_beta_eff.pdf   -- the figure (also copied into the report tree)

What is and is not being tested
------------------------------
Note the scope: gamma0 is held FIXED at 1.0 and only the bath bandwidth is
varied.  This run therefore probes the presence or absence of bath memory, not
the coupling threshold that the reference discusses in terms of V.  No claim
about a threshold in gamma0 can be based on it.

  (i)   Wide band (Markovian): beta_eff(t) approaches beta_bath without
        transient structure.
        -- CONFIRMED.  beta_eff settles to 0.24997 and every feature, both in
        time and across fitting windows, stays within ~1e-4.  (The strict
        monotonicity flag below is False only at that amplitude.)

  (ii)  Finite bandwidth: beta_eff(t) develops large non-monotonic transients,
        including excursions to negative values.
        -- CONFIRMED.  Since this model is quadratic, *effective-temperature
        oscillation and negative effective inverse temperature* need neither
        interaction, disorder, scrambling nor chaos.  Note that beta_eff < 0
        here means the FDT fit returns a negative slope; it does not mean the
        state is a thermodynamic negative-temperature state.

  (iii) Trajectories started from different beta_init cross (an MPC).
        -- NOT ESTABLISHED.  At *every* reported time the spread of beta_eff
        over the fitting windows exceeds the separation between the two
        trajectories (`separation_exceeds_spread_anywhere` is False; spread
        median 0.043, max 1.318, versus separation max 0.164).  The crossing is
        buried in the systematic uncertainty of the thermometer, so this run
        can neither confirm nor exclude an MPC -- and in particular it does NOT
        show that chaos is unnecessary for an MPC, only for (ii).

Meanwhile the spectrally weighted distance D_bath(t) between the nonequilibrium
distribution and the bath Fermi function decreases monotonically in *both*
regimes.  D_bath is itself a constructed distance measure -- it depends on the
frequency window, the normalisation and the spectral weight -- not a state
function, so the contrast with beta_eff is suggestive rather than decisive: the
distribution-level relaxation and the thermometer reading behave differently,
which is the same structure the reference reports for the SYK model, here with
the interaction switched off entirely.

Usage
-----
    python3 run_resonant_level.py
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

from resonant_level import (
    ResonantLevel,
    distance_to_bath,
    effective_beta,
    selftest_equilibrium,
    solve_two_time,
)

EPS_D = 1.0
GAMMA0 = 1.0
BETA_BATH = 0.25
BETA_INITS = (0.02, 0.20)
BANDWIDTHS = {"wide": 20.0, "narrow": 1.0}

T_GRID_MAX = 40.0
N_TIME = 1201
REL_MAX = 16.0
OMEGA_WINDOWS = (2.0, 2.5, 3.0, 3.5, 4.0)  # spread over these gauges the systematics
OMEGA_REFERENCE = 3.0
T_REPORT = np.arange(2.0, 20.01, 0.25)

HERE = Path(__file__).resolve().parent
OUT_ROOT = HERE / "out"
REPORT_FIGS = HERE.parent / "notes_tex" / "report_quantum_mpemba" / "figs"


def _style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 8,
            "axes.labelsize": 8,
            "legend.fontsize": 6.5,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 0.6,
            "lines.linewidth": 1.1,
            "xtick.direction": "in",
            "ytick.direction": "in",
        }
    )


def analyse(model: ResonantLevel, times: np.ndarray, omegas: np.ndarray) -> dict:
    """Run one parameter set and extract beta_eff, its window spread and D_bath."""
    solution = solve_two_time(model, times)
    indices = [int(np.argmin(np.abs(times - T))) for T in T_REPORT]

    reference, low, high, distance = [], [], [], []
    for i in indices:
        values = [
            effective_beta(solution, i, omega_window=w, rel_max=REL_MAX)
            for w in OMEGA_WINDOWS
        ]
        reference.append(
            effective_beta(solution, i, omega_window=OMEGA_REFERENCE, rel_max=REL_MAX)
        )
        low.append(min(values))
        high.append(max(values))
        distance.append(distance_to_bath(solution, i, model, omegas, rel_max=REL_MAX))

    return {
        "beta_eff": np.array(reference),
        "beta_low": np.array(low),
        "beta_high": np.array(high),
        "distance": np.array(distance),
        "completeness_error": solution.completeness_error,
    }


def main() -> None:
    _style()

    check = selftest_equilibrium()
    if not check["passed"]:
        raise RuntimeError(f"equilibrium self-test failed: {check}")

    times = np.linspace(0.0, T_GRID_MAX, N_TIME)
    omegas = np.linspace(-8.0, 8.0, 401)

    results: dict[str, dict[float, dict]] = {}
    for tag, bandwidth in BANDWIDTHS.items():
        results[tag] = {}
        for beta_init in BETA_INITS:
            model = ResonantLevel(
                eps_d=EPS_D,
                gamma0=GAMMA0,
                bandwidth=bandwidth,
                beta_bath=BETA_BATH,
                beta_init=beta_init,
            )
            results[tag][beta_init] = analyse(model, times, omegas)

    # --- figure -----------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(6.1, 2.0))
    colours = {BETA_INITS[0]: "C1", BETA_INITS[1]: "C2"}

    for ax, tag, title in (
        (axes[0], "wide", rf"(a) $\Lambda={BANDWIDTHS['wide']:.0f}$ (Markovian)"),
        (axes[1], "narrow", rf"(b) $\Lambda={BANDWIDTHS['narrow']:.0f}$"),
    ):
        for beta_init in BETA_INITS:
            data = results[tag][beta_init]
            ax.fill_between(
                T_REPORT,
                data["beta_low"],
                data["beta_high"],
                color=colours[beta_init],
                alpha=0.25,
                linewidth=0,
            )
            ax.plot(
                T_REPORT,
                data["beta_eff"],
                color=colours[beta_init],
                label=rf"$\beta_0={beta_init:.2f}$",
            )
        ax.axhline(BETA_BATH, color="0.5", lw=0.6, ls="--")
        ax.set_xlabel(r"time $t$")
        ax.set_ylabel(r"$\beta_{\rm eff}(t)$")
        ax.set_title(title, fontsize=8)
        ax.set_xlim(T_REPORT[0], T_REPORT[-1])

    axes[0].set_ylim(0.0, 0.5)
    axes[1].set_ylim(-1.0, 1.5)
    axes[0].legend(frameon=False, loc="lower right")
    axes[0].text(
        0.96,
        0.62,
        r"$\beta_{\rm bath}$",
        transform=axes[0].transAxes,
        ha="right",
        color="0.4",
        fontsize=6.5,
    )

    ax = axes[2]
    for tag, style in (("wide", "-"), ("narrow", "--")):
        for beta_init in BETA_INITS:
            ax.semilogy(
                T_REPORT,
                np.clip(results[tag][beta_init]["distance"], 1e-8, None),
                style,
                color=colours[beta_init],
                label=rf"$\Lambda={BANDWIDTHS[tag]:.0f}$, $\beta_0={beta_init:.2f}$",
            )
    ax.set_xlabel(r"time $t$")
    ax.set_ylabel(r"$D_{\rm bath}(t)$")
    ax.set_title("(c) distribution-level distance", fontsize=8)
    ax.set_xlim(T_REPORT[0], T_REPORT[-1])
    ax.legend(frameon=False, loc="upper right")

    fig.tight_layout(pad=0.4)

    # --- output -----------------------------------------------------------
    stamp = datetime.now(timezone.utc).astimezone().strftime("%Y-%m-%d_%H%M%S")
    out_dir = OUT_ROOT / f"{stamp}_resonant_level"
    out_dir.mkdir(parents=True, exist_ok=True)

    figure_path = out_dir / "rlm_beta_eff.pdf"
    fig.savefig(figure_path)
    fig.savefig(out_dir / "rlm_beta_eff.png", dpi=200)
    plt.close(fig)

    REPORT_FIGS.mkdir(parents=True, exist_ok=True)
    shutil.copy(figure_path, REPORT_FIGS / "rlm_beta_eff.pdf")

    # Quantify the central claim: trajectory separation vs. window systematics.
    late = T_REPORT >= 3.0
    summary = {}
    for tag in BANDWIDTHS:
        a, b = (results[tag][bi] for bi in BETA_INITS)
        separation = np.abs(a["beta_eff"] - b["beta_eff"])[late]
        spread = np.maximum(
            a["beta_high"] - a["beta_low"], b["beta_high"] - b["beta_low"]
        )[late]
        summary[tag] = {
            "max_trajectory_separation": float(separation.max()),
            "median_window_spread": float(np.median(spread)),
            "max_window_spread": float(spread.max()),
            "separation_exceeds_spread_anywhere": bool(np.any(separation > spread)),
            "beta_eff_monotone": bool(
                np.all(np.diff(a["beta_eff"][late]) >= -1e-6)
                or np.all(np.diff(a["beta_eff"][late]) <= 1e-6)
            ),
            "distance_monotone": bool(
                np.all(np.diff(a["distance"][late]) <= 1e-9)
            ),
            "beta_eff_final": float(a["beta_eff"][-1]),
        }

    metadata = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "script": Path(__file__).name,
        "eps_d": EPS_D,
        "gamma0": GAMMA0,
        "beta_bath": BETA_BATH,
        "beta_inits": list(BETA_INITS),
        "bandwidths": BANDWIDTHS,
        "rel_max": REL_MAX,
        "omega_windows": list(OMEGA_WINDOWS),
        "omega_reference": OMEGA_REFERENCE,
        "equilibrium_selftest": check,
        "completeness_errors": {
            tag: results[tag][BETA_INITS[0]]["completeness_error"]
            for tag in BANDWIDTHS
        },
        "summary": summary,
    }
    (out_dir / "metadata.json").write_text(json.dumps(metadata, indent=2))

    np.savez_compressed(
        out_dir / "results.npz",
        times=T_REPORT,
        **{
            f"{tag}_{str(bi).replace('.', 'p')}_{key}": results[tag][bi][key]
            for tag in BANDWIDTHS
            for bi in BETA_INITS
            for key in ("beta_eff", "beta_low", "beta_high", "distance")
        },
    )

    print(f"equilibrium self-test : {check}")
    for tag, s in summary.items():
        print(f"[{tag}]")
        for k, v in s.items():
            print(f"    {k:<38}: {v}")
    print(f"written to            : {out_dir}")


if __name__ == "__main__":
    main()
