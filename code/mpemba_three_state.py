"""Three-state Markov model: spectral mechanism of the Mpemba effect.

This module reproduces, in the smallest possible discrete setting, the
mechanism of Lu & Raz [PNAS 114, 5083 (2017); arXiv:1609.05271]: for a
Markovian relaxation the approach to the stationary state is controlled by
the overlap a_2 between the initial distribution and the slowest-decaying
eigenmode of the rate matrix.  When |a_2| is non-monotonic in the initial
inverse temperature, a hotter initial state can relax faster than a colder
one -- the Mpemba effect.

Physical setup
--------------
Three energy levels E_1 < E_2 < E_3 in contact with a bath at inverse
temperature beta_bath.  The rate matrix is parametrised as

    R_ij = Gamma_ij * exp(-beta_bath * (E_i - E_j) / 2)     (i != j)
    R_jj = -sum_{i != j} R_ij                               (conservation)

with a symmetric matrix of bare rates Gamma_ij = Gamma_ji >= 0.  This form
satisfies detailed balance with respect to the Gibbs distribution
pi_i ∝ exp(-beta_bath E_i) by construction, since

    R_ij pi_j ∝ Gamma_ij exp(-beta_bath (E_i + E_j) / 2)

is symmetric under i <-> j.

Because of detailed balance R is similar to the symmetric matrix
S = D^{-1} R D with D = diag(sqrt(pi)), so all eigenvalues are real and the
spectral decomposition can be obtained stably from `numpy.linalg.eigh`
applied to S rather than from a general non-symmetric eigensolver.

Conventions
-----------
Eigenvalues are ordered lambda_1 = 0 > lambda_2 >= lambda_3, i.e. mode 2 is
the slowest decaying one.  Note that the sufficient condition of Lu & Raz,
"lambda_2 > lambda_3", refers to this ordering and therefore means
|lambda_2| < |lambda_3|.

Right eigenvectors v_k and left eigenvectors l_k are biorthonormal,
<l_k, v_h> = delta_kh, so that

    p(t) = pi + sum_{k>=2} a_k exp(lambda_k t) v_k ,   a_k = <l_k, p(0)> .

a_2(beta_bath) = 0 holds identically, because l_2 is orthogonal to v_1 = pi.
A *strong* Mpemba effect corresponds to a second root of a_2(beta_0) at some
beta_0* != beta_bath; initial states on opposite sides of that root have
opposite signs of a_2, and |a_2| is non-monotonic between them.

Units
-----
Energies are measured in units of E_2 (so E_2 = 1 by convention below) and
rates in units of the overall bare rate scale; times are then dimensionless.

Author: (SYK_Thermodynamics project)
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Sequence

import numpy as np


# ---------------------------------------------------------------------------
# Model definition
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ThreeStateModel:
    """Parameters of the three-state detailed-balance Markov model.

    Attributes
    ----------
    energies:
        Level energies (E_1, E_2, E_3), assumed increasing.
    gamma:
        Symmetric bare rates (Gamma_12, Gamma_13, Gamma_23), all >= 0.
    beta_bath:
        Inverse temperature of the bath that defines the stationary state.
    """

    energies: tuple[float, float, float]
    gamma: tuple[float, float, float]
    beta_bath: float

    def as_dict(self) -> dict:
        """Return a JSON-serialisable description, for output metadata."""
        return asdict(self)


def gibbs(energies: Sequence[float], beta: float) -> np.ndarray:
    """Return the normalised Gibbs distribution at inverse temperature `beta`.

    The energies are shifted by their minimum before exponentiating, which is
    irrelevant after normalisation but avoids overflow at large `beta`.
    """
    e = np.asarray(energies, dtype=float)
    w = np.exp(-beta * (e - e.min()))
    return w / w.sum()


def build_rate_matrix(model: ThreeStateModel) -> np.ndarray:
    """Build the 3x3 rate matrix R with R_ij the rate from state j to state i.

    Columns sum to zero (probability conservation) and detailed balance holds
    with respect to `gibbs(model.energies, model.beta_bath)`.
    """
    e = np.asarray(model.energies, dtype=float)
    g12, g13, g23 = model.gamma
    bare = np.array(
        [
            [0.0, g12, g13],
            [g12, 0.0, g23],
            [g13, g23, 0.0],
        ]
    )

    # R_ij = Gamma_ij exp(-beta (E_i - E_j) / 2) for i != j.
    half_bias = np.exp(-model.beta_bath * (e[:, None] - e[None, :]) / 2.0)
    rate = bare * half_bias
    np.fill_diagonal(rate, 0.0)
    np.fill_diagonal(rate, -rate.sum(axis=0))
    return rate


# ---------------------------------------------------------------------------
# Spectral decomposition
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Spectrum:
    """Biorthonormal spectral data of a detailed-balance rate matrix.

    Attributes
    ----------
    eigenvalues:
        Real eigenvalues ordered lambda_1 = 0 > lambda_2 >= lambda_3.
    right:
        Right eigenvectors as columns: `right[:, k]` is v_{k+1}.
    left:
        Left eigenvectors as rows: `left[k, :]` is l_{k+1}.
    """

    eigenvalues: np.ndarray
    right: np.ndarray
    left: np.ndarray

    @property
    def gap_condition(self) -> bool:
        """True if the slowest decaying mode is non-degenerate (lambda_2 > lambda_3)."""
        return bool(self.eigenvalues[1] > self.eigenvalues[2] + 1e-12)


def spectral_decomposition(rate: np.ndarray, stationary: np.ndarray) -> Spectrum:
    """Diagonalise `rate` using its detailed-balance symmetrisation.

    With D = diag(sqrt(pi)), the matrix S = D^{-1} R D is symmetric.  Its
    orthonormal eigenvectors u_k give the right and left eigenvectors of R as
    v_k = D u_k and l_k = D^{-1} u_k, which are automatically biorthonormal.
    """
    sqrt_pi = np.sqrt(stationary)
    sym = rate * sqrt_pi[None, :] / sqrt_pi[:, None]
    sym = 0.5 * (sym + sym.T)  # discard round-off asymmetry

    vals, vecs = np.linalg.eigh(sym)
    order = np.argsort(vals)[::-1]  # lambda_1 = 0 first, then descending
    vals = vals[order]
    vecs = vecs[:, order]

    right = vecs * sqrt_pi[:, None]
    left = (vecs / sqrt_pi[:, None]).T

    # Fix the sign of each mode so that the result is reproducible: make the
    # first non-negligible component of every left eigenvector positive.
    for k in range(left.shape[0]):
        idx = np.argmax(np.abs(left[k]))
        if left[k, idx] < 0:
            left[k] *= -1.0
            right[:, k] *= -1.0

    return Spectrum(eigenvalues=vals, right=right, left=left)


def overlap_coefficients(spectrum: Spectrum, distribution: np.ndarray) -> np.ndarray:
    """Return the overlaps a_k = <l_k, p> for all modes k."""
    return spectrum.left @ distribution


def a2_of_beta(
    model: ThreeStateModel, spectrum: Spectrum, betas: np.ndarray
) -> np.ndarray:
    """Return a_2(beta_0) for an array of initial inverse temperatures."""
    return np.array(
        [overlap_coefficients(spectrum, gibbs(model.energies, b))[1] for b in betas]
    )


def a2_unnormalised(
    energies: Sequence[float], spectrum: Spectrum, betas: np.ndarray
) -> np.ndarray:
    """Return a_2(beta_0) up to a strictly positive factor, vectorised over `betas`.

    Since a_2(beta_0) = <l_2, pi(beta_0)> and pi(beta_0) = exp(-beta_0 E) / Z with
    Z > 0, the sign and the roots of a_2 are those of l_2 . exp(-beta_0 E).
    Dropping the normalisation turns the whole beta_0 sweep into a single matrix
    product, which is what makes the parameter scan below affordable.
    """
    e = np.asarray(energies, dtype=float)
    weights = np.exp(-np.outer(np.asarray(betas, dtype=float), e - e.min()))
    return weights @ spectrum.left[1]


# ---------------------------------------------------------------------------
# Time evolution
# ---------------------------------------------------------------------------


def evolve(
    spectrum: Spectrum,
    stationary: np.ndarray,
    initial: np.ndarray,
    times: np.ndarray,
) -> np.ndarray:
    """Propagate `initial` exactly via the spectral decomposition.

    Returns an array of shape (len(times), 3) holding p(t).  The k = 1 mode is
    the stationary state and is added separately, so the result is exact up to
    floating-point error rather than relying on a matrix exponential.
    """
    coeffs = overlap_coefficients(spectrum, initial)
    decay = np.exp(np.outer(times, spectrum.eigenvalues[1:]))
    return stationary[None, :] + decay * coeffs[1:][None, :] @ spectrum.right[:, 1:].T


def l1_distance(trajectory: np.ndarray, stationary: np.ndarray) -> np.ndarray:
    """Return ||p(t) - pi||_1 along a trajectory."""
    return np.abs(trajectory - stationary[None, :]).sum(axis=1)


# ---------------------------------------------------------------------------
# Search for a strong Mpemba effect
# ---------------------------------------------------------------------------


def strong_mpemba_root(
    model: ThreeStateModel,
    spectrum: Spectrum,
    beta_lo: float = 1e-3,
    beta_hi: float | None = None,
    n_scan: int = 4001,
) -> float | None:
    """Locate a root of a_2(beta_0) strictly between `beta_lo` and the bath value.

    a_2 vanishes identically at beta_0 = beta_bath.  A *second* root inside
    (0, beta_bath) is the signature of a strong Mpemba effect: initial states
    on either side of it have opposite signs of a_2, so |a_2| is necessarily
    non-monotonic there.  Returns the root, or None if there is no sign change.
    """
    if beta_hi is None:
        beta_hi = 0.999 * model.beta_bath

    betas = np.linspace(beta_lo, beta_hi, n_scan)
    values = a2_unnormalised(model.energies, spectrum, betas)

    sign_changes = np.nonzero(np.sign(values[:-1]) * np.sign(values[1:]) < 0)[0]
    if sign_changes.size == 0:
        return None

    # Bisection on the first bracketing interval.
    lo, hi = betas[sign_changes[0]], betas[sign_changes[0] + 1]
    f_lo = values[sign_changes[0]]
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        f_mid = a2_unnormalised(model.energies, spectrum, np.array([mid]))[0]
        if f_lo * f_mid <= 0:
            hi = mid
        else:
            lo, f_lo = mid, f_mid
    return 0.5 * (lo + hi)


def scan_for_mpemba(
    e3_values: np.ndarray,
    g13_values: np.ndarray,
    g23_values: np.ndarray,
    beta_bath: float,
    g12: float = 1.0,
) -> list[dict]:
    """Grid-scan model parameters for a strong Mpemba effect.

    E_1 = 0 and E_2 = 1 fix the energy origin and scale; the remaining
    parameters (E_3, Gamma_13, Gamma_23) are scanned with Gamma_12 held fixed.
    Returns the list of parameter sets that admit a strong Mpemba root,
    each annotated with the root location and the spectral gap ratio.
    """
    hits: list[dict] = []
    for e3 in e3_values:
        for g13 in g13_values:
            for g23 in g23_values:
                model = ThreeStateModel(
                    energies=(0.0, 1.0, float(e3)),
                    gamma=(g12, float(g13), float(g23)),
                    beta_bath=beta_bath,
                )
                pi = gibbs(model.energies, model.beta_bath)
                spectrum = spectral_decomposition(build_rate_matrix(model), pi)
                if not spectrum.gap_condition:
                    continue
                root = strong_mpemba_root(model, spectrum)
                if root is None:
                    continue
                hits.append(
                    {
                        "energies": model.energies,
                        "gamma": model.gamma,
                        "beta_bath": model.beta_bath,
                        "beta_star": root,
                        "lambda_2": float(spectrum.eigenvalues[1]),
                        "lambda_3": float(spectrum.eigenvalues[2]),
                        "gap_ratio": float(
                            spectrum.eigenvalues[2] / spectrum.eigenvalues[1]
                        ),
                    }
                )
    return hits
