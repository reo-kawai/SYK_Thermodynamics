"""Resonant-level model: dynamical mechanism of Mpemba-like temperature crossings.

This module solves the quench dynamics of a single fermionic level coupled to a
fermionic bath *exactly*, and extracts the effective inverse temperature with
the same fluctuation-dissipation (FDT) prescription that Wang, Su & Wang
[arXiv:2410.06669] apply to the strongly coupled SYK model.

The point of the exercise is that this model is quadratic: there is no
interaction, no disorder, no scrambling and no Lyapunov exponent.  A feature
of beta_eff(t) that survives here therefore cannot require quantum chaos.  The
calculation below establishes this for effective-temperature oscillations, but
not for trajectory crossings, whose separation stays below the fitting-window
systematic uncertainty.

Model
-----
    H = eps_d d^dag d + sum_k eps_k c_k^dag c_k + sum_k (V_k d^dag c_k + h.c.)

The bath enters only through the hybridisation function

    Gamma(omega) = pi sum_k |V_k|^2 delta(omega - eps_k)
                 = Gamma_0 Lambda^2 / (omega^2 + Lambda^2)      (Lorentzian)

whose Kramers-Kronig partner is, in closed form,

    Delta^R(omega) = Gamma_0 Lambda / (omega + i Lambda) .

Lambda is the bath bandwidth and sets the decay time tau_b ~ 1 / Lambda of the
retarded bath kernel.  For Lambda -> infinity, Delta^R -> -i Gamma_0, so that
the retarded equation becomes time-local.  The thermal lesser kernel still has
a scale set by beta_bath, and a GKSL/Lindblad reduction requires an additional
weak-coupling time-scale separation.  Finite Lambda introduces an explicit
memory kernel.  Sweeping (Gamma_0, Lambda) therefore interpolates continuously
between a wide-band weak-memory regime and a finite-memory strong-coupling one;
it does not by itself prove CP indivisibility.

Exact solution
--------------
Because Delta^R is a simple rational function, the retarded Green's function

    G^R(omega) = 1 / (omega - eps_d - Delta^R(omega))
               = (omega + i Lambda) / [(omega - z_1)(omega - z_2)]

has exactly two poles z_1, z_2 in the lower half plane, so the propagator
amplitude u(t) defined by G^R(t) = -i theta(t) u(t) is a sum of two
exponentials,

    u(t) = sum_j c_j exp(-i z_j t) ,   c_j = (z_j + i Lambda) / (z_j - z_{j'}) ,

with sum_j c_j = u(0) = 1.  This is evaluated analytically below; no ODE is
integrated.

The quench is of the standard partitioned type: for t < 0 the level is
decoupled and occupied with n_0 = 1 / (1 + exp(beta_0 eps_d)), the bath is in
equilibrium at beta_bath, and the coupling is switched on at t = 0.  The exact
two-time correlators then follow from the Cini formula

    G^<(t,t') = G^R(t,0) G^<(0,0) G^A(0,t')
              + int_0^t dt_1 int_0^{t'} dt_2 G^R(t,t_1) Sigma^<(t_1,t_2) G^A(t_2,t')

which, with Sigma^<(tau) = i int (domega/2pi) 2 Gamma(omega) f(omega)
exp(-i omega tau), reduces to

    G^<(t,t') = i [ n_0 u(t) u(t')^*
                    + int (domega/2pi) 2 Gamma(omega) f(omega) A(t,omega) A(t',omega)^* ]
    G^>(t,t') = -i [ (1 - n_0) u(t) u(t')^*
                    + int (domega/2pi) 2 Gamma(omega) (1-f(omega)) A(t,omega) A(t',omega)^* ]

with the memory kernel

    A(t,omega) = int_0^t ds u(t-s) exp(-i omega s)
               = sum_j c_j [exp(-i omega t) - exp(-i z_j t)] / (i (z_j - omega)) ,

again analytic.  The identity

    u(t) u(t')^* + int (domega/2pi) 2 Gamma(omega) A(t,omega) A(t',omega)^* = u(t-t')

for t > t' expresses completeness and is used as a numerical check.

Effective temperature
---------------------
Following the reference, the diagonal slice of the two-time correlator is
Wigner transformed over *positive relative time only*,

    G(omega,T) = int_0^infty dt' exp(i omega t') G(T + t'/2, T - t'/2) ,

which in a quench is necessarily cut off at t' = 2T because the dynamics start
at t = 0.  In equilibrium the FDT reads

    Im[G^>(omega) + G^<(omega)] / Im[G^R(omega)] = +tanh(beta omega / 2) ,

so beta_eff(T) is obtained from the slope of the left-hand side at omega = 0.
`effective_beta` extracts that slope by a least-squares fit through the origin
over a small symmetric frequency window; `EQUILIBRIUM_SIGN` records the
overall sign convention, which is fixed once and for all by
`selftest_equilibrium` rather than assumed.

Units
-----
eps_d = 1 sets the energy unit; Gamma_0, Lambda, omega and 1/beta are measured
in the same unit and times in its inverse.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np

# Sign and factor of the FDT slope.  Half-Fourier transforming
#     G^> - G^< = -i int (domega/2pi) A(omega) exp(-i omega tau)
#     G^> + G^< = -i int (domega/2pi) A(omega) tanh(beta omega / 2) exp(-i omega tau)
# over positive relative time gives Im -> -A(omega)/2 and
# -A(omega) tanh(beta omega / 2) / 2 respectively, so the ratio is
# +tanh(beta omega / 2) and its slope at omega = 0 is +beta/2.  Hence
# beta = 2 x slope.  Confirmed numerically by `selftest_equilibrium`.
EQUILIBRIUM_SIGN = 2.0


@dataclass(frozen=True)
class ResonantLevel:
    """Parameters of the resonant-level quench.

    Attributes
    ----------
    eps_d:
        Energy of the localised level (the energy unit; keep it non-zero, since
        at eps_d = 0 particle-hole symmetry pins the occupation to 1/2 and no
        temperature dynamics remains).
    gamma0:
        Coupling strength, i.e. the peak height of the hybridisation function.
    bandwidth:
        Lorentzian half-width Lambda of the bath band.  `numpy.inf` selects the
        strict wide-band limit of the retarded kernel.
    beta_bath:
        Inverse temperature of the bath.
    beta_init:
        Inverse temperature used to prepare the decoupled level at t < 0.
    """

    eps_d: float
    gamma0: float
    bandwidth: float
    beta_bath: float
    beta_init: float

    @property
    def initial_occupation(self) -> float:
        """Occupation of the decoupled level in equilibrium at `beta_init`."""
        return float(1.0 / (1.0 + np.exp(self.beta_init * self.eps_d)))

    def as_dict(self) -> dict:
        """Return a JSON-serialisable description, for output metadata."""
        d = asdict(self)
        d["bandwidth"] = None if np.isinf(self.bandwidth) else self.bandwidth
        return d


# ---------------------------------------------------------------------------
# Retarded propagator
# ---------------------------------------------------------------------------


def propagator_poles(model: ResonantLevel) -> tuple[np.ndarray, np.ndarray]:
    """Return the poles z_j and residue weights c_j of the retarded propagator.

    For finite bandwidth the denominator of G^R is the quadratic

        (omega - eps_d)(omega + i Lambda) - Gamma_0 Lambda ,

    whose two roots lie in the lower half plane.  In the wide-band limit the
    single pole eps_d - i Gamma_0 is returned with unit weight.
    """
    if np.isinf(model.bandwidth):
        return (
            np.array([model.eps_d - 1j * model.gamma0]),
            np.array([1.0 + 0.0j]),
        )

    lam = model.bandwidth
    # omega^2 + omega (i Lambda - eps_d) - (i eps_d Lambda + Gamma_0 Lambda)
    coeffs = [
        1.0 + 0.0j,
        1j * lam - model.eps_d,
        -(1j * model.eps_d * lam + model.gamma0 * lam),
    ]
    roots = np.roots(coeffs)
    weights = np.array(
        [
            (roots[j] + 1j * lam) / (roots[j] - roots[1 - j])
            for j in range(2)
        ]
    )
    return roots, weights


def u_amplitude(times: np.ndarray, poles: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Return u(t) = sum_j c_j exp(-i z_j t), with u(0) = 1."""
    return np.exp(-1j * np.outer(times, poles)) @ weights


def memory_kernel(
    times: np.ndarray,
    omegas: np.ndarray,
    poles: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    """Return A(t, omega) on a (len(times), len(omegas)) grid.

    A(t,omega) = sum_j c_j [exp(-i omega t) - exp(-i z_j t)] / (i (z_j - omega)).
    The removable singularity at omega = z_j never occurs for real omega because
    every pole has a strictly negative imaginary part.
    """
    t = times[:, None, None]
    w = omegas[None, :, None]
    z = poles[None, None, :]
    c = weights[None, None, :]
    num = np.exp(-1j * w * t) - np.exp(-1j * z * t)
    return (c * num / (1j * (z - w))).sum(axis=2)


# ---------------------------------------------------------------------------
# Bath functions
# ---------------------------------------------------------------------------


def hybridisation(omegas: np.ndarray, model: ResonantLevel) -> np.ndarray:
    """Return the Lorentzian hybridisation Gamma(omega)."""
    if np.isinf(model.bandwidth):
        return np.full_like(omegas, model.gamma0)
    lam = model.bandwidth
    return model.gamma0 * lam**2 / (omegas**2 + lam**2)


def fermi(omegas: np.ndarray, beta: float) -> np.ndarray:
    """Return the Fermi-Dirac occupation, evaluated without overflow.

    `numpy.where` would evaluate both branches and overflow on the discarded
    one, so the exponent is clipped instead.
    """
    x = np.clip(beta * np.asarray(omegas, dtype=float), -700.0, 700.0)
    return 0.5 * (1.0 - np.tanh(0.5 * x))


def frequency_grid(
    model: ResonantLevel, max_time: float, tolerance: float = 1e-6
) -> np.ndarray:
    """Return a symmetric frequency grid resolving both the bath and the phases.

    Two independent requirements fix the grid, and getting either wrong is
    silently fatal -- a range that is too small does not produce visible noise,
    it produces a smooth and completely wrong answer.

    * Range.  At large |omega| the bath integrand behaves as
      2 Gamma(omega) |A|^2 / 2pi ~ Gamma_0 Lambda^2 / (pi omega^4), since
      Gamma ~ Gamma_0 Lambda^2 / omega^2 and |A|^2 ~ 1 / omega^2.  Truncating at
      W therefore costs about Gamma_0 Lambda^2 / (3 pi W^3), so the cutoff is
      chosen as W = (Gamma_0 Lambda^2 / (3 pi tolerance))^{1/3}.  This 1/W^3 law
      was verified numerically against the completeness identity.
    * Spacing.  The kernel A(t, omega) carries the phase exp(-i omega t) with t
      up to `max_time`; d_omega = pi / (3 max_time) was found sufficient, the
      residual error being dominated by the range rather than the spacing.

    The wide-band limit is a special case: Gamma(omega) does not decay at all,
    the tail integral is only logarithmically controlled, and no finite grid is
    accurate.  Use a large finite `bandwidth` to represent the weak-memory regime
    instead; `solve_two_time` refuses the infinite case for this reason.
    """
    base = max(
        10.0 * abs(model.eps_d),
        10.0 * model.gamma0,
        20.0 / max(model.beta_bath, 1e-12),
        20.0 / max(model.beta_init, 1e-12),
        1.0,
    )
    tail = (
        model.gamma0 * model.bandwidth**2 / (3.0 * np.pi * tolerance)
    ) ** (1.0 / 3.0)
    scale = max(base, tail)

    d_omega = np.pi / (3.0 * max(max_time, 1.0))
    n_omega = int(2.0 * scale / d_omega) + 1
    n_omega = max(n_omega, 2001) | 1  # odd, so that omega = 0 is on the grid
    return np.linspace(-scale, scale, n_omega)


# ---------------------------------------------------------------------------
# Two-time correlators
# ---------------------------------------------------------------------------


@dataclass
class TwoTimeSolution:
    """Exact two-time correlators of the quenched resonant level.

    Attributes
    ----------
    times:
        Uniform time grid starting at t = 0.
    lesser, greater:
        Matrices G^<(t_i, t_j) and G^>(t_i, t_j).
    u:
        Propagator amplitude u(t_i).
    occupation:
        n(t) = -i G^<(t,t).
    completeness_error:
        Max violation of the completeness identity; a pure numerical check.
    """

    times: np.ndarray
    lesser: np.ndarray
    greater: np.ndarray
    u: np.ndarray
    occupation: np.ndarray
    completeness_error: float


def solve_two_time(
    model: ResonantLevel, times: np.ndarray, max_gigabytes: float = 2.0
) -> TwoTimeSolution:
    """Solve the quench exactly on the given uniform time grid.

    Raises if the required grid would exceed `max_gigabytes`, since the memory
    grows as (number of times) x (number of frequencies) and it is easy to
    request a combination of long times and wide bands that silently exhausts
    the machine.
    """
    if np.isinf(model.bandwidth):
        raise ValueError(
            "the wide-band limit has a non-decaying Gamma(omega) whose bath "
            "integral no grid resolves; use a large finite bandwidth instead"
        )

    poles, weights = propagator_poles(model)
    omegas = frequency_grid(model, float(times[-1]))

    needed = 16.0 * len(times) * len(omegas) / 1024**3
    if needed > max_gigabytes:
        raise MemoryError(
            f"kernel would need {needed:.1f} GB "
            f"({len(times)} times x {len(omegas)} frequencies); "
            "shorten the time grid or reduce the bandwidth"
        )

    d_omega = omegas[1] - omegas[0]

    u = u_amplitude(times, poles, weights)
    kernel = memory_kernel(times, omegas, poles, weights)

    gamma = hybridisation(omegas, model)
    f_bath = fermi(omegas, model.beta_bath)
    measure = 2.0 * gamma * d_omega / (2.0 * np.pi)

    # int (domega/2pi) 2 Gamma f A(t) A(t')^*  and its particle-hole partner.
    filled = (kernel * (measure * f_bath)) @ kernel.conj().T
    empty = (kernel * (measure * (1.0 - f_bath))) @ kernel.conj().T

    coherent = np.outer(u, u.conj())
    n0 = model.initial_occupation

    lesser = 1j * (n0 * coherent + filled)
    greater = -1j * ((1.0 - n0) * coherent + empty)

    # Completeness: u(t) u(t')^* + int (domega/2pi) 2 Gamma A A^* = u(t - t')
    total = coherent + (kernel * measure) @ kernel.conj().T
    lower = np.tril_indices(len(times))
    reference = u_amplitude(
        times[lower[0]] - times[lower[1]], poles, weights
    )
    completeness_error = float(np.max(np.abs(total[lower] - reference)))

    return TwoTimeSolution(
        times=times,
        lesser=lesser,
        greater=greater,
        u=u,
        occupation=np.real(-1j * np.diag(lesser)),
        completeness_error=completeness_error,
    )


# ---------------------------------------------------------------------------
# Wigner transform and effective temperature
# ---------------------------------------------------------------------------


def wigner_slice(
    solution: TwoTimeSolution,
    index: int,
    omegas: np.ndarray,
    rel_max: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Wigner transform the correlators at central time `solution.times[index]`.

    Returns (G^> + G^<)(omega, T) and G^R(omega, T).  On a uniform grid the
    Wigner coordinates (T + t'/2, T - t'/2) with t' = 2 k dt land exactly on the
    matrix entries (index + k, index - k), so no interpolation is needed.

    Three separate limits cap the relative time t':

    * t' <= 2T, the intrinsic restriction of a quench started at t = 0;
    * t' <= 2 (t_grid_max - T), an artefact of the finite grid, which is why the
      grid should extend well beyond the largest central time of interest;
    * t' <= `rel_max`, an explicit cap.

    Imposing the same `rel_max` at every central time matters for more than
    tidiness.  Truncating the half-Fourier transform at t'_max makes the result
    ring in frequency with period 2 pi / t'_max, so an unconstrained window
    would give every central time a different systematic error and the apparent
    time dependence of beta_eff would be dominated by that artefact.  Choose
    `rel_max` well beyond the correlator decay time 1 / Gamma_0.
    """
    dt = solution.times[1] - solution.times[0]
    k_max = min(index, len(solution.times) - 1 - index)
    if rel_max is not None:
        k_max = min(k_max, int(rel_max / (2.0 * dt)))
    if k_max < 2:
        return (
            np.full_like(omegas, np.nan, dtype=complex),
            np.full_like(omegas, np.nan, dtype=complex),
        )

    k = np.arange(k_max + 1)
    rel = 2.0 * k * dt
    d_rel = 2.0 * dt

    keldysh_rel = (
        solution.greater[index + k, index - k] + solution.lesser[index + k, index - k]
    )
    # G^R(t, t') = theta(t - t') [G^> - G^<](t, t'), which for this quadratic
    # model equals -i u(t - t'); both forms agree and the matrix form is used
    # so that the same numerical object enters numerator and denominator.
    retarded_rel = (
        solution.greater[index + k, index - k] - solution.lesser[index + k, index - k]
    )

    phase = np.exp(1j * np.outer(omegas, rel))
    weights = np.full(k.shape, d_rel)
    weights[0] *= 0.5  # trapezoid, half weight at t' = 0
    weights[-1] *= 0.5

    return phase @ (keldysh_rel * weights), phase @ (retarded_rel * weights)


def effective_beta(
    solution: TwoTimeSolution,
    index: int,
    omega_window: float,
    rel_max: float | None = None,
    n_omega: int = 41,
) -> float:
    """Extract beta_eff at one central time by fitting the FDT form.

    Half-Fourier transformed, the equilibrium FDT reads

        Im[G^> + G^<](omega, T) / Im[G^R](omega, T) = tanh(beta omega / 2) ,

    and beta_eff is defined as the best fit of that form over a symmetric
    frequency window -- the same "best fit to the FDT" that the reference uses,
    and the reason its beta_eff can come out negative far from equilibrium.

    Fitting tanh rather than its slope matters: a straight line through the
    origin systematically underestimates beta by about (beta W)^2 / 20 over a
    window of half-width W, which is a per-cent-level bias at the parameters
    used here and would be mistaken for physics.
    """
    omegas = np.linspace(-omega_window, omega_window, n_omega)
    omegas = omegas[np.abs(omegas) > 1e-12]
    keldysh, retarded = wigner_slice(solution, index, omegas, rel_max=rel_max)
    if not np.all(np.isfinite(keldysh)):
        return np.nan

    ratio = np.imag(keldysh) / np.imag(retarded)

    # Linear estimate first; it is the exact answer as the window shrinks and
    # provides a starting point that keeps the tanh fit out of spurious minima.
    slope = float(np.dot(omegas, ratio) / np.dot(omegas, omegas))
    beta0 = EQUILIBRIUM_SIGN * slope

    from scipy.optimize import least_squares

    fit = least_squares(
        lambda p: np.tanh(0.5 * p[0] * omegas) - ratio,
        x0=[beta0],
        method="lm",
    )
    return float(fit.x[0])


def distribution_function(
    solution: TwoTimeSolution,
    index: int,
    omegas: np.ndarray,
    rel_max: float | None = None,
) -> np.ndarray:
    """Return the nonequilibrium distribution f_eff(omega, T).

    Half-Fourier transforming G^< gives Im F[G^<] = A(omega) f(omega) / 2, while
    Im F[G^> - G^<] = Im G^R = -A(omega) / 2, so

        f_eff(omega, T) = -Im[G^<](omega, T) / Im[G^R](omega, T) .

    Unlike beta_eff this is a property of the state at each frequency, with no
    fitting window and no assumed functional form.
    """
    dt = solution.times[1] - solution.times[0]
    k_max = min(index, len(solution.times) - 1 - index)
    if rel_max is not None:
        k_max = min(k_max, int(rel_max / (2.0 * dt)))
    if k_max < 2:
        return np.full_like(omegas, np.nan, dtype=float)

    k = np.arange(k_max + 1)
    rel = 2.0 * k * dt
    weights = np.full(k.shape, 2.0 * dt)
    weights[0] *= 0.5
    weights[-1] *= 0.5
    phase = np.exp(1j * np.outer(omegas, rel))

    lesser_rel = solution.lesser[index + k, index - k]
    retarded_rel = (
        solution.greater[index + k, index - k] - solution.lesser[index + k, index - k]
    )
    lesser_w = phase @ (lesser_rel * weights)
    retarded_w = phase @ (retarded_rel * weights)
    return -np.imag(lesser_w) / np.imag(retarded_w)


def distance_to_bath(
    solution: TwoTimeSolution,
    index: int,
    model: ResonantLevel,
    omegas: np.ndarray,
    rel_max: float | None = None,
) -> float:
    """Return the spectrally weighted distance of f_eff from the bath function,

        D_bath(T) = int d_omega A_+(omega,T) |f_eff - f_bath| / int d_omega A_+(omega,T) .

    This is the resonant-level analogue of the distance function that the
    reference reports as decreasing monotonically while its beta_eff
    oscillates.  Comparing the two is the whole point: one is a property of the
    state, the other of the thermometer.

    The nonnegative spectral weight A_+(omega,T) = max(-2 Im G^R, 0) is
    essential rather than cosmetic.  A finite-time Wigner transform can have
    small negative spectral lobes, which are clipped exactly as in the code
    below.  For a narrow band the level has almost no spectral weight at
    |omega| >> Lambda, where nothing constrains f_eff and the bath cannot
    thermalise anything; an unweighted integral is then dominated by that
    empty region and saturates at a large value that says nothing about
    relaxation.
    """
    f_eff = distribution_function(solution, index, omegas, rel_max=rel_max)
    if not np.all(np.isfinite(f_eff)):
        return np.nan

    _, retarded = wigner_slice(solution, index, omegas, rel_max=rel_max)
    spectral = -2.0 * np.imag(retarded)
    spectral = np.clip(spectral, 0.0, None)
    norm = float(np.trapezoid(spectral, omegas))
    if norm <= 0.0:
        return np.nan

    f_bath = fermi(omegas, model.beta_bath)
    return float(np.trapezoid(spectral * np.abs(f_eff - f_bath), omegas) / norm)


def effective_beta_trace(
    solution: TwoTimeSolution,
    omega_window: float,
    rel_max: float | None = None,
    stride: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (central times, beta_eff) over the whole trajectory."""
    indices = np.arange(0, len(solution.times), stride)
    betas = np.array(
        [
            effective_beta(solution, int(i), omega_window, rel_max=rel_max)
            for i in indices
        ]
    )
    return solution.times[indices], betas


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------


def selftest_equilibrium(tolerance: float = 0.05) -> dict:
    """Check the machinery against a case whose answer is known.

    At weak coupling the level equilibrates with the bath, so beta_eff at late
    times must approach beta_bath.  This simultaneously fixes the sign
    convention `EQUILIBRIUM_SIGN` and validates the Wigner/FDT extraction.
    """
    model = ResonantLevel(
        eps_d=1.0, gamma0=0.1, bandwidth=10.0, beta_bath=1.0, beta_init=4.0
    )
    times = np.linspace(0.0, 200.0, 2001)
    solution = solve_two_time(model, times)
    index = len(times) // 2
    beta_late = effective_beta(solution, index, omega_window=1.0, rel_max=80.0)
    return {
        "completeness_error": solution.completeness_error,
        "beta_late": beta_late,
        "beta_bath": model.beta_bath,
        "relative_error": abs(beta_late - model.beta_bath) / model.beta_bath,
        "passed": bool(
            solution.completeness_error < 1e-5
            and abs(beta_late - model.beta_bath) < tolerance * model.beta_bath
        ),
    }
