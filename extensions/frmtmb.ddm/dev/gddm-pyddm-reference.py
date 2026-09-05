# /// script
# requires-python = "==3.13.*"
# dependencies = [
#     "pyddm==0.9.0",
#     "numpy==2.5.2",
#     "scipy==1.18.1",
# ]
# ///
"""Freeze a PyDDM reference for the generalized components of gddm().

Run once, by hand, to regenerate the fixture that
tests/testthat/test-gddm-reference.R reads:

    uv run extensions/frmtmb.ddm/dev/gddm-pyddm-reference.py

ON WINDOWS, point UV_CACHE_DIR somewhere short first, for example
C:\\uvcache. uv names the environment after this script, and under a
deep cache directory the path of one scipy file inside it reaches 260
characters, which is MAX_PATH: the file installs, listdir shows it, and
the import system cannot see it, so scipy fails to import a module that
is on disk. A short cache directory is the whole fix.

Nothing in the package's tests runs Python. This writes three CSVs under
tests/testthat/fixtures/ and the test asserts their metadata, so a
regenerated fixture that moved cannot pass silently.

WHY PyDDM. gddm()'s generalized components, leak, exponential and linear
collapse, coherence drift, start-point variability and lapse, are
checked by nothing outside the solver itself. The analytic Wiener check
covers only a constant drift and fixed bounds, the gradient tests prove
derivative consistency rather than correctness, and parameter recovery
is circular because gddm_simulate() draws from the solver's own density
by design. PyDDM (Shinn, Lam and Murray 2020) is the canonical
implementation of exactly those components, and it is an unrelated
discretization: a fixed grid in the original coordinate with the moving
bound sandwiched between the two nodes around it, where gddm() changes
variable so that the bound is pinned and the grid never moves.

WHY DETERMINISM. The fixture is frozen, so nothing here may draw a
random number, read a clock, or iterate to a tolerance that depends on
machine arithmetic. Every grid, parameter and dataset below is written
out; the only search is the resolution ladder, which is a fixed list.
"""

from __future__ import annotations

import hashlib
import importlib.metadata as md
import platform
import sys
from pathlib import Path

import numpy as np
import pyddm

# --------------------------------------------------------------------
# Conventions shared by every case
# --------------------------------------------------------------------

# The modeled window. gddm()'s renormalization and its lapse floor are
# both statements about t_max, so this number is part of the model and
# the test asserts the fixture carries it.
T_DUR = 3.0

# Response times are compared on this grid. The step is a multiple of
# every dt in the ladder below AND of the gddm() steps the test uses,
# so neither side has to interpolate a curve it did not compute.
GRID_STEP = 0.02

# Where the comparison starts, as a decision time. Both solvers are
# wrong in the leading edge for the same structural reason, an implicit
# scheme spreads a little mass everywhere at once where the true
# first-passage density is exponentially small, so a comparison that
# started at zero would measure a shared artefact rather than either
# solver. 0.2 s is what gddm_control()'s help page already claims
# accuracy from.
GRID_FROM = 0.20
GRID_TO = 2.00

# Grid points where the reference density at either wall falls below
# this are dropped. Below it the density carries no likelihood weight
# in any realistic dataset, both solvers' RELATIVE error grows without
# bound, and PyDDM stops stepping entirely once the surviving mass is
# under 1e-4, leaving exact zeros.
DENS_FLOOR = 1e-6

# The reference's own granted tolerance, in log density. The ladder
# stops at the first level whose change from the previous level is
# under this, and the level it stopped at and the change it achieved
# are both written to the fixture.
REF_TOL = 5e-4

# The resolution ladder, halving both steps together. Every dt divides
# GRID_STEP and every non-decision time used below; every dx divides
# every bound half-separation and every start point used below. That is
# not decoration: PyDDM rounds the bound UP to a multiple of dx
# (Model.x_domain), snaps a point start to round(x0/dx) (ICPoint), and
# shifts by int(ndt/dt) whole bins (OverlayNonDecision), so a step that
# does not divide these solves a quietly different model.
LADDER = [
    (1.25e-3, 2.5e-3),
    (6.25e-4, 1.25e-3),
    (3.125e-4, 6.25e-4),
    (1.5625e-4, 3.125e-4),
    (7.8125e-5, 1.5625e-4),
]

FIXTURES = (Path(__file__).resolve().parents[1]
            / "tests" / "testthat" / "fixtures")
PREFIX = "gddm-pyddm-"


# --------------------------------------------------------------------
# The case matrix
# --------------------------------------------------------------------
#
# Each case is written in gddm()'s OWN parameters. The translation to
# PyDDM happens in one place, build_model(), so the mapping the test
# documents and the mapping the fixture was built with are the same
# text.
#
#   gddm()                                  PyDDM
#   dx = a(x,t) dt + dW                     NoiseConstant(noise=1)
#   bounds at +/- B(t), B(0) = bs/2         B = bs/2
#   bias in (0,1) above the lower bound     ICPoint(x0=(bs/2)(2 bias-1))
#   a = mu                                  DriftConstant(drift=mu)
#   a = mu - leak*x                         DriftLinear(drift=mu, x=-leak, t=0)
#   B = (bs/2) exp(-t/tau)                  BoundCollapsingExponential(
#                                             B=bs/2, tau=1/tau)
#   B = (bs/2)(1 - kappa t/t_max)           BoundCollapsingLinear(
#                                             B=bs/2, t=(bs/2) kappa/t_max)
#   a = sign(C) mu (|C|/cmax)^alpha         DriftConstant(drift=that number)
#   uniform start, half width sz            ICRange(sz=sz*bs/2)
#   ndt                                     OverlayNonDecision(nondectime=ndt)
#
# The leak sign is the trap the help page warns about: gddm()'s `leak`
# is the paper's l, POSITIVE for leaky integration, and PyDDM's
# DriftLinear x coefficient is its negative.

CASES = [
    # the calibration case: constant drift, fixed bounds, unbiased start
    dict(id="constant", mu=1.2, bs=2.0, bias=0.5, ndt=0.25),
    # a biased point start and a long non-decision time
    dict(id="bias_ndt", mu=0.8, bs=2.0, bias=0.65, ndt=0.40),
    # leaky integration: positive leak pulls the accumulator back
    dict(id="leak", mu=1.5, leak=1.2, bs=2.0, bias=0.5, ndt=0.25),
    # unstable integration: the other sign, which a sign error hides
    dict(id="unstable", mu=0.6, leak=-0.8, bs=2.0, bias=0.5, ndt=0.25),
    dict(id="exp_collapse", mu=1.2, bs=2.0, tau=1.5, bias=0.5, ndt=0.25),
    dict(id="lin_collapse", mu=1.2, bs=2.0, kappa=0.6, bias=0.5, ndt=0.25),
    # the coherence nonlinearity at cmax and at cmax/2, which is the
    # pair that pins alpha
    dict(id="coh_high", mu=2.0, alpha=1.3, cmax=0.512, coh=0.512,
         bs=2.0, bias=0.5, ndt=0.25),
    dict(id="coh_low", mu=2.0, alpha=1.3, cmax=0.512, coh=0.256,
         bs=2.0, bias=0.5, ndt=0.25),
    dict(id="sz_uniform", mu=1.0, bs=2.0, bias=0.5, sz=0.2, ndt=0.25),
    dict(id="lapse", mu=1.2, bs=2.0, bias=0.5, ndt=0.25, lapse=0.05),
]


def drift_value(c):
    """The constant part of the drift, in gddm()'s parameters."""
    if "coh" in c:
        cc, cmax = c["coh"], c["cmax"]
        return np.sign(cc) * c["mu"] * (abs(cc) / cmax) ** c["alpha"]
    return c["mu"]


def build_model(c, dt, dx, with_lapse=False):
    half = c["bs"] / 2.0
    if "leak" in c:
        drift = pyddm.DriftLinear(drift=drift_value(c), x=-c["leak"], t=0.0)
    else:
        drift = pyddm.DriftConstant(drift=drift_value(c))
    if "tau" in c:
        bound = pyddm.BoundCollapsingExponential(B=half, tau=1.0 / c["tau"])
    elif "kappa" in c:
        bound = pyddm.BoundCollapsingLinear(
            B=half, t=half * c["kappa"] / T_DUR)
    else:
        bound = pyddm.BoundConstant(B=half)
    if "sz" in c:
        # ICRange is centred on the domain and carries no bias, so a
        # case that uses it must be unbiased. The point-start bias is
        # carried by bias_ndt instead.
        assert c["bias"] == 0.5, "ICRange cannot express a biased start"
        ic = pyddm.ICRange(sz=c["sz"] * half)
    else:
        ic = pyddm.ICPoint(x0=half * (2.0 * c["bias"] - 1.0))
    ov = pyddm.OverlayNonDecision(nondectime=c["ndt"])
    if with_lapse:
        ov = pyddm.OverlayChain(overlays=[
            ov, pyddm.OverlayUniformMixture(umixturecoef=c["lapse"])])
    return pyddm.Model(drift=drift, noise=pyddm.NoiseConstant(noise=1.0),
                       bound=bound, IC=ic, overlay=ov,
                       dx=dx, dt=dt, T_dur=T_DUR)


def solver_name(c):
    """Which PyDDM solver this case gets, and why.

    Analytic wherever PyDDM has one, which is a constant bound or a
    linear collapse with a point start (Anderson 1960). Otherwise
    Crank-Nicolson, which PyDDM offers only for a bound that does not
    move: Model.can_solve_cn() returns FALSE the moment the bound
    depends on time. The exponential collapse therefore gets backward
    Euler, which is FIRST order in dt, and its ladder has to run
    further as a result.
    """
    m = build_model(c, LADDER[0][0], LADDER[0][1])
    if m.has_analytical_solution():
        return "analytical"
    return "cn" if m.can_solve_cn() else "implicit"


def solve_raw(c, dt, dx, how, with_lapse=False):
    """Defective per-bin masses at the two walls, on PyDDM's own grid."""
    m = build_model(c, dt, dx, with_lapse=with_lapse)
    if how == "analytical":
        s = m.solve_analytical()
    elif how == "cn":
        s = m.solve_numerical_cn()
    else:
        s = m.solve_numerical(method="implicit")
    return s, m


def on_grid(s, m, grid):
    """Densities at the comparison times, read without interpolation."""
    idx = np.rint(grid / m.dt).astype(int)
    assert np.allclose(m.t_domain()[idx], grid, atol=1e-9), \
        "the comparison grid is not a subset of PyDDM's time grid"
    return s.pdf("correct")[idx], s.pdf("error")[idx]


def case_grid(c):
    n = int(round((GRID_TO - GRID_FROM) / GRID_STEP))
    return np.round(c["ndt"] + GRID_FROM + GRID_STEP * np.arange(n + 1), 10)


def ll_times(c, grid):
    """The fixed dataset's response times.

    Deliberately ODD multiples of 0.005 so that they fall BETWEEN the
    nodes of every gddm() grid the test uses, which are multiples of
    0.01. The curve comparison reads gddm()'s density at its own nodes;
    this reads it through the interpolation the likelihood actually
    performs, which is a different code path.
    """
    lo, hi = grid[0], grid[-1]
    want = lo + (hi - lo) * np.array([0.05, 0.18, 0.31, 0.44, 0.60, 0.80])
    t = np.round(np.round((want - 0.005) / 0.01) * 0.01 + 0.005, 10)
    t = np.unique(t)
    assert np.all(t > lo - 1e-12) and np.all(t < hi + 1e-12)
    return t


# --------------------------------------------------------------------
# Reference values for one case
# --------------------------------------------------------------------

def reference(c):
    how = solver_name(c)
    grid = case_grid(c)
    lt = ll_times(c, grid)
    both = np.concatenate([grid, lt])

    moments = {}

    def at(dt, dx):
        s, m = solve_raw(c, dt, dx, how)
        up, lo = on_grid(s, m, both)
        # Conditional mean DECISION time at each wall, over the whole
        # window rather than the comparison grid. The Euler-Maruyama
        # cross-check in the R test needs a functional both solvers and
        # a path simulation can produce, and this is one. The trapezoid
        # rule, not the sum, because that is the mass of the piecewise
        # linear density a likelihood interpolating between nodes
        # actually reads.
        td = m.t_domain() - c["ndt"]
        for nm, v in (("upper", s.choice_upper), ("lower", s.choice_lower)):
            w = v.copy()
            w[0] *= 0.5
            w[-1] *= 0.5
            moments["mean_dt_" + nm] = float(np.sum(td * w) / np.sum(w))
        # gddm() divides each condition's pair of densities by their own
        # total mass, so its likelihood is explicitly conditional on a
        # response inside the window (gd_densities(), renormalize).
        # PyDDM keeps the undecided mass instead. Renormalizing here is
        # what makes the two comparable, and the undecided mass is
        # written to the fixture rather than hidden.
        tot = s.prob("correct") + s.prob("error")
        return np.array([up, lo]) / tot, tot, s.prob("correct") / tot

    prev = None
    residual = float("nan")
    used = 0
    cur = None
    for level, (dt, dx) in enumerate(LADDER):
        cur = at(dt, dx)
        used = level
        if prev is not None:
            k = (cur[0][0] > DENS_FLOOR) & (cur[0][1] > DENS_FLOOR)
            residual = float(np.max(np.abs(
                np.log(cur[0][:, k]) - np.log(prev[0][:, k]))))
            if residual < REF_TOL:
                break
        prev = cur
    dens, tot, p_up = cur
    dt, dx = LADDER[used]

    ng = len(grid)
    keep = (dens[0][:ng] > DENS_FLOOR) & (dens[1][:ng] > DENS_FLOOR)
    d_grid = dens[:, :len(grid)][:, keep]
    d_ll = dens[:, len(grid):]

    meta = dict(solver=how, level=used, dt=dt, dx=dx, T_dur=T_DUR,
                ref_residual_dlog=residual, prob_upper=p_up,
                prob_lower=1.0 - p_up, mass_decided=tot,
                prob_undecided=1.0 - tot,
                n_grid=int(keep.sum()), t_first=float(grid[keep][0]),
                t_last=float(grid[keep][-1]), grid_step=GRID_STEP,
                dens_floor=DENS_FLOOR, **moments)

    # Where PyDDM also has a closed form, say how far its own numerics
    # are from it. That number is the honest size of the discretization
    # this reference is made of, and for the linear collapse it is where
    # PyDDM's grid-snapped moving bound shows up.
    if how == "analytical":
        alt = "cn" if build_model(c, dt, dx).can_solve_cn() else "implicit"
        sn, mn = solve_raw(c, dt, dx, alt)
        un, ln = on_grid(sn, mn, grid)
        tn = sn.prob("correct") + sn.prob("error")
        num = np.array([un, ln]) / tn
        meta["pyddm_numeric_vs_analytic_dlog"] = float(np.max(np.abs(
            np.log(num[:, keep]) - np.log(dens[:, :len(grid)][:, keep]))))

    if "lapse" in c:
        # gddm() renormalizes and THEN mixes, with a lapse density of
        # 0.5*lapse/t_max at each wall (gd_densities()). PyDDM's
        # OverlayUniformMixture adds 0.5*c/len(t_domain) of MASS to each
        # bin, so as a density its floor is 0.5*c/(T_dur + dt): it
        # spreads the lapse over one bin more than the window, and it
        # mixes into a defective density rather than a normalized one.
        # The mixture FORM is checked against PyDDM's implementation
        # here; the fixture then carries gddm()'s convention.
        lp = c["lapse"]
        sl, ml = solve_raw(c, dt, dx, how, with_lapse=True)
        ul, ll_ = on_grid(sl, ml, grid)
        s0, m0 = solve_raw(c, dt, dx, how)
        u0, l0 = on_grid(s0, m0, grid)
        nbin = len(m0.t_domain())
        manual = np.array([u0, l0]) * (1 - lp) + 0.5 * lp / (nbin * m0.dt)
        meta["lapse_formula_max_abs_diff"] = float(
            np.max(np.abs(manual - np.array([ul, ll_]))))
        meta["lapse_pyddm_density"] = 0.5 * lp / (nbin * m0.dt)
        meta["lapse_gddm_density"] = 0.5 * lp / T_DUR
        d_grid = d_grid * (1 - lp) + 0.5 * lp / T_DUR
        d_ll = d_ll * (1 - lp) + 0.5 * lp / T_DUR

    # The fixed dataset: the first response time at the upper wall, the
    # next at the lower, alternating, so both walls carry likelihood.
    upper = np.arange(len(lt)) % 2 == 0
    ll_dens = np.where(upper, d_ll[0], d_ll[1])
    meta["loglik"] = float(np.sum(np.log(ll_dens)))
    meta["n_obs"] = int(len(lt))
    return grid[keep], d_grid, lt, upper, ll_dens, meta


# --------------------------------------------------------------------
# Write it out
# --------------------------------------------------------------------

def fmt(x):
    return "%.12e" % x


def main():
    FIXTURES.mkdir(parents=True, exist_ok=True)
    dens_rows = ["case,t,dens_upper,dens_lower"]
    ll_rows = ["case,rt,upper,dens"]
    meta_rows = ["case,key,value"]
    for c in CASES:
        g, d, lt, upper, ld, meta = reference(c)
        for i in range(len(g)):
            dens_rows.append("%s,%s,%s,%s" % (c["id"], fmt(g[i]),
                                              fmt(d[0][i]), fmt(d[1][i])))
        for i in range(len(lt)):
            ll_rows.append("%s,%s,%d,%s" % (c["id"], fmt(lt[i]),
                                            int(upper[i]), fmt(ld[i])))
        # every gddm() parameter, so the test can assert it is building
        # the model the fixture was built from
        for k in sorted(c):
            if k == "id":
                continue
            meta_rows.append("%s,par_%s,%s" % (c["id"], k, fmt(c[k])))
        meta_rows.append("%s,drift_value,%s" % (c["id"], fmt(drift_value(c))))
        for k in sorted(meta):
            v = meta[k]
            meta_rows.append("%s,%s,%s" % (
                c["id"], k, v if isinstance(v, str) else fmt(v)))
        print("%-13s %-11s level %d dt=%g dx=%g residual=%.2e n=%d"
              % (c["id"], meta["solver"], meta["level"], meta["dt"],
                 meta["dx"], meta["ref_residual_dlog"], meta["n_grid"]))

    dens_txt = "\n".join(dens_rows) + "\n"
    ll_txt = "\n".join(ll_rows) + "\n"
    # The hash covers the payload only. The test recomputes it, so a
    # regenerated fixture whose numbers moved cannot slip past the
    # metadata assertions.
    h = hashlib.sha256((dens_txt + ll_txt).encode("utf-8")).hexdigest()
    for k, v in [("pyddm", md.version("pyddm")),
                 ("numpy", md.version("numpy")),
                 ("scipy", md.version("scipy")),
                 ("python", platform.python_version()),
                 ("T_dur", fmt(T_DUR)), ("grid_step", fmt(GRID_STEP)),
                 ("grid_from", fmt(GRID_FROM)), ("grid_to", fmt(GRID_TO)),
                 ("dens_floor", fmt(DENS_FLOOR)), ("ref_tol", fmt(REF_TOL)),
                 ("n_cases", str(len(CASES))), ("content_sha256", h)]:
        meta_rows.append("_fixture,%s,%s" % (k, v))
    (FIXTURES / (PREFIX + "density.csv")).write_text(
        dens_txt, encoding="utf-8", newline="\n")
    (FIXTURES / (PREFIX + "loglik.csv")).write_text(
        ll_txt, encoding="utf-8", newline="\n")
    (FIXTURES / (PREFIX + "meta.csv")).write_text(
        "\n".join(meta_rows) + "\n", encoding="utf-8", newline="\n")
    print("sha256", h)
    print("wrote", FIXTURES)
    return 0


if __name__ == "__main__":
    sys.exit(main())
