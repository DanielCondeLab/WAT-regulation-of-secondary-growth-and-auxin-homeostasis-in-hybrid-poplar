"""
Intracellular auxin rheostat ODE model for auxin homeostasis in poplar 
Model consts of four compartments: Apoplast, Cytoplasm, Vacuole.

"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from dataclasses import dataclass
import csv
import os
import warnings


# ---IaaH protonated fraction---
pKa_IAA = 4.75

def f_IAAH(pH: float) -> float:
    return 1.0 / (1.0 + 10.0 ** (pH - pKa_IAA))


def mm(Vmax: float, Km: float, S: float) -> float:
    return Vmax * S / (Km + S)

# ---Model global parameters---

SCALE = 0.00091

V_LAX  = SCALE * 2.97
V_PIN  = SCALE * 3.77
V_ABCB = SCALE * 2.77
V_WAT1 = SCALE * 6.96

# calibrated rate constants  
K_INACT = 0.00458    # active to inactive conversion in vacuole
K_DEG   = 0.00042    # degradation of inactive in vacuole
K_GH3   = 0.00011    # GH3 conjugation in cytoplasm
SYNTH   = 9.748e-6   # synthesis rate

# WAT1 Km for inactive auxin conjugates
KM_WAT1_INACT = 0.62   # uM

# experimental references (nM)
WT_ACTIVE   = 105.0
WT_INACTIVE = 10.8
WAT1_ACTIVE   = 8.4
WAT1_INACTIVE = 22.0

@dataclass
class Params:
    #Auxin transport model with active/inactive pools.
    Vmax_LAX:      float = V_LAX
    Km_LAX:        float = 1.0
    Vmax_PIN:      float = V_PIN
    Km_PIN:        float = 10.0
    Vmax_ABCB_PM:  float = V_ABCB
    Km_ABCB_PM:    float = 1.5
    P_IAH_PM:      float = 0.0008
    Vmax_ABCB_vac: float = V_ABCB
    Km_ABCB_vac:   float = 1.0
    Vmax_WAT1:     float = V_WAT1
    Km_WAT1:       float = KM_WAT1_INACT
    P_IAH_tono:    float = 0.00002
    pH_apo:        float = 5.3
    pH_cyt:        float = 7.3
    pH_vac:        float = 5.5
    k_inact:       float = K_INACT
    k_GH3:         float = K_GH3
    k_deg:         float = K_DEG
    synth:         float = SYNTH


# ---ODEs section---

def ode(t: float, y: np.ndarray, p: Params) -> np.ndarray:
    """
    State y = [A_cyt, A_vac_act, A_vac_inact, A_apo]  (all uM).
    WAT1 exports inactive auxin which becomes active in cytoplasm.
    """
    C, Va, Vi, A = y
    C = max(C, 0); Va = max(Va, 0); Vi = max(Vi, 0); A = max(A, 0)

    f_apo = f_IAAH(p.pH_apo)
    f_cyt = f_IAAH(p.pH_cyt)
    f_vac = f_IAAH(p.pH_vac)

    # PM transport fluxes
    J_LAX     = mm(p.Vmax_LAX,      p.Km_LAX,      A)
    J_PIN     = mm(p.Vmax_PIN,      p.Km_PIN,      C)
    J_ABCB_PM = mm(p.Vmax_ABCB_PM,  p.Km_ABCB_PM,  C)
    J_pass_PM = p.P_IAH_PM * (A * f_apo - C * f_cyt)

    # tonoplast transport fluxes
    J_ABCB_vac = mm(p.Vmax_ABCB_vac, p.Km_ABCB_vac, C)
    J_pass_ton = p.P_IAH_tono * (C * f_cyt - Va * f_vac)

    # metabolic fluxes
    J_inact = p.k_inact * Va               # active2inactive  (in vacuole)
    J_WAT1  = mm(p.Vmax_WAT1, p.Km_WAT1, Vi)  # WAT1 exports INACTIVE from vac
    # WAT1 flux delivers inactive auxin to cytoplasm and becomes active
    J_GH3 = p.k_GH3 * C                    # GH3 conjugation (cytoplasm)
    J_deg = p.k_deg * Vi                   # degradation of inactive (vacuole)

    # ODE RHS
    dC = (+ J_LAX - J_PIN - J_ABCB_PM + J_pass_PM
          - J_ABCB_vac + J_WAT1 - J_pass_ton
          + p.synth - J_GH3)

    dVa = (+ J_ABCB_vac + J_pass_ton - J_inact)

    dVi = (+ J_inact - J_WAT1 - J_deg)

    dA = (- J_LAX + J_PIN + J_ABCB_PM - J_pass_PM)

    return np.array([dC, dVa, dVi, dA])


# ---Solver section---
def steady_state(p: Params, guess: np.ndarray = None) -> np.ndarray:
    """Find [A_cyt, A_vac_act, A_vac_inact, A_apo] at steady state."""
    if guess is not None:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            sol = fsolve(lambda y: ode(0, y, p), guess,
                         xtol=1e-10, maxfev=4000)
        return np.maximum(sol, 0.0)

    # try multiple initial guesses, pick best (lowest residual)
    guesses = [
        np.array([0.01, 0.005, 0.005, 0.05]),
        np.array([0.05, 0.02, 0.01, 0.05]),
        np.array([0.001, 0.001, 0.02, 0.005]),
        np.array([0.02, 0.01, 0.01, 0.02]),
    ]
    best_sol, best_res = None, np.inf
    for g in guesses:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)
                sol = fsolve(lambda y: ode(0, y, p), g,
                             xtol=1e-10, maxfev=4000)
            sol = np.maximum(sol, 0.0)
            res = np.sum(np.abs(ode(0, sol, p)))
            if res < best_res:
                best_res = res
                best_sol = sol
        except Exception:
            continue
    if best_sol is None:
        best_sol = np.array([0.0, 0.0, 0.0, 0.0])
    return best_sol

# ---Simulation routine---
def simulate(p: Params, y0: np.ndarray, t_end: float,
             n: int = 500) -> dict:
    t_eval = np.linspace(0, t_end, n)
    sol = solve_ivp(ode, (0, t_end), y0, args=(p,),
                    method='LSODA', t_eval=t_eval,
                    rtol=1e-9, atol=1e-12)
    return {
        't':    sol.t,
        'C':    sol.y[0],   # A_cyt
        'Va':   sol.y[1],   # A_vac_active
        'Vi':   sol.y[2],   # A_vac_inactive
        'A':    sol.y[3],   # A_apo
    }

# ---Flux definitions---
def fluxes_at_state(C, Va, Vi, A, p: Params) -> dict:
    """Return all individual fluxes (uM/s) for a given 4-state vector."""
    f_apo = f_IAAH(p.pH_apo)
    f_cyt = f_IAAH(p.pH_cyt)
    f_vac = f_IAAH(p.pH_vac)
    return {
        'J_LAX':        mm(p.Vmax_LAX,      p.Km_LAX,      A),
        'J_PIN':        mm(p.Vmax_PIN,      p.Km_PIN,      C),
        'J_ABCB_PM':    mm(p.Vmax_ABCB_PM,  p.Km_ABCB_PM,  C),
        'J_pass_PM':    p.P_IAH_PM * (A * f_apo - C * f_cyt),
        'J_ABCB_vac':   mm(p.Vmax_ABCB_vac, p.Km_ABCB_vac, C),
        'J_pass_ton':   p.P_IAH_tono * (C * f_cyt - Va * f_vac),
        'J_inact':      p.k_inact * Va,
        'J_WAT1':       mm(p.Vmax_WAT1,     p.Km_WAT1,     Vi),
        'J_GH3':        p.k_GH3 * C,
        'J_deg':        p.k_deg * Vi,
        'J_synth':      p.synth,
    }

def summarize(p: Params, label: str = ""):
    ss = steady_state(p)
    C, Va, Vi, A = ss
    active_nM  = (C + Va + A) * 1000
    inactive_nM = Vi * 1000
    total_nM   = (C + Va + Vi + A) * 1000
    header = f"  {label}" if label else ""
    print(f"{header:<25s}  Active = {active_nM:7.1f} nM  "
          f"Inactive = {inactive_nM:6.1f} nM  "
          f"Total = {total_nM:7.1f} nM  "
          f"|  Apo={A*1000:6.1f}  Cyt={C*1000:6.1f}  VacAct={Va*1000:6.1f}  VacInact={Vi*1000:6.1f}")
    return active_nM, inactive_nM, ss


# ---Main routine---

if __name__ == '__main__':
    p_wt = Params()

    print("=" * 70)
    print("ACTIVE / INACTIVE AUXIN POOL MODEL")
    print("=" * 70)
    print(f"  Vmax_LAX={p_wt.Vmax_LAX:.5f}  Vmax_PIN={p_wt.Vmax_PIN:.5f}  "
          f"Vmax_ABCB={p_wt.Vmax_ABCB_PM:.5f}  Vmax_WAT1={p_wt.Vmax_WAT1:.5f}")
    print(f"  k_inact={p_wt.k_inact:.4f}  k_GH3={p_wt.k_GH3:.4f}  "
          f"k_deg={p_wt.k_deg:.6f}  synth={p_wt.synth:.2e}")
    print()

    # ── WT ──
    active_wt, inactive_wt, ss_wt = summarize(p_wt, "Wild type")
    print(f"    target:  Active = {WT_ACTIVE:.1f} nM,  Inactive = {WT_INACTIVE:.1f} nM")

    # ── wat1 mutant ──
    p_wat1 = Params(Vmax_WAT1=0.0)
    active_mut, inactive_mut, _ = summarize(p_wat1, "wat1 mutant")
    print(f"    target:  Active = {WAT1_ACTIVE:.1f} nM,  Inactive = {WAT1_INACTIVE:.1f} nM")

    # ── error ──
    err_wt = abs(active_wt - WT_ACTIVE) + abs(inactive_wt - WT_INACTIVE)
    err_mut = abs(active_mut - WAT1_ACTIVE) + abs(inactive_mut - WAT1_INACTIVE)
    print(f"\n  Fit error: WT = {err_wt:.1f},  wat1 = {err_mut:.1f},  total = {err_wt+err_mut:.1f}")

    # ── redistribution time course ──
    print("\n" + "=" * 70)
    print("TIME COURSE — Redistribution (start: all active IAA in apoplast)")
    y0 = np.array([0.0, 0.0, 0.0, 0.105])  # all 105 nM active in apoplast
    res = simulate(p_wt, y0, t_end=7200, n=500)
    for t_sec in [0, 600, 1800, 3600, 7200]:
        i = np.argmin(np.abs(res['t'] - t_sec))
        act = (res['C'][i] + res['Va'][i] + res['A'][i]) * 1000
        inact = res['Vi'][i] * 1000
        print(f"  t={t_sec:5d}s  Apo={res['A'][i]*1000:6.1f}  "
              f"Cyt={res['C'][i]*1000:6.1f}  VacAct={res['Va'][i]*1000:6.1f}  "
              f"VacInact={res['Vi'][i]*1000:6.1f}  "
              f"Act={act:6.1f}  Inact={inact:6.1f}")

    # ── CSV export ──
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
    os.makedirs(out_dir, exist_ok=True)

    with open(os.path.join(out_dir, "redistribution_timecourse.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["t_s", "A_apo_nM", "A_cyt_nM", "A_vac_act_nM",
                     "A_vac_inact_nM", "active_nM", "inactive_nM", "total_nM"])
        for i in range(len(res['t'])):
            t_s = res['t'][i]
            apo, cyt, va, vi = res['A'][i]*1000, res['C'][i]*1000, res['Va'][i]*1000, res['Vi'][i]*1000
            act_nM = cyt + va + apo
            inact_nM = vi
            w.writerow([t_s, apo, cyt, va, vi, act_nM, inact_nM, act_nM + inact_nM])

    # Flux time series
    flux_names = ['J_LAX', 'J_PIN', 'J_ABCB_PM', 'J_pass_PM',
                  'J_ABCB_vac', 'J_pass_ton', 'J_inact', 'J_WAT1',
                  'J_GH3', 'J_deg', 'J_synth']
    with open(os.path.join(out_dir, "fluxes_timecourse.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["t_s"] + flux_names)
        for i in range(len(res['t'])):
            t_s = res['t'][i]
            f = fluxes_at_state(res['C'][i], res['Va'][i], res['Vi'][i], res['A'][i], p_wt)
            w.writerow([t_s] + [f[name] for name in flux_names])

    with open(os.path.join(out_dir, "model_fit.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["genotype", "active_nM", "inactive_nM", "total_nM",
                     "A_apo_nM", "A_cyt_nM", "A_vac_act_nM", "A_vac_inact_nM"])
        C, Va, Vi, A = ss_wt
        w.writerow(["WT", active_wt, inactive_wt, (C+Va+Vi+A)*1000,
                     A*1000, C*1000, Va*1000, Vi*1000])
        C2, Va2, Vi2, A2 = steady_state(p_wat1)
        w.writerow(["wat1", active_mut, inactive_mut, (C2+Va2+Vi2+A2)*1000,
                     A2*1000, C2*1000, Va2*1000, Vi2*1000])

    # Parameter table with units and references
    with open(os.path.join(out_dir, "parameters.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["parameter", "value", "unit", "source"])
        w.writerows([
            ["SCALE", SCALE, "uM/s per a.u.", "grid-search fit to WT/wat1 data"],
            ["Vmax_AUX1/LAX", round(p_wt.Vmax_LAX, 6), "uM/s", "SCALE * 2.97; Yang 2006, Plant J 47:521-536"],
            ["Vmax_PIN", round(p_wt.Vmax_PIN, 6), "uM/s", "SCALE * 3.77; Petrasek 2006, Science 312:914-918"],
            ["Vmax_ABCB1/19", round(p_wt.Vmax_ABCB_PM, 6), "uM/s", "SCALE * 2.77; Geisler 2005, Plant J 44:179-194"],
            ["Vmax_ABCB_vac", round(p_wt.Vmax_ABCB_vac, 6), "uM/s", "SCALE * 2.77; estimated"],
            ["Vmax_WAT1", round(p_wt.Vmax_WAT1, 6), "uM/s", "SCALE * 6.96; Ranocha 2013, Nature Comm 4:2348"],
            ["Km_AUX1/LAX", p_wt.Km_LAX, "uM", "Yang 2006, Plant J 47:521-536"],
            ["Km_PIN", p_wt.Km_PIN, "uM", "Petrasek 2006, Science 312:914-918"],
            ["Km_ABCB1/19", p_wt.Km_ABCB_PM, "uM", "Geisler 2005, Plant J 44:179-194"],
            ["Km_ABCB_vac", p_wt.Km_ABCB_vac, "uM", "estimated (ABCB family affinity)"],
            ["Km_WAT1", p_wt.Km_WAT1, "uM", "grid-search fit (inactive IAA conjugates); lit. Km 50-200 uM for active IAA"],
            ["P_IAH_PM", p_wt.P_IAH_PM, "1/s", "Gutknecht & Walter 1980, J Membr Biol 56:65-72"],
            ["P_IAH_tono", p_wt.P_IAH_tono, "1/s", "estimated (tonoplast less permeable than PM)"],
            ["k_inact", p_wt.k_inact, "1/s", "grid-search fit (active -> inactive in vacuole)"],
            ["k_GH3", p_wt.k_GH3, "1/s", "grid-search fit; lit. GH3 level = 0.52 a.u."],
            ["k_deg", p_wt.k_deg, "1/s", "grid-search fit (inactive degradation in vacuole)"],
            ["synth", p_wt.synth, "uM/s", "grid-search fit (basal synthesis)"],
            ["pH_apoplast", p_wt.pH_apo, "", "Grignon & Sentenac 1991"],
            ["pH_cytoplasm", p_wt.pH_cyt, "", "standard plant cell physiology"],
            ["pH_vacuole", p_wt.pH_vac, "", "standard plant cell physiology"],
            ["pKa(IAA)", pKa_IAA, "", "Rubery & Sheldrake 1974, Planta 118:101-121"],
            ["TOTAL_IAA_target", 105+10.8, "nM", "experimental measurement (WT active+inactive)"],
            ["PIN_expression", 3.77, "a.u.", "transciptome data"],
            ["AUX1-LAX_expression", 2.97, "a.u.", "transciptome data"],
            ["WAT1_expression", 6.96, "a.u.", "transciptome data"],
            ["ABCB1_expression", 2.77, "a.u.", "transciptome data"],
            ["GH3_expression", 0.52, "a.u.", "transciptome data"],
        ])

    print(f"\nCSVs written to {out_dir}/")

    # ── Transport mutant simulations ──
    mutants = {
        'wt':              {},
        'wat1':            {'Vmax_WAT1': 0.0},
        'wat1_quarter':    {'Vmax_WAT1': V_WAT1 * 0.25},
        'wat1_third':      {'Vmax_WAT1': V_WAT1 / 3.0},
        'wat1_half':       {'Vmax_WAT1': V_WAT1 * 0.5},
        'wat1_sixth':      {'Vmax_WAT1': V_WAT1 / 6.0},
        'wat1_5pct':       {'Vmax_WAT1': V_WAT1 * 0.05},
        'pin':             {'Vmax_PIN': 0.0},
        'aux1_lax':        {'Vmax_LAX': 0.0},
        'abcb_double':     {'Vmax_ABCB_PM': 0.0, 'Vmax_ABCB_vac': 0.0},
        'gh3':             {'k_GH3': 0.0},
    }

    print("\n" + "=" * 70)
    print("TRANSPORT MUTANT SIMULATIONS")
    print("=" * 70)

    for name, kw in mutants.items():
        p = Params(**kw)
        act, inact, ss = summarize(p, name)
        C, Va, Vi, A = ss

        # time course — start with active in apoplast matching WT total
        y0 = np.array([0.0, 0.0, 0.0, 0.105])
        res = simulate(p, y0, t_end=7200, n=300)

        # save concentrations
        fn = os.path.join(out_dir, f"{name}_timecourse.csv")
        with open(fn, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["t_s", "A_apo_nM", "A_cyt_nM", "A_vac_act_nM",
                         "A_vac_inact_nM", "active_nM", "inactive_nM", "total_nM"])
            for i in range(len(res['t'])):
                t_s = res['t'][i]
                apo, cyt, va, vi = res['A'][i]*1000, res['C'][i]*1000, res['Va'][i]*1000, res['Vi'][i]*1000
                act_nM = cyt + va + apo
                w.writerow([t_s, apo, cyt, va, vi, act_nM, vi, act_nM + vi])

        # save fluxes
        fn_f = os.path.join(out_dir, f"{name}_fluxes.csv")
        flux_names = ['J_LAX', 'J_PIN', 'J_ABCB_PM', 'J_pass_PM',
                      'J_ABCB_vac', 'J_pass_ton', 'J_inact', 'J_WAT1',
                      'J_GH3', 'J_deg', 'J_synth']
        with open(fn_f, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["t_s"] + flux_names)
            for i in range(len(res['t'])):
                t_s = res['t'][i]
                f = fluxes_at_state(res['C'][i], res['Va'][i], res['Vi'][i], res['A'][i], p)
                w.writerow([t_s] + [f[n] for n in flux_names])

        print(f"  -> {fn}  +  {fn_f}")

    print("\nDone.")
