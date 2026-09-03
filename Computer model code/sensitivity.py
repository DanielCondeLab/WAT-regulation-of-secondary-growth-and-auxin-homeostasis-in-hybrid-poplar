"""
---Sobol Global Sensitivity Analysis — Active / Inactive Auxin Model---
Outcomes:
  Y1 = WT_active   — total active auxin in wild type (nM)
  Y2 = WT_inactive — inactive auxin in wild type (nM)
  Y3 = WT_total    — total auxin (active+inactive) in WT (nM)
  Y4 = M_ratio     — wat1_active / WT_active  (drop ratio)
"""
import numpy as np
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol as sobol_analyze
from auxin_transport_model import Params, steady_state
import csv, os, warnings
warnings.filterwarnings('ignore')

problem = {
    'num_vars': 14,
    'names': [
        'Vmax_LAX', 'Vmax_PIN', 'Vmax_ABCB', 'Vmax_WAT1',
        'Km_LAX', 'Km_PIN', 'Km_ABCB_PM', 'Km_ABCB_vac', 'Km_WAT1',
        'P_IAH_PM',
        'k_inact', 'k_GH3', 'k_deg', 'synth',
    ],
    'bounds': [
        [0.0005,  0.008],    # Vmax_LAX
        [0.0005,  0.008],    # Vmax_PIN
        [0.0005,  0.008],    # Vmax_ABCB
        [0.001,   0.02],     # Vmax_WAT1
        [0.3,     3.0],      # Km_LAX
        [3.0,    30.0],      # Km_PIN
        [0.5,     5.0],      # Km_ABCB_PM
        [0.3,     3.0],      # Km_ABCB_vac
        [0.1,    10.0],      # Km_WAT1
        [1e-4,    5e-3],     # P_IAH_PM
        [0.001,   0.02],     # k_inact
        [0.00001, 0.001],    # k_GH3
        [0.0001,  0.005],    # k_deg
        [1e-6,    5e-5],     # synth
    ],
}

def evaluate(X_row):
    (Vmax_LAX, Vmax_PIN, Vmax_ABCB, Vmax_WAT1,
     Km_LAX, Km_PIN, Km_ABCB_PM, Km_ABCB_vac, Km_WAT1,
     P_IAH_PM, k_inact, k_GH3, k_deg, synth) = X_row

    try:
        p = Params(
            Vmax_LAX=Vmax_LAX, Vmax_PIN=Vmax_PIN,
            Vmax_ABCB_PM=Vmax_ABCB, Vmax_ABCB_vac=Vmax_ABCB,
            Vmax_WAT1=Vmax_WAT1,
            Km_LAX=Km_LAX, Km_PIN=Km_PIN,
            Km_ABCB_PM=Km_ABCB_PM, Km_ABCB_vac=Km_ABCB_vac,
            Km_WAT1=Km_WAT1, P_IAH_PM=P_IAH_PM,
            k_inact=k_inact, k_GH3=k_GH3, k_deg=k_deg, synth=synth)
        ss = steady_state(p)
        C, Va, Vi, A = ss
        wt_act  = max(0, (C + Va + A) * 1000)
        wt_inact = max(0, Vi * 1000)
        wt_tot  = wt_act + wt_inact

        p2 = Params(
            Vmax_LAX=Vmax_LAX, Vmax_PIN=Vmax_PIN,
            Vmax_ABCB_PM=Vmax_ABCB, Vmax_ABCB_vac=Vmax_ABCB,
            Vmax_WAT1=0.0, Km_LAX=Km_LAX, Km_PIN=Km_PIN,
            Km_ABCB_PM=Km_ABCB_PM, Km_ABCB_vac=Km_ABCB_vac,
            Km_WAT1=Km_WAT1, P_IAH_PM=P_IAH_PM,
            k_inact=k_inact, k_GH3=k_GH3, k_deg=k_deg, synth=synth)
        ss2 = steady_state(p2)
        C2, Va2, Vi2, A2 = ss2
        mut_act = max(0, (C2 + Va2 + A2) * 1000)
        m_ratio = mut_act / wt_act if wt_act > 0 else np.nan

        return (wt_act, wt_inact, wt_tot, m_ratio)
    except Exception:
        return (np.nan,) * 4

if __name__ == '__main__':
    N = 1024
    print(f"Generating {N} base samples x {problem['num_vars']} params = {N*(problem['num_vars']+2)} runs ...")
    param_values = sobol_sample.sample(problem, N, calc_second_order=False)

    print("Evaluating model ...")
    Y = np.zeros((param_values.shape[0], 4))
    for i in range(param_values.shape[0]):
        Y[i, :] = evaluate(param_values[i, :])
        if (i + 1) % 5000 == 0:
            print(f"  {i+1}/{param_values.shape[0]} done")

    valid = ~np.isnan(Y).any(axis=1)
    Y_clean = Y[valid]
    param_clean = param_values[valid]
    n_dropped = Y.shape[0] - Y_clean.shape[0]
    print(f"  Dropped {n_dropped} failed runs ({n_dropped/Y.shape[0]*100:.1f}%)")

    output_names = ['WT_active (nM)', 'WT_inactive (nM)', 'WT_total (nM)', 'M_ratio (wat1/WT)']

    print("\n" + "=" * 80)
    print("SOBOL FIRST-ORDER (S1) AND TOTAL-ORDER (ST) SENSITIVITY INDICES")
    print("=" * 80)

    Si_cache = []
    for j, oname in enumerate(output_names):
        print(f"\n--- {oname} ---")
        Si = sobol_analyze.analyze(problem, Y_clean[:, j], calc_second_order=False, print_to_console=False)
        Si_cache.append(Si)

        idx = np.argsort(Si['ST'])[::-1]
        print(f"  {'Parameter':<20s}  {'S1':>8s}  {'ST':>8s}")
        print(f"  {'-'*20}  {'-'*8}  {'-'*8}")
        for i in idx[:10]:
            print(f"  {problem['names'][i]:<20s}  {max(0, Si['S1'][i]):8.4f}  {max(0, Si['ST'][i]):8.4f}")
        print(f"\n  Top:")
        for rk, i in enumerate(idx[:5]):
            print(f"    {rk+1}. {problem['names'][i]:<20s} ST={Si['ST'][i]:.4f}")

    # Summary matrix
    print("\n" + "=" * 80)
    print("ST SUMMARY MATRIX")
    print("=" * 80)
    print(f"  {'Parameter':<20s}", end="")
    for on in output_names:
        print(f"  {on[:14]:>14s}", end="")
    print()
    for pi, pn in enumerate(problem['names']):
        print(f"  {pn:<20s}", end="")
        for j in range(4):
            print(f"  {max(0, Si_cache[j]['ST'][pi]):14.4f}", end="")
        print()

    # CSV
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "sobol_indices.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["parameter", "output", "S1", "ST"])
        for j, on in enumerate(output_names):
            for i, pn in enumerate(problem['names']):
                w.writerow([pn, on, max(0, Si_cache[j]['S1'][i]), max(0, Si_cache[j]['ST'][i])])

    with open(os.path.join(out_dir, "sobol_st_matrix.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["parameter"] + output_names)
        for i, pn in enumerate(problem['names']):
            w.writerow([pn] + [f"{max(0, Si_cache[j]['ST'][i]):.4f}" for j in range(4)])

    print("\nDone.")
