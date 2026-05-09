import os
import sys
import numpy as np
import matplotlib.pyplot as plt

try:
    import mdtraj as md
except Exception:
    md = None


def main():
    # Inputs
    prep_dir = os.path.dirname(os.path.abspath(__file__))
    prob_path = os.path.join(prep_dir, "prob.npy")
    # User PDB path passed as arg or default
    pdb_path = sys.argv[1] if len(sys.argv) > 1 else \
        "/Users/eric/Desktop/Joshua_Project/AF-P10636-3-F1-model_v6.pdb"
    # Residue pair (0-based indices on CA)
    i = int(sys.argv[2]) if len(sys.argv) > 2 else 50
    j = int(sys.argv[3]) if len(sys.argv) > 3 else 150

    if not os.path.exists(prob_path):
        print("prob.npy not found at:", prob_path)
        sys.exit(1)

    prob = np.load(prob_path)
    if prob.ndim != 3 or i >= prob.shape[0] or j >= prob.shape[1]:
        print("Unexpected prob.npy shape or indices out of range:", prob.shape)
        sys.exit(1)

    p_ij = prob[i, j]
    K = p_ij.shape[0]
    dist_bins = np.linspace(0.0, 32.0, K, endpoint=False)

    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    ax.plot(dist_bins, p_ij, lw=2, color="#4C72B0", label="AF prior")

    # Overlay PDB CA distance if possible
    if md is not None and os.path.exists(pdb_path):
        try:
            t = md.load(pdb_path)
            ca = t.topology.select("name CA")
            if len(ca) > max(i, j):
                d = md.compute_distances(t, [(ca[i], ca[j])]).ravel()[0] * 10.0
                ax.axvline(d, color="#DD8452", lw=2, ls="--", label=f"PDB CA ≈ {d:.2f} Å")
        except Exception as e:
            print("Could not overlay PDB distance:", repr(e))

    ax.set_xlabel("Distance (Å)")
    ax.set_ylabel("Probability")
    ax.set_title(f"P(d | i={i}, j={j}) prior + PDB overlay")
    ax.legend()
    plt.tight_layout()
    out_path = os.path.join(prep_dir, "bayesian_demo.png")
    plt.savefig(out_path, dpi=200)
    print("Saved:", out_path)


if __name__ == "__main__":
    main()


