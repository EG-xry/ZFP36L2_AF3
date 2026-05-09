import os
import sys
import numpy as np
import matplotlib.pyplot as plt

try:
    import mdtraj as md
except Exception:
    md = None


def gaussian(x, mu, sigma):
    s2 = sigma * sigma
    return np.exp(-0.5 * (x - mu) * (x - mu) / s2)


def main():
    prep_dir = os.path.dirname(os.path.abspath(__file__))
    prob_path = os.path.join(prep_dir, "prob.npy")
    pdb_path = sys.argv[1] if len(sys.argv) > 1 else \
        "/Users/eric/Desktop/Joshua_Project/AF-P10636-3-F1-model_v6.pdb"
    i = int(sys.argv[2]) if len(sys.argv) > 2 else 50
    j = int(sys.argv[3]) if len(sys.argv) > 3 else 150
    sigma = float(sys.argv[4]) if len(sys.argv) > 4 else 2.0  # Å noise

    if not os.path.exists(prob_path):
        print("prob.npy not found at:", prob_path)
        sys.exit(1)

    prob = np.load(prob_path)
    if prob.ndim != 3 or i >= prob.shape[0] or j >= prob.shape[1]:
        print("Unexpected prob.npy shape or indices out of range:", prob.shape)
        sys.exit(1)

    p_prior = prob[i, j].astype(float)
    K = p_prior.shape[0]
    # approximate bin centers in [0, 32) Å matching AF-distance map convention
    bin_centers = np.linspace(0.0, 32.0, K, endpoint=False)

    d_obs = None
    if md is not None and os.path.exists(pdb_path):
        try:
            t = md.load(pdb_path)
            ca = t.topology.select("name CA")
            if len(ca) > max(i, j):
                d_obs = md.compute_distances(t, [(ca[i], ca[j])]).ravel()[0] * 10.0
        except Exception as e:
            print("Could not compute PDB distance:", repr(e))

    if d_obs is None:
        print("No PDB distance available; using a synthetic observation at the prior mean.")
        d_obs = np.sum(p_prior * bin_centers) / np.sum(p_prior)

    # Likelihood over bins for observed distance with Gaussian noise
    likelihood = gaussian(bin_centers, d_obs, sigma)
    # Posterior ∝ prior × likelihood
    p_post = p_prior * likelihood
    if p_post.sum() > 0:
        p_post /= p_post.sum()

    fig, ax = plt.subplots(1, 1, figsize=(7, 4))
    ax.plot(bin_centers, p_prior, lw=2, color="#4C72B0", label="Prior")
    ax.plot(bin_centers, p_post, lw=2, color="#55A868", label="Posterior")
    ax.axvline(d_obs, color="#DD8452", lw=2, ls="--", label=f"Obs (PDB) ≈ {d_obs:.2f} Å")
    ax.set_xlabel("Distance (Å)")
    ax.set_ylabel("Probability")
    ax.set_title(f"Bayesian update for pair (i={i}, j={j}), σ={sigma:.1f} Å")
    ax.legend()
    plt.tight_layout()
    out_path = os.path.join(prep_dir, "bayesian_posterior_demo.png")
    plt.savefig(out_path, dpi=200)
    print("Saved:", out_path)


if __name__ == "__main__":
    main()


