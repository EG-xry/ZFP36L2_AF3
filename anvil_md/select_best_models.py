#!/usr/bin/env python3
"""
Select the BEST AlphaFold3 model for each ZFP36L2 variant across all runs.

Strategy: Use AF3's ranking_score (which combines ipTM and pTM) to select
the single best model across all runs × 5 models per run.

For MD simulations, using the best-ranked model is the standard approach:
- It gives the most confident protein-RNA interface geometry
- ipTM directly measures interface quality (protein↔RNA contact prediction)
- Higher ranking_score = better starting structure for MD
- This avoids biasing MD results with a poor starting conformation

The ranking_score in AF3 = 0.8*ipTM + 0.2*pTM (interface-weighted)

Outputs:
- Converts best CIF → PDB using OpenMM or gemmi
- Creates a summary CSV of selections
- Copies best PDBs to a clean directory for upload to Anvil
"""

import json
import os
import sys
import csv
import shutil
from pathlib import Path

# Configuration
DATA_DIR = Path(__file__).resolve().parent.parent / "Data"
OUTPUT_DIR = Path(__file__).resolve().parent / "best_models"
OUTPUT_DIR.mkdir(exist_ok=True)

VARIANTS = {
    "WT": "fold_zfp36l2_wt_nonamer_wt_2zn",
    "Y154H": "fold_zfp36l2_y154h_nonamer_wt_2zn",
    "R160W": "fold_zfp36l2_r160w_nonamer_wt_2zn",
    "R184C": "fold_zfp36l2_r184c_nonamer_wt_2zn",
    "G204D": "fold_zfp36l2_g204d_nonamer_wt_2zn",
    "C206F": "fold_zfp36l2_c206f_nonamer_wt_2zn",
    "S185N": "fold_zfp36l2_s185n_nonamer_wt_2zn",
    "T187S": "fold_zfp36l2_t187s_nonamer_wt_2zn",
    "D210E": "fold_zfp36l2_d210e_nonamer_wt_2zn",
    "C208A": "fold_zfp36l2_deleterious_c208a",
}


def cif_to_pdb(cif_path, pdb_path):
    """Convert mmCIF to PDB format using gemmi (preferred) or fallback."""
    try:
        import gemmi
        doc = gemmi.cif.read(str(cif_path))
        st = gemmi.make_structure_from_block(doc[0])
        st.write_pdb(str(pdb_path))
        return True
    except ImportError:
        pass

    # Fallback: try OpenMM's PDBxFile
    try:
        from openmm.app import PDBxFile, PDBFile
        pdbx = PDBxFile(str(cif_path))
        with open(pdb_path, "w") as f:
            PDBFile.writeFile(pdbx.topology, pdbx.positions, f)
        return True
    except Exception:
        pass

    # Last resort: just copy the CIF and note it needs conversion on Anvil
    print(f"  WARNING: Could not convert CIF→PDB. Copying CIF as-is.")
    shutil.copy2(cif_path, str(pdb_path).replace(".pdb", ".cif"))
    return False


def find_best_model(variant, base_name):
    """Find the best model across all runs for a variant."""
    var_dir = DATA_DIR / variant
    if not var_dir.is_dir():
        return None

    best_score = -1
    best_info = None

    for run_dir in sorted(var_dir.iterdir()):
        if not run_dir.is_dir():
            continue

        for model_idx in range(5):
            conf_file = run_dir / f"{base_name}_summary_confidences_{model_idx}.json"
            if not conf_file.exists():
                continue

            try:
                with open(conf_file) as f:
                    conf = json.load(f)

                ranking = conf.get("ranking_score", 0)
                iptm = conf.get("iptm", 0)
                ptm = conf.get("ptm", 0)

                # Get protein-RNA pair ipTM
                pair_iptm = None
                cpi = conf.get("chain_pair_iptm", [])
                if len(cpi) > 0 and len(cpi[0]) > 3:
                    pair_iptm = cpi[0][3]

                if ranking > best_score:
                    best_score = ranking
                    cif_path = run_dir / f"{base_name}_model_{model_idx}.cif"
                    best_info = {
                        "variant": variant,
                        "run": run_dir.name,
                        "model": model_idx,
                        "ranking_score": ranking,
                        "iptm": iptm,
                        "ptm": ptm,
                        "pair_iptm_prot_rna": pair_iptm,
                        "cif_path": str(cif_path),
                    }
            except Exception:
                pass

    return best_info


def main():
    print("=" * 70)
    print("Selecting BEST AF3 model per variant (by ranking_score)")
    print("ranking_score = 0.8*ipTM + 0.2*pTM (interface-weighted)")
    print("=" * 70)
    print()

    results = []

    for variant, base_name in sorted(VARIANTS.items()):
        print(f"[{variant}] Scanning all runs...")
        info = find_best_model(variant, base_name)

        if info is None:
            print(f"  ⚠ No data found for {variant}")
            continue

        print(f"  Best: {info['run']} model_{info['model']}")
        print(f"    ranking={info['ranking_score']:.4f}  ipTM={info['iptm']:.4f}  pTM={info['ptm']:.4f}")

        # Convert CIF to PDB
        cif_path = Path(info["cif_path"])
        pdb_name = f"{variant}_best_model.pdb"
        pdb_path = OUTPUT_DIR / pdb_name

        if cif_path.exists():
            success = cif_to_pdb(cif_path, pdb_path)
            if success:
                print(f"  ✓ Converted to PDB: {pdb_path}")
                info["pdb_output"] = str(pdb_path)
            else:
                print(f"  ⚠ CIF copied (needs conversion on Anvil)")
                info["pdb_output"] = str(pdb_path).replace(".pdb", ".cif")
        else:
            print(f"  ⚠ CIF file not found: {cif_path}")
            info["pdb_output"] = "MISSING"

        results.append(info)
        print()

    # Write summary CSV
    csv_path = OUTPUT_DIR / "best_model_selection.csv"
    if results:
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        print(f"\nSummary saved: {csv_path}")

    # Print comparison: run1_model0 vs best
    print("\n" + "=" * 70)
    print("COMPARISON: run1_model0 (currently uploaded) vs BEST model")
    print("=" * 70)
    print(f"{'Variant':<10} {'Run1 ranking':<15} {'Best ranking':<15} {'Best Run':<45} {'Δ':<8}")
    print("-" * 95)

    for info in results:
        variant = info["variant"]
        base_name = VARIANTS[variant]
        # Get run1 model0 score
        run1_dir = DATA_DIR / variant
        run1_score = None
        for d in sorted(run1_dir.iterdir()):
            if d.is_dir() and "(1)" in d.name:
                conf1 = d / f"{base_name}_summary_confidences_0.json"
                if conf1.exists():
                    with open(conf1) as f:
                        c = json.load(f)
                    run1_score = c.get("ranking_score", 0)
                break

        best = info["ranking_score"]
        if run1_score is not None:
            delta = best - run1_score
            print(f"{variant:<10} {run1_score:<15.4f} {best:<15.4f} {info['run']:<45} {delta:+.4f}")
        else:
            print(f"{variant:<10} {'N/A':<15} {best:<15.4f} {info['run']:<45}")

    print()
    print("RECOMMENDATION: Use the BEST models for MD simulation.")
    print("Upload command:")
    print(f'  scp "{OUTPUT_DIR}/"*.pdb YOUR_USER@anvil.rcac.purdue.edu:~/zfp36l2_md/inputs/')
    print()


if __name__ == "__main__":
    main()
