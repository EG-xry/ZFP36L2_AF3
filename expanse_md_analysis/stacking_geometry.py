#!/usr/bin/env python3
"""Aromatic pi-stacking geometry core (unit-agnostic; used in Angstrom).

Stacking = two ring planes that are both CLOSE (centroid distance < d_max) and
near-PARALLEL (inter-plane angle < angle_max). Plane normal via best-fit (SVD);
angle uses |cos| so antiparallel normals (parallel planes) read as 0 deg.
"""
from __future__ import annotations
import numpy as np


def ring_centroid(xyz) -> np.ndarray:
    return np.asarray(xyz, float).mean(axis=0)


def ring_normal(xyz) -> np.ndarray:
    xyz = np.asarray(xyz, float)
    _, _, vt = np.linalg.svd(xyz - xyz.mean(axis=0), full_matrices=False)
    n = vt[-1]
    return n / np.linalg.norm(n)


def interplane_angle_deg(n1, n2) -> float:
    n1 = np.asarray(n1, float); n2 = np.asarray(n2, float)
    n1 = n1 / np.linalg.norm(n1); n2 = n2 / np.linalg.norm(n2)
    c = min(1.0, max(0.0, abs(float(np.dot(n1, n2)))))
    return float(np.degrees(np.arccos(c)))


def is_stacked(ring_p, ring_b, d_max: float = 5.0, angle_max_deg: float = 40.0) -> bool:
    d = float(np.linalg.norm(ring_centroid(ring_p) - ring_centroid(ring_b)))
    ang = interplane_angle_deg(ring_normal(ring_p), ring_normal(ring_b))
    return bool(d < d_max and ang < angle_max_deg)


# ---- batched-over-frames versions (X: (F, N_atoms, 3)) ----

def ring_centroids_frames(X) -> np.ndarray:
    return np.asarray(X, float).mean(axis=1)


def ring_normals_frames(X) -> np.ndarray:
    X = np.asarray(X, float)
    c = X.mean(axis=1, keepdims=True)
    _, _, vt = np.linalg.svd(X - c, full_matrices=False)   # batched over frames
    n = vt[:, -1, :]
    return n / np.linalg.norm(n, axis=1, keepdims=True)


def stacking_series(prot_frames, base_frames_list,
                    d_max: float = 5.0, angle_max_deg: float = 40.0) -> np.ndarray:
    """Per-frame boolean: does the protein ring stack with AT LEAST ONE base ring?

    prot_frames: (F, Np, 3);  base_frames_list: list of (F, Nb, 3).
    """
    cp = ring_centroids_frames(prot_frames)      # (F,3)
    npn = ring_normals_frames(prot_frames)       # (F,3)
    F = cp.shape[0]
    stacked_any = np.zeros(F, dtype=bool)
    for b in base_frames_list:
        cb = ring_centroids_frames(b); nb = ring_normals_frames(b)
        d = np.linalg.norm(cp - cb, axis=1)
        ang = np.degrees(np.arccos(np.clip(np.abs(np.sum(npn * nb, axis=1)), 0.0, 1.0)))
        stacked_any |= (d < d_max) & (ang < angle_max_deg)
    return stacked_any


def stacking_occupancy(prot_frames, base_frames_list,
                       d_max: float = 5.0, angle_max_deg: float = 40.0) -> float:
    """Fraction of frames where the protein ring stacks with at least one base ring."""
    return float(np.mean(stacking_series(prot_frames, base_frames_list,
                                         d_max=d_max, angle_max_deg=angle_max_deg)))
