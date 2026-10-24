#!/usr/bin/env python3
"""TDD tests for the aromatic pi-stacking geometry core (unit-agnostic; Angstrom in use).

Run:  /Users/eric/micromamba/envs/omm/bin/python expanse_md_analysis/test_stacking_geometry.py
(Self-contained runner; no pytest dependency.)
"""
import numpy as np
from stacking_geometry import (
    ring_centroid, ring_normal, interplane_angle_deg, is_stacked,
    ring_centroids_frames, ring_normals_frames, stacking_occupancy, stacking_series,
)


def hexagon(center, e1, e2, radius=1.4):
    """Regular hexagon (6 atoms) centred at `center`, spanning unit vectors e1,e2."""
    center = np.asarray(center, float); e1 = np.asarray(e1, float); e2 = np.asarray(e2, float)
    ang = np.deg2rad(np.arange(0, 360, 60))
    return center + radius * (np.outer(np.cos(ang), e1) + np.outer(np.sin(ang), e2))


def approx(a, b, tol=1e-6):
    return abs(float(a) - float(b)) < tol


def test_ring_centroid_is_geometric_center():
    ring = hexagon([1, 2, 3], [1, 0, 0], [0, 1, 0])
    assert np.allclose(ring_centroid(ring), [1, 2, 3], atol=1e-9)


def test_ring_normal_of_xy_plane_ring_is_z_axis():
    ring = hexagon([0, 0, 0], [1, 0, 0], [0, 1, 0])
    n = ring_normal(ring)
    assert approx(abs(n[2]), 1.0)          # normal parallel to z (sign arbitrary)


def test_interplane_angle_parallel_is_zero():
    assert approx(interplane_angle_deg([0, 0, 1], [0, 0, 1]), 0.0)


def test_interplane_angle_antiparallel_is_zero():
    # opposite normals still describe PARALLEL planes -> stacking angle 0
    assert approx(interplane_angle_deg([0, 0, 1], [0, 0, -1]), 0.0)


def test_interplane_angle_perpendicular_is_ninety():
    assert approx(interplane_angle_deg([0, 0, 1], [1, 0, 0]), 90.0)


def test_stacked_true_for_close_parallel_rings():
    p = hexagon([0, 0, 0], [1, 0, 0], [0, 1, 0])
    b = hexagon([0, 0, 3.5], [1, 0, 0], [0, 1, 0])
    assert is_stacked(p, b, d_max=5.0, angle_max_deg=40.0) is True


def test_not_stacked_when_too_far():
    p = hexagon([0, 0, 0], [1, 0, 0], [0, 1, 0])
    b = hexagon([0, 0, 8.0], [1, 0, 0], [0, 1, 0])
    assert is_stacked(p, b, d_max=5.0, angle_max_deg=40.0) is False


def test_not_stacked_when_perpendicular():
    p = hexagon([0, 0, 0], [1, 0, 0], [0, 1, 0])      # normal ~ z
    b = hexagon([0, 0, 3.5], [0, 0, 1], [0, 1, 0])    # normal ~ x -> planes perpendicular
    assert is_stacked(p, b, d_max=5.0, angle_max_deg=40.0) is False


def test_batched_normals_match_per_frame():
    frames = np.stack([hexagon([0, 0, 0], [1, 0, 0], [0, 1, 0]),
                       hexagon([0, 0, 0], [1, 0, 0], [0, 0, 1])])  # (2,6,3)
    n = ring_normals_frames(frames)
    assert approx(abs(n[0][2]), 1.0)   # frame0 normal ~ z
    assert approx(abs(n[1][1]), 1.0)   # frame1 normal ~ y


def test_batched_centroids():
    frames = np.stack([hexagon([1, 0, 0], [1, 0, 0], [0, 1, 0]),
                       hexagon([0, 2, 0], [1, 0, 0], [0, 1, 0])])
    c = ring_centroids_frames(frames)
    assert np.allclose(c, [[1, 0, 0], [0, 2, 0]], atol=1e-9)


def test_stacking_occupancy_counts_any_base_over_frames():
    # protein ring parallel in both frames at origin
    prot = np.stack([hexagon([0, 0, 0], [1, 0, 0], [0, 1, 0]),
                     hexagon([0, 0, 0], [1, 0, 0], [0, 1, 0])])
    # base1: always far (8A) -> never stacks
    base1 = np.stack([hexagon([0, 0, 8.0], [1, 0, 0], [0, 1, 0]),
                      hexagon([0, 0, 8.0], [1, 0, 0], [0, 1, 0])])
    # base2: frame0 stacked (3.5A parallel); frame1 perpendicular -> not stacked
    base2 = np.stack([hexagon([0, 0, 3.5], [1, 0, 0], [0, 1, 0]),
                      hexagon([0, 0, 3.5], [0, 0, 1], [0, 1, 0])])
    occ = stacking_occupancy(prot, [base1, base2], d_max=5.0, angle_max_deg=40.0)
    assert approx(occ, 0.5)


def test_stacking_series_returns_per_frame_bool():
    prot = np.stack([hexagon([0, 0, 0], [1, 0, 0], [0, 1, 0]),
                     hexagon([0, 0, 0], [1, 0, 0], [0, 1, 0])])
    base2 = np.stack([hexagon([0, 0, 3.5], [1, 0, 0], [0, 1, 0]),      # frame0 stacked
                      hexagon([0, 0, 3.5], [0, 0, 1], [0, 1, 0])])     # frame1 perpendicular
    s = stacking_series(prot, [base2], d_max=5.0, angle_max_deg=40.0)
    assert list(np.asarray(s, bool)) == [True, False]


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for t in tests:
        try:
            t(); print(f"PASS  {t.__name__}")
        except Exception as e:
            failed += 1; print(f"FAIL  {t.__name__}: {type(e).__name__}: {e}")
    print(f"\n{len(tests)-failed}/{len(tests)} passed")
    raise SystemExit(1 if failed else 0)
