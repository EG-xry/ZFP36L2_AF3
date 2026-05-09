#!/usr/bin/env bash
set -euo pipefail

echo "== Binary checks =="
for cmd in python3 docker mpirun plumed gmx; do
  if command -v "$cmd" >/dev/null 2>&1; then
    echo "[OK] $cmd -> $(command -v "$cmd")"
  else
    echo "[MISSING] $cmd"
  fi
done

check_python_modules() {
  local py_bin="$1"
  local label="$2"
  if [ ! -x "$py_bin" ]; then
    echo "[SKIP] $label ($py_bin not found)"
    return
  fi
  echo
  echo "== Python module checks ($label) =="
  "$py_bin" - <<'PY'
import importlib.util
mods = ["openmm", "openmmplumed", "mpi4py", "mdtraj", "pandas", "numpy", "Bio"]
for m in mods:
    ok = importlib.util.find_spec(m) is not None
    print(f"{'[OK]' if ok else '[MISSING]'} {m}")
PY
}

check_python_modules "$(command -v python3)" "python3"
check_python_modules "/Users/eric/micromamba/envs/omm/bin/python" "micromamba:omm"
