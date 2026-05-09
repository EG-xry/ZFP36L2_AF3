# Simulation Tuning Notes (AlphaFold-IDP)

## Key knobs to address NaN / ∞ energy failures

| Area | Parameters to adjust | Notes |
| --- | --- | --- |
| **Integrator** | `LangevinMiddleIntegrator(temp, friction, dt)` | We already dropped `dt` from 5 fs to 2 fs. If instabilities persist, try 1 fs or raise the friction constant (>0.01 ps⁻¹). |
| **Minimisation** | `minimizeEnergy(tol, maxIterations)` in `simulate_em.py` | Currently `tol=1 kJ mol⁻¹ nm⁻¹`, `maxIterations=1000`. Increase iterations or decrease tolerance for tougher systems. |
| **Velocity initialisation** | `context.setVelocitiesToTemperature` | Keep this after loading checkpoints so each replica starts with Maxwell-Boltzmann velocities. |
| **PLUMED PBMetaD** | `HEIGHT`, `PACE`, `SIGMA_MIN/MAX`, `BIASFACTOR` | Large `HEIGHT` or tight `SIGMA` can create very stiff bias forces. Consider lowering `HEIGHT` (e.g. 0.3) or increasing `PACE` to deposit hills less frequently. |
| **Metainference noise model** | `SIGMA0`, `SIGMA_MIN/MAX`, `MC_STEPS`, `MC_CHUNKSIZE` | If the Bayesian noise shrinks too aggressively, forces can spike. Widening `SIGMA_MIN` or lowering `MC_STEPS` slows the adaptation. |
| **Checkpointing / restart frequency** | `StateDataReporter` stride & manual checkpoints | More frequent checkpoints make it easier to restart from before an instability. |

## Debug workflow

1. **Energy drift check** – watch `stats_*.csv` to ensure potential energy/temperature plateau near physical values.
2. **Short validation run** – run a 2–5 minute `timeout` test after each parameter tweak (as done with `short_run.log`).
3. **Coefficient sweeps** – change only one knob at a time (timestep, hill height, etc.) so the effect is attributable.
4. **Archive outputs** – copy `output_*.dcd`, `COLVAR.*`, `BAYES*`, and `stats_*.csv` before reruns; they are overwritten.

> This file is informational only; no simulation parameters were changed here. Refer to `simulate.py`, `simulate_em.py`, and `plumed.dat` when applying any of the suggestions above.

