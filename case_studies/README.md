# ThermOpt case studies

This directory contains the cycle cases and research programs formerly under
`examples/`. Each case retains its original subdirectory layout, YAML inputs,
scripts, reference material and existing local results. The existing
`partial_evaporation_orc_offdesign_prompt.md` is also preserved here.

Available case families include:

- `sCO2_recuperated/`, `sCO2_recompression/`, `sCO2_comparison/`.
- `rankine_steam/`, `rankine_organic/`, `rankine_wet2dry/`,
  `rankine_transcritical_vs_trilateral/`.
- `heat_pump_CO2/`, `refrigeration_air/`, `refrigeration_R134a/`.
- `PTES_CO2/`, `demo_hx_node_snap/`, `HX_development/`.

Read the individual case instructions before running it. Scientific sweeps,
historical figure preparation and input-dependent programs have different
requirements; relocating them does not establish that every full study runs.
Keep generated outputs, tables and logs inside the case's `results/` directory.

For the verified six-way recompression solver/backend comparison, run from the
repository root after installing the locked dependencies:

```powershell
.\.venv\Scripts\python.exe case_studies/sCO2_recompression/compare_solvers.py --no-show
```

This compares controlled SLSQP, Pounce interior point and Pounce SQP using
CoolProp HEOS and a runner-defined 128x128 native bicubic table. It preserves
the input YAML and saves figures, histories and an ASCII timing comparison in
`sCO2_recompression/results/solver_comparison/`.

Most `run_optimization.py` launchers expect their case directory as the working
directory. Standalone comparison runners resolve their own input paths; their
case README documents the supported invocation.

Tests keep their frozen physical inputs under `tests/data/cases/`. Historical
benchmark JSON, checkpoints and audits retain the original `examples/` paths
and numerical hashes. Those records describe the run that produced them; active
source code and current instructions use `case_studies/`.

The relocation was verified on 16 November 2026: all 660 tracked files were
recognized as renames, all original local files were retained, and physical
YAML/data/figure bytes were unchanged. The full suite passed 543 tests with four
unchanged strict expected failures. The maintained recuperated source/bicubic
smoke profile also passed its numerical checks from the new layout; this is
path/workflow verification, not a new between-phase performance acceptance.
