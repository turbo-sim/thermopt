"""Snapshot the best matching saved coarse-grid state; never run optimization."""

from copy import deepcopy
import hashlib
import json
import shutil

import run_fluid_comparison as run

SOURCE_VERSIONS = ("rankine_wet2dry_v3", "rankine_wet2dry_v2")
RUN_FAMILIES = ("butane", "cyclopentane", "novec")
DESTINATION = run.HERE / "seed_results"


def schema(config):
    config = deepcopy(config)
    for variable in config["design_variables"].values():
        variable.pop("value", None)
    return json.dumps(config, sort_keys=True, default=lambda value: value.tolist())


def prepare_family(family):
    output = DESTINATION / f"sensitivity_{family}"
    if (output / "manifest.csv").exists():
        return
    original = run.HERE.parent / SOURCE_VERSIONS[0] / "results" / f"sensitivity_{family}"
    table = run.pd.read_csv(original / "sweep_results.csv", float_precision="round_trip")
    output.mkdir(parents=True, exist_ok=True)
    manifest = []
    for row in table.itertuples():
        if row.source_min_C == 140:
            continue
        candidates = []
        expected = schema(run.th.read_configuration_file(original / row.case / "optimal_solution.yaml"))
        for version in SOURCE_VERSIONS:
            folder = run.HERE.parent / version / "results" / f"sensitivity_{family}" / row.case
            if not (folder / "result.json").exists():
                continue
            saved = json.loads((folder / "result.json").read_text(encoding="utf-8"))
            if not saved["row"]["success"]:
                continue
            config = run.th.read_configuration_file(folder / "optimal_solution.yaml")
            if schema(config) != expected:
                continue
            candidates.append((saved["row"]["system_efficiency_pct"], version, folder, config, saved))
        if not candidates:
            raise RuntimeError(f"No matching saved optimum for {family}/{row.case}")
        score, version, folder, config, saved = max(candidates, key=lambda candidate: candidate[0])
        assert all(config["fixed_parameters"][name]["num_elements"] == 25 for name in ("heater", "cooler"))
        problem = run.th.ThermodynamicCycleProblem(config, out_dir=str(output))
        problem.fitness(problem.x0)
        assert run.physical_cycle_valid(problem) and run.violation(problem) <= 1e-5
        assert abs(100 * float(problem.cycle_data["energy_analysis"]["system_efficiency"]) - score) < 1e-7
        destination = output / row.case
        destination.mkdir(exist_ok=True)
        shutil.copy2(folder / "optimal_solution.yaml", destination / "optimal_solution.yaml")
        shutil.copy2(folder / "result.json", destination / "coarse_result.json")
        source = dict(case=row.case, source_version=version,
                      source_path=str(folder.relative_to(run.HERE.parent)),
                      coarse_system_efficiency_pct=float(score),
                      coarse_fixed_reference_efficiency_pct=(100000 * saved["row"]["specific_system_power_kW_per_kg_s"]
                                                            / run.fixed_reference_heat_J_per_kg()),
                      seed_sha256=hashlib.sha256((destination / "optimal_solution.yaml").read_bytes()).hexdigest())
        (destination / "source.json").write_text(json.dumps(source, indent=2), encoding="utf-8")
        manifest.append(source)
    assert len(manifest) == 120
    run.pd.DataFrame(manifest).to_csv(output / "manifest.csv", index=False)
    print(f"Prepared {family}: 120 exact saved coarse-grid seeds.", flush=True)


if __name__ == "__main__":
    for family in RUN_FAMILIES:
        prepare_family(family)
