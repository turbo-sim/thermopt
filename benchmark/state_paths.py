"""Time complete source PTES cycles, including nondimensional state paths.

This complements run.py's source/bicubic profiles: PTES exercises the radial
machinery branch whose public return changes from lists to native batched states.
Shared workload definitions currently come from the regression helpers, pending
Phase 7's consolidation. No thermodynamic equations are duplicated here.
"""

import argparse
import hashlib
import importlib.metadata as metadata
import os
from pathlib import Path
import platform
import statistics
import sys
from time import perf_counter

from run import ROOT, cpu_affinity, save, sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=("ptes", "ptes_turbo"), default="ptes_turbo")
    parser.add_argument("--samples", type=int, default=12)
    parser.add_argument("--cpu-affinity", type=lambda value: int(value, 0))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.samples < 4:
        parser.error("Use at least four complete-cycle samples")
    if args.output.exists():
        parser.error("Existing results are never overwritten")
    affinity = cpu_affinity(args.cpu_affinity)
    os.environ.update(
        JAX_ENABLE_X64="true",
        JAX_PLATFORM_NAME="cpu",
        MPLBACKEND="Agg",
        OMP_NUM_THREADS="1",
    )
    import jax
    import jaxprop as jp

    if not Path(jp.__file__).resolve().is_relative_to(Path(sys.prefix).resolve()):
        raise RuntimeError("Use installed JaxProp, without checkout shadowing")
    if not jax.config.x64_enabled or jax.default_backend() != "cpu":
        raise RuntimeError("This benchmark requires CPU/x64 execution")
    sys.path.insert(0, str(ROOT / "tests"))
    from regression_support import make_problem, snapshot

    def source_hash():
        return hashlib.sha256(
            b"".join(p.read_bytes() for p in sorted((ROOT / "thermopt").rglob("*.py")))
        ).hexdigest()

    before = source_hash()
    script_hash = sha(__file__)
    input_hash = sha(ROOT / f"tests/data/cases/{args.case}.yaml")
    helper_hash = sha(ROOT / "tests/regression_support.py")
    lock_hash = sha(ROOT / "poetry.lock")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    start = perf_counter()
    problem = make_problem(
        args.case, args.output.parent / (args.output.stem + "_artifacts")
    )
    setup = perf_counter() - start
    initial = snapshot(problem, problem.x0)
    observations = []
    for _ in range(args.samples):
        start = perf_counter()
        problem.fitness(problem.x0)
        jax.block_until_ready(problem.cycle_data)
        observations.append(perf_counter() - start)
    assert snapshot(problem, problem.x0) == initial
    assert source_hash() == before and sha(__file__) == script_hash
    assert sha(ROOT / f"tests/data/cases/{args.case}.yaml") == input_hash
    assert sha(ROOT / "tests/regression_support.py") == helper_hash
    assert sha(ROOT / "poetry.lock") == lock_hash
    assert cpu_affinity()["actual_mask"] == affinity["actual_mask"]
    save(
        args.output,
        dict(
            case=args.case,
            samples=args.samples,
            execution="complete eager source cycle; all outputs retained and synchronized",
            environment=dict(
                python=platform.python_version(),
                platform=platform.platform(),
                processor=platform.processor(),
                logical_cpus=os.cpu_count(),
                cpu_affinity=affinity,
                x64=jax.config.x64_enabled,
                dependencies={
                    name: metadata.version(name)
                    for name in (
                        "jaxprop",
                        "jax",
                        "jaxlib",
                        "CoolProp",
                        "numpy",
                        "scipy",
                        "pysolver_view",
                    )
                },
            ),
            thermopt_root=str(ROOT),
            jaxprop_path=str(Path(jp.__file__).resolve()),
            source_sha256=before,
            harness_sha256=script_hash,
            helper_sha256=helper_hash,
            input_sha256=input_hash,
            lock_sha256=lock_hash,
            setup_seconds=setup,
            warm_samples_seconds=observations,
            warm_median_seconds=statistics.median(observations),
            numerical_snapshot=initial,
            repeated_numerics_identical=True,
        ),
    )
    print(args.output)


if __name__ == "__main__":
    main()
