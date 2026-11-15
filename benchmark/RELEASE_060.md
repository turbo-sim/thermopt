# Published JaxProp 0.6 adoption checkpoint

Official PyPI metadata confirms 0.6.0 was uploaded on 15 November 2026 at
13:14 UTC. Wheel SHA256:
`0a8dfb68ca4d21cf97b7483e1da776b946cec108a94d50ebb770bb9848127d3c`.
All 43 installed Python files match that published artifact, and match the
prepared checkpoint after line-ending normalization (including ASTs).

ThermOpt requires `jaxprop >=0.6.0,<0.7.0`. The normal `.venv` now uses the
published package without `PYTHONPATH` substitution; its full suite passes
253 tests plus four existing expected failures. Its editable metadata and
dependency check also pass.

The regenerated lock changes only one package version: JaxProp 0.5.7 to 0.6.0.
Poetry also normalizes some requirement metadata; all other locked package
versions remain unchanged. The earlier description of a missing JaxProp lock
entry was incorrect: the incoming lock already contained 0.5.7.

A fresh isolated main-dependency installation succeeds, with all 38 production
dependency versions matching the lock. Pytest and its helpers are installed at
their locked versions. The final full run passes **262 tests**, with **four
unchanged expected failures**, zero failures and zero ordinary skips. Nine
new deterministic tests cover the benchmark comparison protocol without
wall-clock assertions. CI's Poetry 2.2.1 accepts and installs this lock; the lock
was generated with Poetry 2.5.1 in format 2.1. Existing project-metadata
deprecation warnings are retained. Hosted CI and other Python/platform profiles
have not been claimed as executed.

The finished [benchmark entry point](run.py) passes its five-mode `phase`
profile with 12 synchronized samples per backend and 100x100 working-fluid
tables. It verifies numerical agreement before recording timings and checks
that source, protocol and lock identities stay unchanged during the run.
The [portable starting checkpoint](checkpoints/pypi_060_phase.json) retains raw
samples, effective configurations, environment identities and agreement metrics.
This single starting cohort is not a between-phase performance acceptance claim.
The original prepared-wheel timings remain separate historical evidence.

Next is **Phase 2**, removing legacy constructor/YAML compatibility and migrating
maintained configurations to the new vocabulary. Later phases cover states,
components, lifecycle, outputs, examples and CI. See the
[implementation checklist](../z_JAXPROP_060_IMPLEMENTATION_PLAN.md) and
[benchmark instructions](README.md). No cycle JAX port is part of the plan.

[Full adoption manifest](RELEASE_060.json). Detailed logs are retained in ignored
`benchmark/results/pypi_060/`; no unrelated existing changes were committed.
