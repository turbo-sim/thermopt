# JaxProp metadata fix and ThermOpt verification

The Phase 3 source checkpoint is `fb6745a`. JaxProp's corresponding fix is
`52c98b8`: aligned immutable metadata comparisons avoid repeated dictionary
construction during JIT lookup. Mapping order/hash/pickle behavior, real setting
changes, numerical validation and all thermodynamic equations remain intact.
There is no global sharing of mutable evaluators or application-side adapter.

The fix is tested as an **unpublished local wheel**, retaining version 0.6.0.
It is identified by SHA256
`bd9d3cdf7ee73ffb096a615f1ed34abd4ced141a18105f09e29c256101e40d71`.
The official PyPI wheel remains
`0a8dfb68ca4d21cf97b7483e1da776b946cec108a94d50ebb770bb9848127d3c`.
`configuration.py` is the only package source delta; all 43 installed Python files
match the local artifact. The isolated candidate environment matches the locked
numerical dependencies, and the official-wheel environment is preserved.

## Correctness

- 17 new JaxProp metadata tests pass, including collisions, reordered mappings
  and JIT reuse/retracing. The regression fails before the fix.
- JaxProp's full HEOS/REFPROP run passes 3,044 tests, with two optional C checks
  initially skipped. Both pass with the installed compiler: 3,046 distinct cases
  are verified, with no remaining unexecuted test from that run.
- ThermOpt passes 396 tests plus its four existing strict expected failures,
  including all standard optimizations and the native state contracts.
- All six complete-cycle comparison cohorts and the targeted follow-up preserve
  exactly identical physical snapshots. No scientific references are rewritten.

## Performance and remaining observations

All fourteen standalone full-state scalar/batch and compiled-consumer cases meet
their original budgets. The equivalent independent-model control changes from
an earlier 28.8% second-model penalty to 7.0%. All twelve paired application
comparisons with separate model instances meet 1.10 after the fix.

Three balanced fresh-process wheel pairs measure the same current ThermOpt code,
CPU affinity and physical inputs. Nine of ten aggregate workloads meet 1.10.
Nitrogen bicubic initially records 1.117. A separate three-pair repeat uses
36 samples and a fixed Python hash seed, giving 0.999 source and 1.004 bicubic.
The initial flag is retained; changing seed, sample count and time does not
establish a single cause for its shift. Construction/table generation remains
separate from warmed evaluations.

Two older absolute comparisons against Phase 2's stored medians remain above
1.10: all-stream recuperated bicubic (1.127) and nitrogen refrigeration bicubic
(1.220). These historical costs have not disappeared or been waived. The fresh
matched library comparison and its controls support the targeted fix; historical
absolute timing closure remains final verification work. Phase 4 proceeds under
the user's instruction with these measurements explicitly retained.

[METADATA_FIX.json](METADATA_FIX.json) retains all raw samples, original flags,
repeat controls, test hashes and artifact identities. The
[artifact-specific checkpoint](checkpoints/metadata_fix.json) retains the current
three-cohort medians and flags. It is not a clean all-green historical checkpoint.
Use the recorded fixed-wheel environment for matched Phase 4 measurements;
publishing a future JaxProp patch and updating the lock remains a separate step.

The next implementation work is Phase 4's independent component physics checks
and measured polytropic reporting improvements. Full-cycle JAX conversion remains
a later, separately planned sweep.
