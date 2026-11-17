# Upstream follow-up: normalize NumPy static table options

Status on 17 November 2026: reproduced with the official JaxProp 0.6.0 wheel;
prepared locally, not submitted as an external issue. ThermOpt already normalizes
configuration scalars once at setup, so no repeated property-call checks are needed.

The following direct JaxProp construction creates a table but emits Equinox's
`UserWarning`: "A JAX array is being set as static!" The recorded origin is
`jaxprop/factory.py:253`. The observed issue is the warning/static metadata;
this check does not establish incorrect thermodynamic results.

```python
from pathlib import Path
import jax
import jaxprop as jp
import numpy as np

jax.config.update("jax_enable_x64", True)
fluid = jp.Fluid(
    "CO2",
    method="bicubic",
    table_options=dict(
        h_min=np.float64(110e3), h_max=1.3e6,
        p_min=8e6, p_max=39e6,
        N_h=np.int64(16), N_p=16, N_p_sat=4, coarse_step=1,
        table_dir=str(Path("benchmark/results/static_option_tables").resolve()),
    ),
)
```

Suggested fix in JaxProp: normalize appropriate static options to Python scalars
at the constructor boundary, preserving dynamic numerical leaves and native
input validation. Add warning-regression and equivalent-model cache-identity
tests for NumPy/Python options. Assess cold and warm complete-state consumers;
do not add metadata normalization to `get_state`.

Raw local evidence: `benchmark/results/phase8/static_metadata_warning.json`.
This issue is separate from the already tested but unpublished metadata fix
`52c98b8`; neither a local commit nor a prepared wheel establishes PyPI publication.
