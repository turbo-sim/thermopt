# Retired bicubic transition experiment

These three original scripts are preserved byte for byte with their SHA256
identities in `sources.json`. They document the earlier adapter/monkey-patching
investigation and are not maintained launchers or supported interfaces.
The local `.gitattributes` disables text conversion for these original scripts,
so Windows `core.autocrlf` also preserves the recorded byte identities.

Current ThermOpt constructs native JaxProp evaluators directly. Use
`benchmark/native_backends.py`, `benchmark/run.py` and the ordinary native YAML
variants under `case_studies/`. Dated development results and their original
source hashes remain unchanged. Do not use this archive to generate new
acceptance measurements; its original relative paths and old API assumptions
are intentionally not rewritten.
