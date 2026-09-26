# Clean-clone verification

Use Python 3.12 (the CI version). From a terminal:

```sh
git clone https://github.com/ShaneThomas06/single-cell-ibd-immune-atlas.git
cd single-cell-ibd-immune-atlas
python -m venv .venv
```

Activate with `source .venv/bin/activate` on Linux/macOS or
`.\.venv\Scripts\Activate.ps1` in Windows PowerShell, then run:

```sh
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
python -m pip check
python -m ruff check src tests app
python -m pytest -q
immune-atlas-pseudobulk --help
```

Expected: lint and tests pass, and the final command prints usage and exits
successfully. These checks require no study downloads, GPU, or model weights;
installation requires internet access. Tests use small synthetic fixtures and
temporary outputs. They do not regenerate or overwrite tracked scientific results.
CI runs the same checks on every push and pull request.

This is a software smoke check, not a reproduction of the complete study.
For an exact environment record, save `python -m pip freeze > environment.txt`
alongside your own run outputs, and record the repository commit and input checksums.
Version bounds describe supported dependency families, not an exact lockfile.

For study data preparation and donor-aware analysis, follow the README.
The synthetic unit tests do not require downloaded patient data.
