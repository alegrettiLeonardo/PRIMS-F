"""Select an explicit, physically homogeneous OSC envelope study (stdlib only)."""
from __future__ import annotations

import json
from pathlib import Path

from compare_validation import REQUIRED_MODEL_KEYS, configuration_differences, source_signature
from validation_evidence import metadata

DEFAULT_PATTERNS = ("final_D_OSC_C3_BAROTROPIC_MATCHED_LAGGED_ROE_*",
                    "refined_D_OSC_C3_BAROTROPIC_MATCHED_LAGGED_ROE_*")


def select_runs(results_root, patterns=None, run_dirs=None):
    """Explicit selectors replace historical defaults; never silently add them."""
    if patterns and run_dirs:
        raise ValueError("Use either --pattern or --run-dir, not both")
    if run_dirs:
        paths = {Path(p).resolve() for p in run_dirs}
    else:
        paths = {p.resolve() for pattern in (patterns or DEFAULT_PATTERNS)
                 for p in Path(results_root).glob(pattern) if p.is_dir()}
    for path in sorted(paths):
        for filename in ("monitor.csv", "run_metadata.txt", "manifest.json"):
            if not (path / filename).is_file():
                raise ValueError(f"Selected run lacks {filename}: {path}")
    if not paths:
        raise ValueError("No runs matched the requested selectors")
    return sorted(paths)


def cohort_metadata(paths):
    """Prevent a combined refinement fit across physical models/source versions."""
    records = []
    for path in paths:
        pars = metadata(path / "run_metadata.txt")
        manifest = json.loads((path / "manifest.json").read_text())
        record = {"run": path.name, "case": pars.get("regime"),
                  "parameters": pars, "manifest": manifest}
        missing = [key for key in REQUIRED_MODEL_KEYS if key not in pars]
        if missing:
            raise ValueError(f"Missing model metadata for {path.name}: {', '.join(missing)}")
        if record["case"] != "D_OSC":
            raise ValueError(f"Envelope protocol requires D_OSC: {path.name}")
        if source_signature(record) is None:
            raise ValueError(f"Source hashes unavailable for {path.name}")
        if records:
            differences = configuration_differences(records[0], record)
            if source_signature(records[0]) != source_signature(record):
                differences.append("source_sha256")
            if differences:
                raise ValueError(f"Analyze separate cohorts; {path.name} differs in: {', '.join(differences)}")
        records.append(record)
    return records
