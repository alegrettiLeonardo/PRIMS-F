"""Explicit overrides let a repeat campaign use a fresh output tree."""
from pathlib import Path
import os
ROOT=Path(os.environ.get('PRIMSF_BOUNDARY_ROOT','/home/leonardoab/Documentos/PRIMS-F-boundary-audit-20261119')).resolve()
OUT=Path(os.environ.get('PRIMSF_BOUNDARY_OUT','/home/leonardoab/Documentos/PRIMS-F-boundary-results-20261119')).resolve()
