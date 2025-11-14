"""
Helper modules extracted from the legacy bin/feature_main.py script.

Exposes the CLI parsing and pipeline entry point so callers can
import and reuse them programmatically.
"""

from .cli import parse_args  # noqa: F401
from .pipeline import run_pipeline  # noqa: F401
from .transcript_input import filter_low_impact_transcripts  # noqa: F401

