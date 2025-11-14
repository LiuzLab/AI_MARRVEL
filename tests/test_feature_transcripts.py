import pandas as pd
import sys
from pathlib import Path


def _ensure_bin_on_path():
    repo_root = Path(__file__).resolve().parents[1]
    bin_path = repo_root / "bin"
    if str(bin_path) not in sys.path:
        sys.path.insert(0, str(bin_path))


_ensure_bin_on_path()

from feature import filter_low_impact_transcripts  # noqa: E402


def test_filter_low_impact_transcripts_prefers_high_or_moderate():
    transcript_df = pd.DataFrame(
        {
            "#Uploaded_variation": ["var1", "var1", "var2"],
            "IMPACT": ["LOW", "HIGH", "LOW"],
        }
    )

    filtered = filter_low_impact_transcripts(transcript_df)

    assert list(filtered["#Uploaded_variation"]) == ["var1", "var2"]
    assert list(filtered["IMPACT"]) == ["HIGH", "LOW"]


def test_filter_low_impact_transcripts_keeps_low_when_only_option():
    transcript_df = pd.DataFrame(
        {
            "#Uploaded_variation": ["var1", "var2"],
            "IMPACT": ["LOW", "LOW"],
        }
    )

    filtered = filter_low_impact_transcripts(transcript_df)

    assert filtered.equals(transcript_df)
