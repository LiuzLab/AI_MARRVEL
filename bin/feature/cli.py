import argparse
from typing import Optional, Sequence


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("-varFile", "--varFile", help="Provide input variant file")
    parser.add_argument(
        "-inFileType",
        "--inFileType",
        help="Provide type of input file: vcf, vepAnnotTab",
    )
    parser.add_argument(
        "-patientFile", "--patientFile", help="Provide HPO IDs"
    )  # currently not used
    parser.add_argument(
        "-patientFileType",
        "--patientFileType",
        help="Provide type of file: one, two",
    )  # currently not used
    parser.add_argument(
        "-patientHPOsimiOMIM",
        "--patientHPOsimiOMIM",
        help="Provide patient HPO similarity file-OMIM",
    )
    parser.add_argument(
        "-patientHPOsimiHGMD",
        "--patientHPOsimiHGMD",
        help="Provide patient HPO similarity file-HGMD",
    )
    parser.add_argument(
        "-diseaseInh", "--diseaseInh", help="Provide disease Inheritance: AD, AR, XD, XR"
    )
    parser.add_argument(
        "-genomeRef", "--genomeRef", help="Provide genome ref: hg19, hg38"
    )
    parser.add_argument(
        "-enableLIT",
        "--enableLIT",
        action="store_true",
        default=False,
        help="Enable low-impact transcripts",
    )
    return parser


def check_user_args(args: argparse.Namespace) -> None:
    # Placeholder for future validation mirroring the legacy stub.
    return None


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = build_parser()
    args = parser.parse_args(argv)
    check_user_args(args)
    return args
