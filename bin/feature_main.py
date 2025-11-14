#!/usr/bin/env python3.8
"""CLI entry point for the feature pipeline."""

from feature import parse_args, run_pipeline


def main() -> None:
    args = parse_args()
    run_pipeline(args)


if __name__ == "__main__":
    main()

