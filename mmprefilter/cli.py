"""Command line. Analysis choices live in the config file."""

import argparse
import sys

from . import __version__
from .config import STAGES, ConfigError, load_config
from .pipeline import run_config


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="mmprefilter",
        description=(
            "Score MarkerMiner alignments and write the subset to open in Geneious. "
            "Paths, thresholds, and header fixes are read from the config file."
        ),
    )
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command")
    run = sub.add_parser("run", help="Run the stages enabled in a config file.")
    run.add_argument("--config", required=True, help="Path to the JSON config.")
    run.add_argument(
        "--from",
        dest="start",
        default=None,
        choices=STAGES,
        help="First stage to run. Default: the first stage.",
    )
    run.add_argument(
        "--to",
        dest="end",
        default=None,
        choices=STAGES,
        help="Last stage to run. Default: report.",
    )
    args = parser.parse_args(argv)
    if args.command != "run":
        parser.print_help()
        return 2
    try:
        cfg = load_config(args.config)
        return run_config(cfg, start=args.start, end=args.end)
    except ConfigError as exc:
        sys.stderr.write("Error: {0}\n".format(exc))
        return 1


if __name__ == "__main__":
    sys.exit(main())
