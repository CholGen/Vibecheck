import argparse
import subprocess
import sys
from pathlib import Path
from typing import Union

from vibecheck.src.console import console


def add_cli_arguments(
    parser: argparse.ArgumentParser,
    default_tree: Path,
    default_barcodes: Path,
    default_aliases: Path,
    threads: int,
) -> None:

    parser.add_argument("query", nargs="*", help="Query sequences to classify")
    parser.add_argument(
        "-o",
        "--outdir",
        action="store",
        default=".",
        help="Output directory. Default: current working directory",
    )
    parser.add_argument(
        "--outfile",
        action="store",
        default="lineage_report.csv",
        help="Optional output file name. Default: lineage_report.csv",
    )
    parser.add_argument(
        "--tempdir",
        action="store",
        help="Specify where you want the temp stuff to go. Default: $TMPDIR",
    )
    parser.add_argument(
        "--no-temp",
        action="store_true",
        help="Output all intermediate files, for dev purposes.",
    )
    parser.add_argument(
        "--lineage-aliases",
        default=str(default_aliases),
        help="Path to CSV file containing lineage aliases. File must have two columns: `alias` and `lineage`.",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=threads,
        help=f"Number of threads to use when possible. Default: all available cores, "
        f"{threads} detected on this machine",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="store_true",
        help="Prints the version of Vibecheck and exits.",
    )

    sequence_args = parser.add_argument_group("Sequence-based classification")
    read_args = parser.add_argument_group("Read-based classification")

    sequence_args.add_argument(
        "-u",
        "--usher-tree",
        default=str(default_tree),
        help="UShER Mutation Annotated Tree protobuf file to use instead of default tree",
    )
    sequence_args.add_argument(
        "-m",
        "--max-ambiguity",
        type=float,
        default=0.3,
        help="Maximum number of ambiguous bases a sequence can have before its filtered from the analysis. Default: 0.3",
    )

    read_args.add_argument(
        "-b",
        "--barcodes",
        type=str,
        default=str(default_barcodes),
        help="Feather formatted lineage barcodes to use instead of default O1 barcodes",
    )
    read_args.add_argument(
        "-s",
        "--subsample",
        type=float,
        default=0.2,
        help="Fraction of reads to use in classification. Default: 0.2",
    )
    read_args.add_argument(
        "--no-subsample",
        action="store_true",
        help="Do not subsample reads. Default: False",
    )


def run_command(
    command: str,
    log_path: Union[Path, None] = None,
    error_message: str = "Command failed",
) -> None:
    """Runs a shell command, logs output to a file (if provided), and exits on failure.

    Parameters
    ----------
    command: str
        Bash command to run.
    log_path: Union[Path, None], optional
        Path to the log file to write to. If not provided, then stdout and stderr is ignored.
    error_message: str, optional
        Message to display on failure. A generic message is used if not provided.
    """
    result = subprocess.run(
        command, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )

    # Handle errors
    if result.returncode != 0:
        console.print(result.stdout)
        console.log(f"Error: {error_message}. See above for details.")
        sys.exit(result.returncode)

    # Write command output to log file if a log_path is provided
    if log_path:
        with open(log_path, "w") as logfile:
            logfile.write(result.stdout)
