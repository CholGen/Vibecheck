import argparse
import multiprocessing
import os
import sys
from importlib.resources import files

from vibecheck.src import qc, usher_tasks, utilities, freyja_tasks
from vibecheck.src.console import console
from vibecheck import __version__

PACKAGE_DIR = files("vibecheck")
TREE = PACKAGE_DIR / "resources/o1_cholera.no_missing.pb"
REFERENCE = PACKAGE_DIR / "resources/reference.fasta"
BARCODES = PACKAGE_DIR / "resources/o1_barcodes.feather"
ALIASES = PACKAGE_DIR / "resources/o1_aliases.csv"

# Alternative names: Cholin, choline, Colin, vibecheck, vibcheck, vibration, vibrator,
# marlin, vibrato, chool-aid


def main(sysargs=None):
    if sysargs is None:
        sysargs = sys.argv[1:]

    parser = argparse.ArgumentParser(
        prog="vibecheck",
        description=r"""██╗   ██╗██╗██████╗ ███████╗ ██████╗██╗  ██╗███████╗ ██████╗██╗  ██╗
██║   ██║██║██╔══██╗██╔════╝██╔════╝██║  ██║██╔════╝██╔════╝██║ ██╔╝
██║   ██║██║██████╔╝█████╗  ██║     ███████║█████╗  ██║     █████╔╝ 
╚██╗ ██╔╝██║██╔══██╗██╔══╝  ██║     ██╔══██║██╔══╝  ██║     ██╔═██╗ 
 ╚████╔╝ ██║██████╔╝███████╗╚██████╗██║  ██║███████╗╚██████╗██║  ██╗
  ╚═══╝  ╚═╝╚═════╝ ╚══════╝ ╚═════╝╚═╝  ╚═╝╚══════╝ ╚═════╝╚═╝  ╚═╝

        Rapid classification of O1 Vibrio cholerae lineages.""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    threads = multiprocessing.cpu_count()
    cwd = os.getcwd()

    utilities.add_cli_arguments( parser=parser, default_tree=TREE, default_barcodes=BARCODES, default_aliases=ALIASES, threads=threads)

    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)

    args = parser.parse_args(sysargs)

    if args.version:
        print(__version__)
        sys.exit(0)

    console.rule("[bold] Checking input and setting up analysis")

    # Checking requested threads.
    threads = qc.check_threads(args.threads, threads)

    # Check input files
    query_file, use_usher = qc.check_query_file(args.query)

    aliases = {}
    if args.lineage_aliases:
        aliases = qc.check_lineage_aliases(args.lineage_aliases)

    # Check directories
    outdir = qc.setup_outdir(args.outdir, cwd)
    outfile = outdir / args.outfile
    if outfile.exists():
        console.log(
            f"Warning: Specified output file {outfile} already exists and will be "
            f"overwritten."
        )
    tempdir = qc.setup_tempdir(args.tempdir, outdir, args.no_temp)

    if use_usher:
        # Checking UShER pipeline.
        protobuf_file = qc.check_tree(args.usher_tree)
        max_ambiguity = qc.check_max_ambiguity(args.max_ambiguity)
        console.rule("[bold] Classifying input sequences")
        with console.status("Processing...", spinner="bouncingBall"):
            usher_tasks.run_pipeline(
                query_file=query_file,
                protobuf_tree=protobuf_file,
                reference=REFERENCE,
                max_ambiguity=max_ambiguity,
                lineage_aliases=aliases,
                tempdir=tempdir,
                outfile=outfile,
                threads=threads,
            )
    else:
        # Checking Freyja pipeline
        barcodes = qc.check_barcodes(args.barcodes)
        subsample_frac = qc.check_subsampling_frac(args.subsample)
        console.rule("[bold] Classifying input reads")
        with console.status("Processing...", spinner="bouncingBall"):
            freyja_tasks.run_pipeline(
                reads=query_file,
                barcodes=barcodes,
                reference=REFERENCE,
                subsample_fraction=subsample_frac,
                no_subsample=args.no_subsample,
                lineage_aliases=aliases,
                tempdir=tempdir,
                outfile=outfile,
                threads=threads,
            )

    console.rule("[bold] Complete!")
    console.print(f"Lineage assignments are available in {outfile}.")


if __name__ == "__main__":
    main()
