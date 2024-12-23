import argparse
import multiprocessing
import os
import sys
import time
from importlib.resources import files

from phonebook.src import qc, tasks
from phonebook.src.console import console

PACKAGE_DIR = files("phonebook")
TREE = PACKAGE_DIR / "resources/o1_cholera.annotated.pb"
REFERENCE = PACKAGE_DIR / "resources/reference.fasta"


def main(sysargs=None):
    if sysargs is None:
        sysargs = sys.argv[1:]

    parser = argparse.ArgumentParser(
        prog="phonebook",
        description="""\t\t     _               _           _   
\t\t ___| |_ ___ ___ ___| |_ ___ ___| |_ 
\t\t| . |   | . |   | -_| . | . | . | '_|
\t\t|  _|_|_|___|_|_|___|___|___|___|_,_|
\t\t|_|   

        Rapid classification of O1 Vibrio cholerae lineages.""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    threads = multiprocessing.cpu_count()
    cwd = os.getcwd()

    parser.add_argument(
        "query", nargs="*", help="Query fasta files of sequences to classify"
    )
    parser.add_argument(
        "--usher-tree",
        default=TREE,
        help="UShER Mutation Annotated Tree protobuf file to use instead of default tree",
    )
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
        "-t",
        "--threads",
        type=int,
        default=threads,
        help=f"Number of threads to use when possible. Default: all available cores, "
        f"{threads} detected on this machine",
    )

    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)

    args = parser.parse_args(sysargs)

    console.rule("[bold] Checking input and setting up analysis")

    # Checking input files
    query_file = qc.check_query_file(args.query)
    protobuf_file = qc.check_tree(args.usher_tree)

    outdir = qc.setup_outdir(args.outdir, cwd)
    outfile = outdir / args.outfile
    if outfile.exists():
        console.log(
            f"Warning: Specified output file {outfile} already exists and will be "
            f"overwritten."
        )
    tempdir = qc.setup_tempdir(args.tempdir, args.outdir, args.no_temp)

    if args.threads > threads:
        console.log(
            f"Error: Cannot request more cores than your computer has [{threads}]."
        )
        sys.exit(-99)
    threads = args.threads

    console.rule("[bold] Classifying input")
    with console.status("Processing...", spinner="bouncingBall"):
        console.log("Aligning sequences to reference")
        aln = tasks.align_sequences(query_file, REFERENCE, tempdir, threads)
        time.sleep(5)

        console.log("Converting alignment to VCF")
        vcf = tasks.convert_to_vcf(aln, REFERENCE, tempdir)
        time.sleep(5)

        console.log("Placing sequences into global phylogeny")
        results = tasks.classify_usher(vcf, protobuf_file, tempdir, threads)
        time.sleep(5)

        console.log("Parsing Usher results")
        tasks.usher_parsing(results, outfile)
        time.sleep(5)

    console.rule("[bold] Complete!")
    console.print(f"Lineage assignments are available in {outfile}.")


if __name__ == "__main__":
    main()
