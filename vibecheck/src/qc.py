import sys
import tempfile
from pathlib import Path
from typing import Tuple

from vibecheck.src.console import console


def check_threads(requested: int, threads: int) -> int:
    """Checks whether the requested number of threads is possible.
    Exits if not possible.

    Parameters
    ----------
    requested: int
    threads: int
        Number of cores computer has access to.

    Returns
    -------
    int
        Number of threads to use.
    """
    if requested > threads:
        console.log(
            f"Error: Cannot request more cores than your computer has [{threads}]."
        )
        sys.exit(-99)
    elif requested < 1:
        console.log("Error: Must request at least one core using the --threads option.")
        sys.exit(-98)
    return requested


def check_parse_float_fraction(value: float, name: str) -> float:
    """Conforms fraction value is between 0 and 1. Exits if parsing is not possible.

    Parameters
    ----------
    value: float
        Putative value to check.
    name: str
        Name of the value, to be used in error message.

    Returns
    -------
    float
        Value between 0 and 1.
    """
    if value < 0:
        console.log("Error: {name} cannot be negative.")
        sys.exit(-97)
    elif value > 100:
        console.log("Error: {name} must represent a fraction")
        sys.exit(-97)
    elif value > 1:
        value /= 100
    return value


def check_max_ambiguity(max_ambiguity: float) -> float:
    """Checks maximum abiguity value can be parsed to a fraction.

    Parameters
    ----------
    max_ambiguity: float

    Returns
    -------
    float
    """
    return_value = check_parse_float_fraction(max_ambiguity, "Maximum ambiguity")
    return return_value


def check_subsampling_frac(subsampling_fraction: float) -> float:
    """Checks subsampling value can be parsed to a fraction.

    Parameters
    ----------
    subsampling_fraction: float

    Returns
    -------
    float
    """
    return_value = check_parse_float_fraction(
        subsampling_fraction, "Subsampling fraction"
    )
    return return_value


def check_query_file(files: list[str]) -> Tuple[Path, bool]:
    """Confirms the existance of a single query file.
    Parameters
    ----------
    files: list[str]
        list of input files to check.

    Returns
    -------
    Path
        Location of input fasta file.
    bool
        Indicates whether the input represents a fasta file of consensus genomes or paired-end fastq files.
    """
    if len(files) > 1:
        console.log(
            f"[bold red]Error: Too many query (input) fasta files supplied: "
            f"{files}\nPlease supply one only."
        )
        sys.exit(-1)
    try:
        file_path = Path(files[0])
        if files[0].startswith("~"):
            file_path = file_path.expanduser()
        if not file_path.exists():
            console.log(f"Error: No such file: {file_path}")
            sys.exit(-2)
        console.print(f"Using query file: {file_path}")

        return file_path, True

    except IndexError:
        console.log(f"Error: No query file supplied: {files}\nPlease supply one.")
        sys.exit(-3)


def check_tree(usher_tree: str) -> Path:
    """Confirms the existance of protobuf tree

    Parameters
    ----------
    usher_tree: str
        Puntative location of protobuf tree.

    Returns
    -------
    Path
        Location of confirmed protobuf tree.
    """
    usher_path = Path(usher_tree)
    if usher_tree.startswith("~"):
        usher_path = usher_path.expanduser()

    if not usher_path.exists():
        console.log(
            f"Error: {usher_tree} was not found. Please check path supplied to `--usher-tree`."
        )
        sys.exit(-4)

    console.print(f"Using Usher tree: {usher_path}")
    return usher_path


def check_barcodes(barcodes: str) -> Path:
    """Confirms the existance of barcodes.

    Parameters
    ----------
    barcodes: str
        Puntative location of barcodes.

    Returns
    -------
    Path
        Location of confirmed barcodes.
    """
    barcodes_path = Path(barcodes)
    if barcodes.startswith("~"):
        barcodes_path = barcodes_path.expanduser()
    if not barcodes_path.exists():
        console.log(
            f"Error: {barcodes} was not found. Please check path supplied to `--barcodes`."
        )
        sys.exit(-9)
    console.print(f"Using lineage barcodes: {barcodes_path}")
    return barcodes_path


def setup_outdir(outdir: str, cwd: str) -> Path:
    """Identifies and/or creates output directory for analysis. If supplied output
    directory is a relative path, append it to the current working directory to create
    an absolute path. Generate working directory if necessary and possible. Exits if
    directory cannot be created.

    Parameters
    ----------
    outdir: str
        Putative output directory.
    cwd: str
        Current working directory.

    Returns
    -------
    Path
        Location of output directory.
    """
    if outdir:
        outdir_path = Path(outdir)
        if outdir.startswith("~"):
            outdir_path = outdir_path.expanduser()

        if not outdir_path.is_absolute():
            outdir_path = Path(cwd, outdir)
        if not outdir_path.exists():
            try:
                outdir_path.mkdir()
            except FileNotFoundError:
                console.log(
                    f"Error: Unable to locate or create {outdir_path} directory."
                )
                sys.exit(-5)
            except OSError:
                console.log(
                    f"Error: Unable to locate or create {outdir_path} directory."
                )
                sys.exit(-5)
    else:
        outdir_path = Path(cwd)
    console.print(f"Using output directory: {outdir_path}")
    return outdir_path


def setup_tempdir(tempdir: str, outdir: Path, no_temp: bool) -> Path:
    """Identifies and/or creates temporary directory for analysis. If no temporary
    files are to be created, use the output directory as the temporary directory. Else,
    create output directory if necessary. Exits if directory cannot be created.

    Parameters
    ----------
    tempdir: str
        Putative temporary directory.
    outdir: Path
        Established output directory.
    no_temp: bool
        Whether to create temporary directory. If false, use the output directory.

    Returns
    -------
    Path
        Location of temporary directory.
    """
    if no_temp:
        tempdir_path = Path(outdir)
    elif tempdir:
        tempdir_path = Path(tempdir)
        if tempdir.startswith("~"):
            tempdir_path = tempdir_path.expanduser()
        try:
            if not tempdir_path.exists():
                tempdir_path.mkdir()
        except FileNotFoundError:
            console.log(f"Error: cannot create temp directory: {tempdir}.")
            sys.exit(-6)
        except OSError:
            console.log(f"Error: cannot create temp directory: {tempdir}.")
            sys.exit(-6)
    else:
        tempdir_path = Path(tempfile.mkdtemp())
        try:
            if not tempdir_path.exists():
                tempdir_path.mkdir()
        except FileNotFoundError:
            console.log(f"Error: cannot create temp directory: {tempdir}.")
            sys.exit(-7)
        except OSError:
            console.log(f"Error: cannot create temp directory: {tempdir}.")
            sys.exit(-7)

        try:
            with open(tempdir_path / "tests.txt", "w") as f:
                f.write("Test")
        except OSError:
            console.log(f"Error: cannot write to temp directory: {tempdir}.")
            sys.exit(-8)
    console.print(f"Intermediate files will be written to: {tempdir_path}")
    (tempdir_path / "logs").mkdir(exist_ok=True)
    return tempdir_path
