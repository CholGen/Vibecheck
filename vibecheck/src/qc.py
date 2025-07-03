import sys
import tempfile
from pathlib import Path
from typing import Tuple, Union

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


def check_query_file(
    file_paths: list[str],
) -> Tuple[Union[Path, Tuple[Path, Path]], bool]:
    """Validates a list of sequence file paths and determines if they represent FASTA or FASTQ files.

    Parameters
    ----------
    file_paths : List[str]
        List of file paths (relative or absolute) to sequence files

    Returns
    -------
    Tuple[bool, Union[Path, Tuple[Path, Path]]]
        First element is either:
        - A single Path object for FASTA file
        - A tuple of 2 Path objects for FASTQ files
        Second element is True if the usher pipeline should be used to classify sequences in a fasta file, or false if
        the freyja pipeline should be used to classify reads in fastq files.

    Raises
    ------
    SystemExit
        If validation fails due to:
        - Empty input list
        - Mix of FASTA and FASTQ files
        - Multiple FASTA files
        - More than 2 FASTQ files
        - Files don't exist
    """
    if not file_paths:
        console.log("Error: No query file supplied. Please supply one.")
        sys.exit(-3)

    # Expand and validate all paths
    expanded_paths = set()
    for file_path in file_paths:
        path = Path(file_path).expanduser().resolve()
        if not path.exists():
            console.log(f"Error: query file does not exist: {path}")
            sys.exit(-2)
        expanded_paths.add(path)

    # Define valid extensions
    fasta_extensions = {".fasta", ".fa", ".fasta.gz", ".fa.gz"}
    fastq_extensions = {".fastq", ".fq", ".fastq.gz", ".fq.gz"}

    # Check file extensions
    fasta_files = {
        path
        for path in expanded_paths
        if any(str(path).lower().endswith(ext) for ext in fasta_extensions)
    }
    fastq_files = {
        path
        for path in expanded_paths
        if any(str(path).lower().endswith(ext) for ext in fastq_extensions)
    }

    # Validate file combinations
    if fasta_files and fastq_files:
        console.log("Error: Cannot mix FASTA and FASTQ files for query")
        sys.exit(-4)

    if len(fasta_files) > 1:
        console.log(
            f"Error: More than 1 fasta file provided: {fasta_files}. Please concatenate fasta files or supply only one."
        )
        sys.exit(-5)

    if len(fastq_files) > 2:
        console.log(f"Error: More than 2 fastq files provided: {fastq_files}")
        sys.exit(-6)

    if len(expanded_paths) != len(fasta_files) + len(fastq_files):
        invalid_files = expanded_paths - (fasta_files | fastq_files)
        console.log(
            f"Error: File extensions not recognized: {', '.join(str(p) for p in invalid_files)}. Only fasta ({fasta_extensions}) and fastq ({fastq_extensions}) files are recognized."
        )
        sys.exit(-7)

    # Return tuple of (is_fasta, file_paths)
    if fasta_files:
        use_file = next(iter(fasta_files))
        console.print(f"Using query fasta file: {str(use_file)}")
        return use_file, True  # Return single Path for FASTA
    else:
        use_file = tuple(sorted(fastq_files))
        if len(use_file) == 1:
            use_file = (*use_file, None)
        console.print(f"Using query fastq files: {list(map(str, use_file))}")
        return use_file, False  # Return tuple of Paths for FASTQ


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
