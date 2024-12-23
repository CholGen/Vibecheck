import sys
import tempfile
from pathlib import Path

from phonebook.src.console import console


def check_query_file(files: list[str]) -> Path:
    """Confirms the existance of a single query file.
    Parameters
    ----------
    files: list[str]
        list of input files to check.

    Returns
    -------
    Path
        Location of input fasta file.
    """
    if len(files) > 1:
        console.log(
            f"[bold red]Error: Too many query (input) fasta files supplied: "
            f"{files}\nPlease supply one only."
        )
        sys.exit(-1)
    try:
        file_path = Path(files[0])
        if not file_path.exists():
            console.log(f"Error: No such file: {file_path}")
            sys.exit(-1)
        console.print(f"Query file found: {file_path}")

        return file_path

    except IndexError:
        console.log(f"Error: No query file supplied: {files}\nPlease supply one.")
        sys.exit(-1)


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
    usher_tree = Path(usher_tree)

    if not usher_tree.exists():
        console.log(
            f"Error: {usher_tree} was not found. Please check path supplied to `--usher-tree`."
        )
        sys.exit(-2)

    console.print(f"Usher tree found: {usher_tree}")
    return usher_tree


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
        outdir_path = Path(outdir).expanduser()
        if not outdir_path.is_absolute():
            outdir_path = Path(cwd, outdir)
        if not outdir_path.exists():
            try:
                outdir_path.mkdir()
            except FileNotFoundError:
                console.log(
                    f"Error: Unable to locate or create {outdir_path} directory."
                )
                sys.exit(-3)
    else:
        outdir_path = Path(cwd)
    console.print(f"Output directory created: {outdir_path}")
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
        tempdir = Path(outdir)
    elif tempdir:
        tempdir = Path(tempdir)
        try:
            if not tempdir.exists():
                tempdir.mkdir()
        except FileNotFoundError:
            console.log(f"Error: cannot create temp directory: {tempdir}.")
            sys.exit(-4)
    else:
        tempdir = Path(tempfile.mkdtemp())
        try:
            if not tempdir.exists():
                tempdir.mkdir()
        except FileNotFoundError:
            console.log(f"Error: cannot create temp directory: {tempdir}.")
            sys.exit(-4)
        try:
            with open(tempdir / "tests.txt", "w") as f:
                f.write("Test")
        except OSError:
            console.log(f"Error: cannot write to temp directory: {tempdir}.")
            sys.exit(-4)
    console.print(f"Intermediate files will be written to: {tempdir}")
    (tempdir / "logs").mkdir()
    return tempdir
