import sys
import tempfile
from pathlib import Path

from vibecheck.src.console import console


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
        if files[0].startswith("~"):
            file_path = file_path.expanduser()
        if not file_path.exists():
            console.log(f"Error: No such file: {file_path}")
            sys.exit(-2)
        console.print(f"Using query file: {file_path}")

        return file_path

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
