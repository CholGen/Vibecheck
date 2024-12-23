import sys
import tempfile
from pathlib import Path

from phonebook.src.console import console


def check_query_file(files: list[str]) -> Path:
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
    usher_tree = Path(usher_tree)

    if not usher_tree.exists():
        console.log(
            f"Error: {usher_tree} was not found. Please check path supplied to `--usher-tree`."
        )
        sys.exit(-2)

    console.print(f"Usher tree found: {usher_tree}")
    return usher_tree


def setup_outdir(outdir: str, cwd: str) -> Path:
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


def setup_tempdir(tempdir: str, outdir: str, no_temp: bool) -> Path:
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
            with open(tempdir / "test.txt", "w") as f:
                f.write("Test")
        except OSError:
            console.log(f"Error: cannot write to temp directory: {tempdir}.")
            sys.exit(-4)
    console.print(f"Intermediate files will be written to: {tempdir}")
    (tempdir / "logs").mkdir()
    return tempdir
