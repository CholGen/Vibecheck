import subprocess
import sys
from pathlib import Path
from typing import Union

from vibecheck.src.console import console


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
