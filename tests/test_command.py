import shutil
from pathlib import Path
from unittest.mock import patch

import pytest

from phonebook.command import main

expected_resultA = "hash,lineage,conflict,usher_note\nAfrica|KEN|ERR037738|T10|2010-01-01,T10,0.0,Usher placements: T10(1/1)\n"
expected_resultB = "hash,lineage,conflict,usher_note\nAfrica|ZAF|ERS14903183|T15|2023-01-01,T15,0.0,Usher placements: T15(1/1)\n"
expected_resultC = "hash,lineage,conflict,usher_note\nAfrica|TZA|SAMN19110428|T13|2017-01-01,T13,0.0,Usher placements: T13(1/1)\n"


@patch("argparse.ArgumentParser.print_help")
def test_empty_args_exit(mock_print_help):
    sys_argv = []
    with pytest.raises(SystemExit) as e:
        main(sys_argv)
    assert e.value.code == -1
    mock_print_help.assert_called_once()


def test_too_many_threads():
    sys_argv = ["--threads", "99", "tests/example_fasta/ERR037738.fasta"]
    with pytest.raises(SystemExit) as e:
        main(sys_argv)
    assert e.value.code == -99


def test_zero_threads():
    sys_argv = ["--threads", "0", "tests/example_fasta/ERR037738.fasta"]
    with pytest.raises(SystemExit) as e:
        main(sys_argv)
    assert e.value.code == -98


def test_odd_threads():
    sys_argv = ["--threads", "foo", "tests/example_fasta/ERR037738.fasta"]
    with pytest.raises(SystemExit):
        main(sys_argv)


@pytest.fixture
def inputA():
    outdir = Path("tests/inputA")
    sys_argv = [
        "--outdir",
        str(outdir),
        "--threads",
        "4",
        "tests/example_fasta/ERR037738.fasta",
    ]
    main(sys_argv)
    yield outdir

    if outdir.exists():
        shutil.rmtree(outdir)


@pytest.fixture
def inputB():
    outdir = Path("tests/inputB")
    sys_argv = [
        "--outdir",
        str(outdir),
        "--threads",
        "4",
        "tests/example_fasta/ERS14903183.fasta",
    ]
    main(sys_argv)
    yield outdir

    if outdir.exists():
        shutil.rmtree(outdir)


@pytest.fixture
def inputC():
    outdir = Path("tests/inputC")
    sys_argv = [
        "--outdir",
        str(outdir),
        "--threads",
        "4",
        "tests/example_fasta/SAMN19110428.fasta",
    ]
    main(sys_argv)
    yield outdir

    if outdir.exists():
        shutil.rmtree(outdir)


@pytest.fixture
def inputA_notemp():
    outdir = Path("tests/inputA_notemp")
    sys_argv = [
        "--outdir",
        str(outdir),
        "--no-temp",
        "tests/example_fasta/ERR037738.fasta",
        "--threads",
        "4",
    ]
    main(sys_argv)
    yield outdir

    if outdir.exists():
        shutil.rmtree(outdir)


# Try to refactor this in the future.
# @pytest.mark.parametrize("test_input,expected", [("3+5", 8), ("2+4", 6), ("6*9", 42)])
def test_input(inputA):
    results = inputA / "lineage_report.csv"
    assert results.exists()
    assert results.read_text() == expected_resultA


def test_inputB(inputB):
    results = inputB / "lineage_report.csv"
    assert results.exists()
    assert results.read_text() == expected_resultB


def test_inputC(inputC):
    results = inputC / "lineage_report.csv"
    assert results.exists()
    assert results.read_text() == expected_resultC


def test_inputA_notemp_exists(inputA_notemp):
    expected_files = [
        "lineage_report.csv",
        "logs/usher.txt",
        "logs/minimap2.txt",
        "mapped.sam",
        "alignment.fasta",
        "sequences.aln.vcf",
        "sequences.withref.fa",
        "clades.txt",
    ]
    for file in expected_files:
        file_path = inputA_notemp / file
        assert file_path.exists()
