import os
import shutil
import subprocess
from pathlib import Path
from unittest.mock import patch

import pytest

from vibecheck.command import main

expected_resultAll = (
    "sequence_id,qc_status,qc_notes,lineage,confidence,classification_notes\n"
    "Africa|KEN|ERR037738|T10|2010-01-01,pass,Ambiguous_content:0.01%,T10,1.0,Usher placements: T10(1/1)\n"
    "Africa|TZA|SAMN19110428|T13|2017-01-01,pass,Ambiguous_content:0.01%,T13,1.0,Usher placements: T13(1/1)\n"
    "Africa|ZAF|ERS14903183|T15|2023-01-01,pass,Ambiguous_content:0.04%,T15,1.0,Usher placements: T15(1/1)\n"
)


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
def inputAll():
    print(os.getcwd())
    outdir = Path("tests/inputA/")
    outdir.mkdir(exist_ok=True)
    input_alignment = outdir / "all.fasta"
    subprocess.run(
        # C10.fasta is a multi-contig assembly fixture used by the --assemblies
        # tests, not a single-sequence sample, so it's excluded here.
        "cat tests/example_fasta/ERR037738.fasta "
        "tests/example_fasta/SAMN19110428.fasta "
        f"tests/example_fasta/ERS14903183.fasta > {input_alignment}",
        shell=True,
    )
    sys_argv = [
        "--outdir",
        str(outdir),
        "--threads",
        "4",
        str(input_alignment),
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


def test_input(inputAll):
    results = inputAll / "lineage_report.csv"
    assert results.exists()
    assert results.read_text() == expected_resultAll


def test_inputA_notemp_exists(inputA_notemp):
    expected_files = [
        "lineage_report.csv",
        "logs/usher.txt",
        "logs/minimap2.txt",
        "mapped.sam",
        "alignment.fasta",
        "sequences.aln.vcf",
        "sequences.withref.fa",
        "alignment.filtered.fasta",
        "seq_status.csv",
        "clades.txt",
    ]
    for file in expected_files:
        file_path = inputA_notemp / file
        assert file_path.exists()


@pytest.fixture
def inputFreyha():
    outdir = Path("tests/inputFreyja")
    outdir.mkdir(exist_ok=True)
    sys_argv = [
        "--outdir",
        str(outdir),
        "--threads",
        "4",
        "tests/example_fastqs/OUG-1858.subsample.1.fastq.gz",
        "tests/example_fastqs/OUG-1858.subsample.2.fastq.gz",
    ]
    main(sys_argv)
    yield outdir

    if outdir.exists():
        shutil.rmtree(outdir)


@pytest.fixture
def inputAssemblies():
    outdir = Path("tests/inputAssemblies")
    sys_argv = [
        "--outdir",
        str(outdir),
        "--assemblies",
        "--threads",
        "4",
        "tests/example_fasta/C10.fasta",
    ]
    main(sys_argv)
    yield outdir

    if outdir.exists():
        shutil.rmtree(outdir)


def test_assemblies_pipeline(inputAssemblies):
    results = inputAssemblies / "lineage_report.csv"
    assert results.exists()

    lines = results.read_text().splitlines()
    header = lines[0].split(",")
    rows = [dict(zip(header, line.split(","))) for line in lines[1:]]

    # The 34 contigs in C10.fasta all belong to a single sample, so classifying
    # with --assemblies should collapse them into one row rather than 34.
    assert len(rows) == 1
    assert rows[0]["sample_id"] == "C10"
    assert rows[0]["lineage"] == "T13"


@pytest.mark.skip(reason="Input files are too large")
def test_Freyja_pipeline(inputFreyha):
    results = inputFreyha / "lineage_report.csv"

    expected_result = (
        "sequence_id,lineage,confidence,classification_notes\n"
        "OUG-1858,T13,1.000,Freyja results: T13(100.0%)\n"
    )

    assert results.exists()
    assert results.read_text() == expected_result
