from pathlib import Path
from unittest.mock import patch

import pytest

from phonebook.src.tasks import (
    align_sequences,
    classify_usher,
    convert_to_vcf,
    run_command,
    usher_parsing,
)


def test_run_command_success(tmp_path):
    log_path = tmp_path / "log.txt"
    run_command("echo 'Hello, World!'", log_path)
    assert log_path.read_text().strip() == "Hello, World!"


def test_run_command_failure():
    with pytest.raises(SystemExit) as e:
        run_command("exit 1", error_message="Test failure")
    assert e.value.code == 1


def test_run_command_no_log():
    run_command("echo 'No log file'", None)  # Should not raise an error


# TODO: This basically just tests that the output is created and called
#  "alignment.fasta". That needs to be changed.
@patch("phonebook.src.tasks.run_command")
def test_align_sequences_success(mock_run_command, tmp_path):
    query_sequences = tmp_path / "query.fasta"
    reference = tmp_path / "reference.fasta"
    tempdir = tmp_path / "temp"
    alignment = tempdir / "alignment.fasta"

    query_sequences.touch()
    reference.touch()
    tempdir.mkdir()
    (tempdir / "logs").mkdir()
    alignment.touch()  # Simulate the output file creation

    result = align_sequences(query_sequences, reference, tempdir, threads=1)
    assert result == alignment


@patch("phonebook.src.tasks.run_command", side_effect=SystemExit(1))
def test_align_sequences_failure(mock_run_command, tmp_path):
    query_sequences = tmp_path / "query.fasta"
    reference = tmp_path / "reference.fasta"
    tempdir = tmp_path / "temp"

    query_sequences.touch()
    reference.touch()
    tempdir.mkdir()
    (tempdir / "logs").mkdir()

    with pytest.raises(SystemExit) as e:
        align_sequences(query_sequences, reference, tempdir, threads=1)
    assert e.value.code == 1


@patch("phonebook.src.tasks.run_command")
def test_convert_to_vcf_success(mock_run_command, tmp_path):
    aln = tmp_path / "alignment.fasta"
    reference = tmp_path / "reference.fasta"
    tempdir = tmp_path / "temp"
    vcf = tempdir / "sequences.aln.vcf"

    aln.touch()
    reference.touch()
    tempdir.mkdir()
    vcf.touch()  # Simulate VCF creation

    result = convert_to_vcf(aln, reference, tempdir)
    assert result == vcf


def test_convert_to_vcf_missing_alignment(tmp_path):
    reference = tmp_path / "reference.fasta"
    tempdir = tmp_path / "temp"
    alignment = tmp_path / "alignment.fasta"

    reference.touch()
    tempdir.mkdir()

    with pytest.raises(SystemExit) as e:
        convert_to_vcf(alignment, reference, tempdir)
    assert e.value.code != 0  # Could also be 1


@patch("phonebook.src.tasks.run_command")
def test_classify_usher_success(mock_run_command, tmp_path):
    vcf = tmp_path / "sequences.vcf"
    protobuf_tree = tmp_path / "tree.pb"
    tempdir = tmp_path / "temp"
    clades = tempdir / "clades.txt"

    vcf.touch()
    protobuf_tree.touch()
    tempdir.mkdir()
    clades.touch()  # Simulate clades file creation

    result = classify_usher(vcf, protobuf_tree, tempdir, threads=4)
    assert result == clades


def test_classify_usher_missing_vcf(tmp_path):
    protobuf_tree = tmp_path / "tree.pb"
    tempdir = tmp_path / "temp"

    protobuf_tree.touch()
    tempdir.mkdir()

    with pytest.raises(SystemExit) as e:
        classify_usher(Path("missing.vcf"), protobuf_tree, tempdir, threads=4)
    assert e.value.code != 0


def test_usher_parsing_valid_input(tmp_path):
    results = tmp_path / "results.txt"
    outfile = tmp_path / "parsed_results.csv"

    results.write_text("hash1\tA.28*|A.28(1/10),B.1(6/10)\nhash2\tB.1.1\n")
    usher_parsing(results, outfile)

    expected = (
        "sequence_id,lineage,conflict,usher_note\n"
        "hash1,A.28,0.9,Usher placements: A.28(1/10) B.1(6/10)\n"
        "hash2,B.1.1,,\n"
    )
    assert outfile.read_text() == expected


def test_usher_parsing_empty_input(tmp_path):
    results = tmp_path / "results.txt"
    outfile = tmp_path / "parsed_results.csv"

    results.touch()
    usher_parsing(results, outfile)

    assert outfile.read_text() == "sequence_id,lineage,conflict,usher_note\n"


def test_usher_parsing_malformed_input(tmp_path):
    results = tmp_path / "results.txt"
    outfile = tmp_path / "parsed_results.csv"

    results.write_text("invalid_line_without_tab\n")
    with pytest.raises(SystemExit) as e:  # Or the expected exception
        usher_parsing(results, outfile)
    assert e.value.code == -12
