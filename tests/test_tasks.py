from pathlib import Path
from unittest.mock import patch

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phonebook.src.tasks import (
    align_sequences,
    calculate_ambiquity,
    classify_usher,
    convert_to_vcf,
    run_command,
    sequence_qc,
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
    usher_parsing(results, tmp_path)

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
    usher_parsing(results, tmp_path)

    assert outfile.read_text() == "sequence_id,lineage,conflict,usher_note\n"


def test_usher_parsing_malformed_input(tmp_path):
    results = tmp_path / "results.txt"

    results.write_text("invalid_line_without_tab\n")
    with pytest.raises(SystemExit) as e:  # Or the expected exception
        usher_parsing(results, tmp_path)
    assert e.value.code == -12


def test_calculate_completeness_valid_sequence():
    record = SeqRecord(Seq("ACGTACGTACGT"), id="valid_seq")
    completeness = calculate_ambiquity(record)
    assert completeness == 0


def test_calculate_completeness_with_ambiguity():
    record = SeqRecord(Seq("ACGTNNNNACGT"), id="ambiguous_seq")
    completeness = calculate_ambiquity(record)
    assert completeness == pytest.approx(4 / 12)


@patch("rich.console.Console.log")
def test_calculate_completeness_with_invalid_characters(mock_log):
    record = SeqRecord(Seq("ACGTACGTXXX"), id="invalid_seq")
    completeness = calculate_ambiquity(record)
    assert mock_log.call_count == 1
    mock_log.assert_called_with(
        "Warning: Invalid characters in sequence invalid_seq. Might not be a valid nucleotide sequence."
    )
    # Invalid characters don't contribute to counts
    assert completeness == pytest.approx(3 / 11)


def create_fasta(records, path):
    """Helper function to create a test FASTA file."""
    with open(path, "w") as f:
        SeqIO.write(records, f, "fasta")


def test_sequence_qc_all_pass(tmp_path):
    alignment = tmp_path / "alignment.fasta"
    tempdir = tmp_path / "temp"
    tempdir.mkdir()

    # Create test sequences
    records = [
        SeqRecord(Seq("ACGTACGTACGT"), id="seq1"),
        SeqRecord(Seq("ACGTNNNNACGT"), id="seq2"),
    ]
    create_fasta(records, alignment)

    qc_status, pass_qc = sequence_qc(alignment, tempdir, max_ambiguity=0.5)

    # Check QC status file
    with open(qc_status) as f:
        lines = f.readlines()
        assert len(lines) == 3  # Header + 2 sequences
        assert "seq1,pass" in lines[1]
        assert "seq2,pass" in lines[2]

    # Check filtered FASTA file
    passed_records = list(SeqIO.parse(pass_qc, "fasta"))
    assert len(passed_records) == 2


def test_sequence_qc_some_fail(tmp_path):
    alignment = tmp_path / "alignment.fasta"
    tempdir = tmp_path / "temp"
    tempdir.mkdir()

    # Create test sequences
    records = [
        SeqRecord(Seq("ACGTACGTACGT"), id="seq1"),
        SeqRecord(Seq("NNNNNNNNNNNN"), id="seq2"),
    ]
    create_fasta(records, alignment)

    qc_status, pass_qc = sequence_qc(alignment, tempdir, max_ambiguity=0.5)

    # Check QC status file
    with open(qc_status) as f:
        lines = f.readlines()
        assert len(lines) == 3  # Header + 2 sequences
        assert "seq1,pass" in lines[1]
        assert "seq2,fail" in lines[2]

    # Check filtered FASTA file
    passed_records = list(SeqIO.parse(pass_qc, "fasta"))
    assert len(passed_records) == 1
    assert passed_records[0].id == "seq1"


def test_sequence_qc_all_fail(tmp_path):
    alignment = tmp_path / "alignment.fasta"
    tempdir = tmp_path / "temp"
    tempdir.mkdir()

    # Create test sequences
    records = [
        SeqRecord(Seq("NNNNNNNNNNNN"), id="seq1"),
        SeqRecord(Seq("NNNNNNNNNNNN"), id="seq2"),
    ]
    create_fasta(records, alignment)

    qc_status, pass_qc = sequence_qc(alignment, tempdir, max_ambiguity=0.1)

    # Check QC status file
    with open(qc_status) as f:
        lines = f.readlines()
        assert len(lines) == 3  # Header + 2 sequences
        assert "seq1,fail" in lines[1]
        assert "seq2,fail" in lines[2]

    # Check filtered FASTA file
    passed_records = list(SeqIO.parse(pass_qc, "fasta"))
    assert len(passed_records) == 0


def test_sequence_qc_empty_input(tmp_path):
    alignment = tmp_path / "empty.fasta"
    tempdir = tmp_path / "temp"
    tempdir.mkdir()

    # Create empty FASTA file
    alignment.touch()

    qc_status, pass_qc = sequence_qc(alignment, tempdir, max_ambiguity=0.5)

    # Check QC status file
    with open(qc_status) as f:
        lines = f.readlines()
        assert len(lines) == 1  # Only header

    # Check filtered FASTA file
    passed_records = list(SeqIO.parse(pass_qc, "fasta"))
    assert len(passed_records) == 0
