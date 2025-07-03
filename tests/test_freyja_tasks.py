import shutil

import pandas as pd
import pytest
from pathlib import Path
from unittest.mock import patch, call
from vibecheck.src.freyja_tasks import (
    downsample_reads,
    align_reads,
    generate_depth,
    call_variants,
    freyja_demix,
    parse_freyja_results,
)


@pytest.fixture
def mock_paths():
    return {
        "read1": Path("/test/read1.fastq.gz"),
        "read2": Path("/test/read2.fastq.gz"),
        "tempdir": Path("/test/temp"),
        "reference": Path("/test/reference.fasta"),
        "alignment": Path("/test/temp/alignment.bam"),
        "depth": Path("/test/temp/depth.txt"),
        "variants": Path("/test/temp/variants_filled.vcf"),
        "barcodes": Path("/test/barcodes.feather"),
    }


@pytest.fixture
def mock_file_exists():
    with patch("pathlib.Path.exists", return_value=True):
        yield


def test_freyja_pipeline_dependencies():
    required_tools = ["freyja", "bcftools", "samtools", "minimap2", "seqtk"]
    missing_tools = []

    for tool in required_tools:
        if not shutil.which(tool):
            missing_tools.append(tool)

    if missing_tools:
        pytest.fail(f"Required tools not found in PATH: {', '.join(missing_tools)}")


@patch("vibecheck.src.freyja_tasks.run_command")
def test_downsample_reads(
    mock_run_command,
    mock_paths,
    mock_file_exists,
):
    result = downsample_reads(
        mock_paths["read1"],
        mock_paths["read2"],
        mock_paths["tempdir"],
        fraction=0.2,
    )

    expected_calls = [
        call(
            "seqtk sample -s 42 /test/read1.fastq.gz 0.2 > /test/temp/subsampled_R1.fastq",
            error_message="Subsampling read1 failed",
        ),
        call(
            "seqtk sample -s 42 /test/read2.fastq.gz 0.2 > /test/temp/subsampled_R2.fastq",
            error_message="Subsampling read2 failed",
        ),
    ]
    mock_run_command.assert_has_calls(expected_calls)
    assert result == (
        mock_paths["tempdir"] / "subsampled_R1.fastq",
        mock_paths["tempdir"] / "subsampled_R2.fastq",
    )

@patch("vibecheck.src.freyja_tasks.run_command")
def test_downsample_single_read(
    mock_run_command,
    mock_paths,
    mock_file_exists,
):
    result = downsample_reads(
        mock_paths["read1"],
        None,
        mock_paths["tempdir"],
        fraction=0.2,
    )

    expected_calls = [
        call(
            "seqtk sample -s 42 /test/read1.fastq.gz 0.2 > /test/temp/subsampled_R1.fastq",
            error_message="Subsampling read1 failed",
        )
    ]
    mock_run_command.assert_has_calls(expected_calls)
    assert result == (
        mock_paths["tempdir"] / "subsampled_R1.fastq",
        None
    )


@patch("vibecheck.src.freyja_tasks.run_command")
def test_align_reads(mock_run_command, mock_paths, mock_file_exists):
    result = align_reads(
        mock_paths["read1"],
        mock_paths["read2"],
        mock_paths["reference"],
        mock_paths["tempdir"],
        1,
    )

    expected_calls = [
        call(
            "minimap2 -ax sr -t 1 /test/reference.fasta /test/read1.fastq.gz /test/read2.fastq.gz | samtools view -b - | samtools sort -o /test/temp/alignment.bam -",
            error_message="Alignment of raw reads failed",
        ),
        call(
            "samtools index /test/temp/alignment.bam",
            error_message="Indexing alignment failed",
        ),
    ]
    mock_run_command.assert_has_calls(expected_calls)
    assert result == mock_paths["tempdir"] / "alignment.bam"

@patch("vibecheck.src.freyja_tasks.run_command")
def test_align_reads_single_file(mock_run_command, mock_paths, mock_file_exists):
    result = align_reads(
        mock_paths["read1"],
        None,
        mock_paths["reference"],
        mock_paths["tempdir"],
        1
    )

    expected_calls = [
        call(
            "minimap2 -ax sr -t 1 /test/reference.fasta /test/read1.fastq.gz | samtools view -b - | samtools sort -o /test/temp/alignment.bam -",
            error_message="Alignment of raw reads failed",
        ),
        call(
            "samtools index /test/temp/alignment.bam",
            error_message="Indexing alignment failed",
        ),
    ]
    mock_run_command.assert_has_calls(expected_calls)
    assert result == mock_paths["tempdir"] / "alignment.bam"


@patch("vibecheck.src.freyja_tasks.run_command")
def test_generate_depth(mock_run_command, mock_paths, mock_file_exists):
    result = generate_depth(
        mock_paths["alignment"], mock_paths["reference"], mock_paths["tempdir"]
    )

    mock_run_command.assert_called_once_with(
        "samtools mpileup -aa -A -d 600000 -Q 20 -q 0 -B -f /test/reference.fasta /test/temp/alignment.bam | cut -f1-4 > /test/temp/depth.txt",
        error_message="Generating depth file failed",
    )
    assert result == mock_paths["tempdir"] / "depth.txt"


@patch("vibecheck.src.freyja_tasks.run_command")
def test_call_variants(mock_run_command, mock_paths, mock_file_exists):
    result = call_variants(
        mock_paths["alignment"], mock_paths["reference"], mock_paths["tempdir"]
    )

    expected_calls = [
        call(
            "bcftools mpileup --threads 4 -d 600000 -Q 20 -q 0 -B -a INFO/AD,INFO/ADF,INFO/ADR -Ou -f /test/reference.fasta /test/temp/alignment.bam | bcftools call --threads 4 -mv -Ov --ploidy 1 -o /test/temp/variants.vcf",
            error_message="Variant calling failed",
        ),
        call(
            "bcftools +fill-tags /test/temp/variants.vcf -Ou -o /test/temp/variants_filled.vcf -- -t AF",
            error_message="Filling tags in variants failed",
        ),
    ]
    mock_run_command.assert_has_calls(expected_calls)
    assert result == mock_paths["tempdir"] / "variants_filled.vcf"


@patch("vibecheck.src.freyja_tasks.run_command")
def test_freyja_demix(mock_run_command, mock_paths, mock_file_exists):
    result = freyja_demix(
        mock_paths["variants"],
        mock_paths["depth"],
        mock_paths["barcodes"],
        mock_paths["tempdir"],
    )

    mock_run_command.assert_called_once_with(
        "freyja demix /test/temp/variants_filled.vcf /test/temp/depth.txt --output /test/temp/freyja_results.txt --barcodes /test/barcodes.feather",
        error_message="Freyja demix failed",
    )
    assert result == mock_paths["tempdir"] / "freyja_results.txt"


@pytest.fixture
def freyja_resultsA(tmp_path):
    input_file = tmp_path / "freyja_results.txt"
    output_file = tmp_path / "parsed_results.csv"

    input_text = (
        "mixed_sample.vcf\n"
        "summarized      [('Other', 0.9999999968413253)]\n"
        "lineages        MEASLES-D9 MEASLES-H1\n"
        "abundances      0.79692605 0.20307394\n"
        "resid   214.51679168207156\n"
        "coverage        91.39927016484208\n"
    )

    input_file.write_text(input_text)
    parse_freyja_results(input_file, "testA", output_file)
    return output_file


def test_parse_freyja_results_A(freyja_resultsA):
    expected_result = (
        "sequence_id,lineage,confidence,classification_notes\n"
        "testA,MEASLES-D9,0.604,Freyja results: MEASLES-D9(79.7%) MEASLES-H1(20.3%)\n"
    )
    assert freyja_resultsA.read_text() == expected_result


def test_parse_freyja_results_valid_csv(freyja_resultsA):
    df = pd.read_csv(freyja_resultsA)
    assert df.shape == (1, 4)
    assert df.columns.to_list() == [
        "sequence_id",
        "lineage",
        "confidence",
        "classification_notes",
    ]


@pytest.fixture
def freyja_resultsB(tmp_path):
    input_file = tmp_path / "freyja_results.txt"
    output_file = tmp_path / "parsed_results.csv"

    input_text = (
        "\tOUG-1858.variants-filled.vcf\n"
        "summarized\t[('Other', 0.9999999999873125)]\n"
        "lineages\tT13\n"
        "abundances\t1.00000000\n"
        "resid\t2.661981530311349\n"
        "coverage\t94.37996916326536\n"
    )

    input_file.write_text(input_text)
    parse_freyja_results(input_file, "testB", output_file)
    return output_file


def test_parse_freyja_results_B(freyja_resultsB):
    expected_result = (
        "sequence_id,lineage,confidence,classification_notes\n"
        "testB,T13,1.000,Freyja results: T13(100.0%)\n"
    )
    assert freyja_resultsB.read_text() == expected_result

def test_parse_freyja_results_with_alias(tmp_path):
    input_file = tmp_path / "freyja_results.txt"
    output_file = tmp_path / "parsed_results.csv"

    input_text = (
        "\tOUG-1858.variants-filled.vcf\n"
        "summarized\t[('Other', 0.9999999999873125)]\n"
        "lineages\tT13\n"
        "abundances\t1.00000000\n"
        "resid\t2.661981530311349\n"
        "coverage\t94.37996916326536\n"
    )

    aliases = {"T13" : "Other"}

    input_file.write_text(input_text)
    parse_freyja_results(input_file, "testB", output_file, aliases=aliases)

    expected_result = (
        "sequence_id,lineage,confidence,classification_notes\n"
        "testB,Other,1.000,Freyja results: Other(100.0%)\n"
    )

    assert output_file.read_text() == expected_result


