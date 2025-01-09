import sys
from pathlib import Path
from typing import Tuple

from vibecheck.src.utilities import run_command
from vibecheck.src.console import console


def downsample_reads(
    read1: Path, read2: Path, tempdir: Path, fraction: float = 0.2
) -> Tuple[Path, Path]:
    """Downsamples the input reads using `seqtk`.

    Parameters
    ----------
    read1 : Path
        Location of fastq.gz file contain first set of reads.
    read2 : Path
        Location of fastq.gz file contain second set of reads.
    tempdir : Path
        Location of temporary directory.
    fraction : float, optional
        Fraction of reads to sample (default is 0.2).

    Returns
    -------
    Tuple[Path, Path]
        location of fastq.gz files containing downsampled reads.
    """
    sub_read1 = tempdir / "subsampled_R1.fastq"
    sub_read2 = tempdir / "subsampled_R2.fastq"

    command1 = f"seqtk sample -s 42 {read1} {fraction} > {sub_read1}"
    command2 = f"seqtk sample -s 42 {read2} {fraction} > {sub_read2}"

    run_command(command1, error_message="Subsampling read1 failed")
    run_command(command2, error_message="Subsampling read2 failed")

    if not sub_read1.exists():
        console.log(
            f"Subsampled first set of reads {sub_read1} does not exist. Check {tempdir}"
        )
        sys.exit(-60)
    if not sub_read2.exists():
        console.log(
            f"Subsampled second set of reads {sub_read2} does not exist. Check {tempdir}"
        )
        sys.exit(-61)

    return sub_read1, sub_read2


def align_reads(reference: Path, read1: Path, read2: Path, tempdir: Path) -> Path:
    """Aligns reads to a reference using `minimap2` and processes with `samtools`.

    Parameters
    ----------
    reference : Path
        Location of FASTA file containing the reference sequence.
    read1 : Path
        Location of fastq.gz file contain first set of reads.
    read2 : Path
        Location of fastq.gz file contain second set of reads.
    tempdir : Path
        Location of temporary directory.

    Returns
    -------
    Path
        Location of sorted and indexed BAM file containing aligned reads.
    """
    alignment = tempdir / "alignment.bam"

    minimap_command = f"minimap2 -ax sr {reference} {read1} {read2} | samtools view -b - | samtools sort -o {alignment} -"
    index_command = f"samtools index {alignment}"

    run_command(minimap_command, error_message="Alignment of raw reads failed")
    run_command(index_command, error_message="Indexing alignment failed")

    if not alignment.exists():
        console.log(f"Error: Alignment {alignment} does not exist. Check {tempdir}.")
        sys.exit(-62)

    return alignment


def generate_depth(reference: Path, alignment: Path, tempdir: Path) -> Path:
    """Generates depth information from the alignment file using `samtools mpileup`.

    Parameters
    ----------
    reference : Path
        Location of FASTA file containing reference sequence.
    alignment : Path
        Location of sorted and indexed BAM file containing aligned reads.
    tempdir : Path
        Location of temporary directory.

    Returns
    -------
    Path
        Location of TSV file containing depth information.
    """
    depth = tempdir / "depth.txt"

    mpileup_command = f"samtools mpileup -aa -A -d 600000 -Q 20 -q 0 -B -f {reference} {alignment} | cut -f1-4 > {depth}"
    run_command(mpileup_command, error_message="Generating depth file failed")

    if not depth.exists():
        console.log(f"Depth file {depth} does not exist. Check {tempdir}.")
        console.exit(-63)
    return depth


def call_variants(
    reference: Path, alignment: Path, tempdir: Path, threads: int = 4
) -> Path:
    """Calls variants from the alignment file using bcftools.

    Parameters
    ----------
    reference : Path
        Location of FASTA file containing reference sequence.
    alignment : Path
        Location of BAM file containing alignment of reads.
    tempdir : Path
        Location of temporary directory.
    threads : int, optional
        Number of threads to use for variant calling (default is 4).

    Returns
    -------
    Path
        Location of VCF file containing variants.
    """
    variants = tempdir / "variants.vcf"
    variants_filled = tempdir / "variants_filled.vcf"

    mpileup_command = (
        f"bcftools mpileup --threads {threads} -d 600000 -Q 20 -q 0 -B -a INFO/AD,INFO/ADF,INFO/ADR -Ou "
        f"-f {reference} {alignment} | bcftools call --threads {threads} -mv -Ov --ploidy 1 -o {variants}"
    )
    fill_tags_command = (
        f"bcftools +fill-tags {variants} -Ou -o {variants_filled} -- -t AF"
    )

    run_command(mpileup_command, error_message="Variant calling failed")
    run_command(fill_tags_command, error_message="Filling tags in variants failed")

    if not variants_filled.exists():
        console.log(f"Variants file {variants_filled} does not exist. Check {tempdir}.")
        sys.exit(-65)

    return variants_filled


def freyja_demix(
    variants_filled: Path, depth: Path, barcodes: Path, tempdir: Path
) -> Path:
    """Resolves lineages from sample using `freyja demix`.

    Parameters
    ----------
    variants_filled : Path
        Location of VCF file containing variants.
    depth : Path
        Location of TSV file containing depth information.
    barcodes : Path
        Location of feather file containing lineage barcodes.
    tempdir : Path
        Location of temporary directory.

    Returns
    -------
    Path
        Location of text file containing results from freyja.
    """
    freyja_output = tempdir / "freyja_results.txt"
    demix_command = f"freyja demix {variants_filled} {depth} --output {freyja_output} --barcodes {barcodes}"
    run_command(demix_command, error_message="Freyja demix failed")

    if not freyja_output.exists():
        console.log(f"Freyja output {freyja_output} does not exist. Check {tempdir}.")
        sys.exit(-66)
    return freyja_output
