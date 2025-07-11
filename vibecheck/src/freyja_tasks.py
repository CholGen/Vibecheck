import sys
from pathlib import Path
from typing import Tuple

import numpy as np

from vibecheck.src.utilities import run_command
from vibecheck.src.console import console


def run_pipeline(
    reads: Tuple[Path, Path],
    reference: Path,
    barcodes: Path,
    lineage_aliases: dict[str, str],
    subsample_fraction: float,
    no_subsample: bool,
    outfile: Path,
    tempdir: Path,
    threads: int,
) -> None:

    name = str(reads[0].name).split(".")[0]

    if no_subsample:
        console.log("Aligning reads to reference")
        alignment = align_reads(reads[0], reads[1], reference, tempdir, threads)
    else:
        console.log(f"Sampling {subsample_fraction:.0%} of reads for classification")

        if len( reads ) < 2:
            sub_read1, sub_read2 = downsample_reads( read1=reads[0], read2=None, tempdir=tempdir, fraction=subsample_fraction)
        else:
            sub_read1, sub_read2 = downsample_reads(*reads, tempdir=tempdir, fraction=subsample_fraction)

        console.log("Aligning reads to reference")
        alignment = align_reads(sub_read1, sub_read2, reference, tempdir, threads)

    console.log("Calculating depth of coverage across reference genome.")
    depth = generate_depth(alignment, reference, tempdir)

    console.log("Calling variants between reads and reference.")
    variants_filled = call_variants(alignment, reference, tempdir, threads)

    console.log("Calculating relative lineage abundances using Freyja.")
    freyja_results = freyja_demix(variants_filled, depth, barcodes, tempdir)

    console.log("Parsing Freyja results.")
    parse_freyja_results(freyja_results, name, outfile, lineage_aliases=lineage_aliases)


def downsample_reads(
    read1: Path, read2: Path, tempdir: Path, fraction: float = 0.2
) -> tuple[Path, Path] | tuple[Path, None]:
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
    Tuple[Path, Path] | Tuple[Path, None]
        location of fastq.gz files containing downsampled reads.
    """
    sub_read1 = tempdir / "subsampled_R1.fastq"
    command1 = f"seqtk sample -s 42 {read1} {fraction} > {sub_read1}"
    run_command(command1, error_message="Subsampling read1 failed")

    if not sub_read1.exists():
        console.log(
            f"Subsampled first set of reads {sub_read1} does not exist. Check {tempdir}"
        )
        sys.exit(-60)

    if read2 is not None:
        sub_read2 = tempdir / "subsampled_R2.fastq"
        command2 = f"seqtk sample -s 42 {read2} {fraction} > {sub_read2}"
        run_command(command2, error_message="Subsampling read2 failed")

        if not sub_read2.exists():
            console.log(
                f"Subsampled second set of reads {sub_read2} does not exist. Check {tempdir}"
            )
            sys.exit(-61)

        return sub_read1, sub_read2
    else:
        return sub_read1, None


def align_reads(
    read1: Path, read2: Path | None, reference: Path, tempdir: Path, threads: int
) -> Path:
    """Aligns reads to a reference using `minimap2` and processes with `samtools`.

    Parameters
    ----------
    read1 : Path
        Location of fastq.gz file contain first set of reads.
    read2 : Path | None
        Location of fastq.gz file contain second set of reads.
    reference : Path
        Location of FASTA file containing the reference sequence.
    tempdir : Path
        Location of temporary directory.
    threads : int
        number of threads to use during alignment

    Returns
    -------
    Path
        Location of sorted and indexed BAM file containing aligned reads.
    """
    alignment = tempdir / "alignment.bam"

    # Construct read inputs based on whether paired-end or single-end
    read_inputs = f"{read1} {read2}" if read2 is not None else str(read1)

    # Build the minimap2 command and associated BAM file wrangling.
    minimap_command = (
        f"minimap2 -ax sr -t {threads} {reference} {read_inputs} | "
        f"samtools view -b - | "
        f"samtools sort -o {alignment} -"
    )
    index_command = f"samtools index {alignment}"

    run_command(minimap_command, error_message="Alignment of raw reads failed")
    run_command(index_command, error_message="Indexing alignment failed")

    if not alignment.exists():
        console.log(f"Error: Alignment {alignment} does not exist. Check {tempdir}.")
        sys.exit(-62)

    return alignment


def generate_depth(alignment: Path, reference: Path, tempdir: Path) -> Path:
    """Generates depth information from the alignment file using `samtools mpileup`.

    Parameters
    ----------
    alignment : Path
        Location of sorted and indexed BAM file containing aligned reads.
    reference : Path
        Location of FASTA file containing reference sequence.
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
    alignment: Path, reference: Path, tempdir: Path, threads: int = 4
) -> Path:
    """Calls variants from the alignment file using bcftools.

    Parameters
    ----------
    alignment : Path
        Location of BAM file containing alignment of reads.
    reference : Path
        Location of FASTA file containing reference sequence.
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


def parse_freyja_results(freyja_results: Path, name: str, outfile: Path, lineage_aliases: dict[str, str] ) -> None:
    """Parses Freyja lineage abundance results from a text file. Exits if required
    fields are missing.

    Parameters
    ----------
    freyja_results : Path
        Location of text file containing Freyja results
    name :  str
        Name of sample.
    outfile : Path
        Location to save CSV file containing parsed Freyja results.
    lineage_aliases : dict[str, str]
        Dictionary of aliases mapped to lineage names.
    """
    # Define expected fields and their types
    field_parsers = {
        "lineages": lambda x: x[1:],  # Take everything after 'lineages'
        "abundances": lambda x: list(map(float, x[1:])),  # Convert to float
        "resid": lambda x: float(x[-1]),  # Take last value as float
        "coverage": lambda x: float(x[-1]),
    }

    results = {}

    with open(freyja_results) as f:
        for line in f:
            # Split line into tokens and find matching parser
            tokens = line.strip().split()
            field = next((k for k in field_parsers if line.startswith(k)), None)

            if field:
                results[field] = field_parsers[field](tokens)

    # Verify all required fields were found
    missing_fields = set(field_parsers) - set(results)
    if missing_fields:
        console.log(
            f"Error: The following fields are missing from the Freyja results: {missing_fields}. Please check the output of Freyja {freyja_results}"
        )
        sys.exit(-67)

    new_lineages = []
    for lin in results["lineages"]:
        if lin in lineage_aliases:
            new_lineages.append(lineage_aliases[lin])
        else:
            new_lineages.append(lin)
    results["lineages"] = new_lineages

    # Calculate conflict score, top lineage, and construct a summary
    confidence = np.exp(sum(a * np.log(a) for a in results["abundances"] if a > 0))
    top_lineage = results["lineages"][np.argmax(results["abundances"])]
    summary = "Freyja results: " + " ".join(
        f"{lin}({abun:.1%})"
        for lin, abun in zip(results["lineages"], results["abundances"])
    )

    with open(outfile, "wt") as out_file:
        out_file.write("sequence_id,lineage,confidence,classification_notes\n")
        out_file.write(f"{name},{top_lineage},{confidence:.3f},{summary}\n")
