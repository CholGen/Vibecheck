import sys
from pathlib import Path
import re
from typing import Tuple
import pandas as pd
import numpy as np

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from vibecheck.src.console import console
from vibecheck.src.utilities import run_command

OTHER_IUPAC = {"r", "y", "s", "w", "k", "m", "d", "h", "b", "v"}
VALID_CHARACTERS = [{"a"}, {"c"}, {"g"}, {"t"}, {"n"}, OTHER_IUPAC, {"-"}, {"?"}]


def run_pipeline(
    query_file: Path,
    protobuf_tree: Path,
    reference: Path,
    max_ambiguity: float,
    tempdir: Path,
    outfile: Path,
    threads: int,
):
    console.log("Aligning sequences to reference")
    aln = align_sequences(query_file, reference, tempdir, threads)

    console.log(
        f"Filtering sequences with greater than {max_ambiguity:.0%} ambiguous bases"
    )
    qc_stats, filtered_aln = sequence_qc(aln, tempdir, max_ambiguity)

    console.log("Converting alignment to VCF")
    vcf = convert_to_vcf(filtered_aln, reference, tempdir)

    console.log("Placing sequences into global phylogeny")
    results = classify_usher(vcf, protobuf_tree, tempdir, threads)

    console.log("Parsing Usher results")
    parsed_results = usher_parsing(results, tempdir)

    console.log("Writing results")
    combine_results(parsed_results, qc_stats, outfile)


def align_sequences(
    query_sequences: Path, reference: Path, tempdir: Path, threads: int
) -> Path:
    """Aligns sequences against a reference sequence using minimap2 and gofasta.

    Parameters
    ----------
    query_sequences: Path
        Location of FASTA file containing sequences to align against reference.
    reference: Path
        Location of FASTA file containing reference sequence.
    tempdir: Path
        Location of temporary directory.
    threads: int
        Number of threads to use.

    Returns
    -------
    Path
        Location of multisequence alignment in FASTA format.
    """
    log_path = tempdir / "logs/minimap2.txt"
    samfile = tempdir / "mapped.sam"
    alignment = tempdir / "alignment.fasta"

    align_command = f"minimap2 -a -x asm20 --sam-hit-only --secondary=no --score-N=0 -T {threads} {reference} {query_sequences} -o {samfile}"
    msa_command = f"gofasta sam toMultiAlign -s {samfile} -t {threads} --reference {reference} > {alignment}"

    # Run alignment with minimap2
    run_command(align_command, log_path, "Alignment with minimap2 failed")

    # Run conversion to multiple sequence alignment
    run_command(msa_command, log_path, "Conversion of alignment with gofasta failed")

    if not alignment.exists():
        console.log(
            f"Error: Alignment failed to generate the output {alignment}. Check {tempdir}."
        )
        sys.exit(-10)

    return alignment


def calculate_ambiquity(record: SeqRecord) -> float:
    """Calculates the proportion of a sequence that is ambiguous characters. Here
    ambiguous characters are anything besides ATCG.

    Parameters
    ----------
    record: Bio.SeqRecord.SeqRecord

    Returns
    -------
    float
    """
    seq = record.seq.lower()
    seq_length = len(seq)
    counts = []

    for v in VALID_CHARACTERS:
        counts.append(sum(map(lambda x: seq.count(x), v)))
    invalid_nucleotides = seq_length - sum(counts)

    if invalid_nucleotides > 0:
        console.log(
            f"Warning: Invalid characters in sequence {record.id}. Might not be a valid nucleotide sequence."
        )

    return 1 - (sum(counts[:4]) / seq_length)


def sequence_qc(
    alignment: Path, tempdir: Path, max_ambiguity: float = 0.3
) -> Tuple[Path, Path]:
    """Loads the aligned fasta file, calculates the percentage of bases that are N,
    if it's greater than max_ambiguity, indicate that qc is failed, writes result to
    file, but not to output fasta.

    Parameters
    ----------
    alignment: Path
        Location of aligned fasta file.
    tempdir: Path
        Location of temporary directory.
    max_ambiguity: float
        Maximum number of ambiguous bases a seqeuence can have before it is removed.

    Returns
    ----------
    Path
        Location of CSV file containing QC results
    Path
        Location of fasta alignment file, containing only sequences that passed QC.
    """

    pass_qc = tempdir / "alignment.filtered.fasta"
    qc_status = tempdir / "seq_status.csv"

    with open(qc_status, "w") as fw, open(pass_qc, "w") as fw_pass:
        fw.write("sequence_id,qc_status,qc_notes\n")

        total_input = 0
        total_pass = 0

        for record in SeqIO.parse(alignment, "fasta"):
            total_input += 1

            proportion_N = calculate_ambiquity(record)
            if proportion_N > max_ambiguity:
                fw.write(f"{record.id},fail,Ambiguous_content:{proportion_N:.2%}\n")
            else:
                total_pass += 1

                fw.write(f"{record.id},pass,Ambiguous_content:{proportion_N:.2%}\n")
                fw_pass.write(f">{record.id}\n{record.seq}\n")

    return qc_status, pass_qc


def convert_to_vcf(aln: Path, reference: Path, tempdir: Path) -> Path:
    """Converts a multi-sequence alignment in FASTA format to VCF format.

    Parameters
    ----------
    aln: Path
        Location of FASTA file containing multi-sequence alignment.
    reference: Path
        Location of FASTA file containing reference sequence.
    tempdir: Path
        Location of temporary directory.

    Returns
    -------
    Path
        Location of VCF file containing alignment.
    """
    seqs_ref = tempdir / "sequences.withref.fa"
    vcf = tempdir / "sequences.aln.vcf"

    generate_ref_aln_command = (
        f"cat {reference} > {seqs_ref} && cat {aln} >> {seqs_ref}"
    )
    generate_vcf_command = f"faToVcf {seqs_ref} {vcf}"

    # Concatenate reference sequence and alignment
    run_command(
        generate_ref_aln_command,
        error_message="Concatenation of reference sequence and alignment failed",
    )

    # Convert alignment to VCF
    run_command(
        generate_vcf_command, error_message="Conversion of alignment to VCF failed"
    )

    if not vcf.exists():
        console.log(f"Error: VCF file {vcf} does not exist. Check {tempdir}.")
        sys.exit(-11)

    return vcf


def classify_usher(vcf: Path, protobuf_tree: Path, tempdir: Path, threads: int) -> Path:
    """Place sequences into annotated phylogenetic tree using Usher.

    Parameters
    ----------
    vcf: Path
        Location of VCF file containing query sequences to classify.
    protobuf_tree: Path
        Location of protobuf tree to place sequences into.
    tempdir: Path
        Location of temporary directory.
    threads: int
        Number of threads to use.

    Returns
    -------
    Path
        Location of TSV file containing lineage assignments for query sequences.
    """
    clades = tempdir / "clades.txt"
    usher_log = tempdir / "logs/usher.txt"

    usher_command = f"usher -n -D -i {protobuf_tree} -v {vcf} -T {threads} -d {tempdir}"

    run_command(
        usher_command,
        log_path=usher_log,
        error_message="Classification with Usher failed",
    )

    if not clades.exists():
        console.log(f"Error: Usher results {clades} does not exist. Check {tempdir}.")
        sys.exit(-12)

    return clades


def usher_parsing(results: Path, tempdir: Path) -> Path:
    """Parses results from Usher classification into columnar format. Basically just
    splits lineage histogram produced when there are multiple parsimonious placements.

    Parameters
    ----------
    results: Path
        Location of results from Usher classification.
    tempdir: Path
        Location of temporary directory.
    """

    outfile = tempdir / "parsed_results.csv"

    with open(outfile, "w") as fw:
        fw.write("sequence_id,lineage,confidence,classification_notes\n")

        with open(results, "r") as f:
            for line in f:

                try:
                    name, lineage_histogram = line.rstrip("\n").split("\t")
                except ValueError:
                    console.log(
                        f"Error: Usher result file is malformed. Please check the file {results}."
                    )
                    sys.exit(-12)
                if "*|" in lineage_histogram:
                    # example: A.28*|A.28(1/10),B.1(6/10),B.1.511(1/10),B.1.518(2/10)
                    lineage, histogram = lineage_histogram.split("*|")
                    histo_list = [i for i in histogram.split(",") if i]
                    confidence = 0.0
                    if len(histo_list) > 1:
                        for lin_counts in histo_list:
                            m = re.match(
                                r"([A-Z0-9.]+)\(([0-9]+)\/([0-9]+)\)", lin_counts
                            )
                            if m:
                                place_count, total = [
                                    int(m.group(2)),
                                    int(m.group(3)),
                                ]
                                confidence += (place_count / total) * np.log(
                                    place_count / total
                                )
                    confidence = np.exp(confidence)
                    histogram_note = "Usher placements: " + " ".join(histo_list)
                else:
                    lineage = lineage_histogram
                    confidence = 1.0
                    histogram_note = ""

                fw.write(f"{name},{lineage},{confidence},{histogram_note}\n")
    return outfile


def combine_results(usher_results: Path, qc_results: Path, outfile: Path) -> None:
    usher_pd = pd.read_csv(usher_results)
    qc_pd = pd.read_csv(qc_results)
    results = qc_pd.merge(usher_pd, how="outer", on="sequence_id")
    results.to_csv(outfile, index=False)
