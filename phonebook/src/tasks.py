import sys
from pathlib import Path
import subprocess
import re
from typing import Tuple, Union

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from phonebook.src.console import console

OTHER_IUPAC = {"r", "y", "s", "w", "k", "m", "d", "h", "b", "v"}
VALID_CHARACTERS = [{"a"}, {"c"}, {"g"}, {"t"}, {"n"}, OTHER_IUPAC, {"-"}, {"?"}]


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
    seq = record.seq.lower()
    l = len(seq)
    counts = []

    for v in VALID_CHARACTERS:
        counts.append(sum(map(lambda x: seq.count(x), v)))
    invalid_nucleotides = l - sum(counts)

    if invalid_nucleotides > 0:
        console.log(
            f"Warning: Invalid characters in sequence {record.id}. Might not be a valid nucleotide sequence."
        )

    return 1 - (sum(counts[:4]) / l)


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
    """Converts a multisequence alignment in FASTA format to VCF format.

    Parameters
    ----------
    aln: Path
        Location of FASTA file containing multisequence alignment.
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
    generate_vcf_command = f"faToVcf -includeNoAltN {seqs_ref} {vcf}"

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


def usher_parsing(results: Path, outfile: Path) -> None:
    """Parses results from Usher classification into columnar format. Basically just
    splits lineage histogram produced when there are multiple parsimonious placements.

    Parameters
    ----------
    results: Path
        Location of results from Usher classification.
    outfile: Path
        Location of desired output file.
    """
    with open(outfile, "w") as fw:
        fw.write("sequence_id,lineage,conflict,usher_note\n")

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
                    conflict = 0.0
                    if len(histo_list) > 1:
                        selected_count = 0
                        total = 0
                        for lin_counts in histo_list:
                            m = re.match(
                                r"([A-Z0-9.]+)\(([0-9]+)/([0-9]+)\)", lin_counts
                            )
                            if m:
                                lin, place_count, total = [
                                    m.group(1),
                                    int(m.group(2)),
                                    int(m.group(3)),
                                ]
                                if lin == lineage:
                                    selected_count = place_count
                                    break
                        conflict = (total - selected_count) / total
                    histogram_note = "Usher placements: " + " ".join(histo_list)
                else:
                    lineage = lineage_histogram
                    conflict = ""
                    histogram_note = ""

                fw.write(f"{name},{lineage},{conflict},{histogram_note}\n")
