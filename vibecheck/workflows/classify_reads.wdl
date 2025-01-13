version 1.0

workflow classify_cholera_reads {
    meta {
        description: "Rapidly assigns O1 Vibrio cholerae sequencing data to transmission lineages using variant frequency demixing."
        author: "Nathaniel L. Matteson"
        email: "nmatteson@bwh.harvard.edu"
        allowNestedInputs: true
    }
    parameter_meta {
        read1: {
            description: "Forward reads",
            patterns: ["*.fastq","*.fq", "*.fastq.gz", "*.fq.gz"]
        }
        read2: {
            description: "Reverse reads",
            patterns: ["*.fastq","*.fq", "*.fastq.gz", "*.fq.gz"]
        }
        subsampling_fraction: {
            description: "Fraction of reads to use in classification. Default: 0.2"
        }
        skip_subsampling: {
            description: "Whether to skip subsampling of reads."
        }
    }

    call vibecheck_reads
    output {
        File vibecheck_reads_report                 = vibecheck_reads.lineage_report
        String vibecheck_reads_lineage              = vibecheck_reads.top_lineage
        Float vibecheck_reads_confidence            = vibecheck_reads.confidence
        String vibecheck_reads_classification_notes = vibecheck_reads.classification_notes
        String vibecheck_reads_version              = vibecheck_reads.version
    }
}


task vibecheck_reads {
    meta {
        description: "Rapidly assigns O1 Vibrio cholerae sequencing data to transmission lineages using variant frequency demixing."
    }
    input {
        File read1
        File read2
        File? lineage_barcodes
        Float? subsampling_fraction
        Boolean skip_subsampling=false
        String docker="watronfire/vibecheck:latest"
    }
    Int disk_size = 16
    command <<<
        set -ex

        vibecheck -v | tee VERSION

        vibecheck "~{read1}" "~{read2}"  \
            --outdir . \
            ~{"--barcodes " + lineage_barcodes} \
            ~{"--subsampling_fraction " + subsampling_fraction} \
            ~{true='--no-detect' false='' skip_subsampling}

        python3 <<CODE
        import csv
        with open( "lineage_report.csv", "rt" ) as csv_file:
            reader = csv.DictReader( csv_file )
            line = next( reader )
            for key in ["lineage", "confidence", "classification_notes"]:
                with open( key.upper(), "wt" ) as outf:
                    outf.write(line[key])
        CODE
    >>>
    runtime {
        docker: docker
        memory: "3 GB"
        cpu:    2
        disks:  "local-disk " + disk_size + " HDD"
        disk: disk_size + " GB" # TES
        dx_instance_type: "mem1_ssd1_v2_x2"
    }
    output {
        File lineage_report         = "lineage_report.csv"
        String top_lineage          = read_string("LINEAGE")
        Float confidence            = read_float("CONFIDENCE")
        String classification_notes = read_string("CLASSIFICATION_NOTES")
        String version              = read_string("VERSION")
    }
}
