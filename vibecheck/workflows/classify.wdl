version 1.0

workflow classify_cholera {
    meta {
        description: "Rapidly assigns O1 Vibrio cholerae genome sequences to transmission lineages using phylogenetic placement."
        author: "Nathaniel L. Matteson"
        email: "nmatteson@bwh.harvard.edu"
        allowNestedInputs: true
    }
    parameter_meta {
        query: {
            description: "Genome sequence in FASTA format",
            patterns: ["*.fasta","*.fa"]
        }
        max_ambiguiety: {
            description: "Maximum number of ambiguous bases a sequence can have before its filtered from the analysis"
        }
    }

    call vibecheck
    output {
        File vibecheck_report                   = vibecheck.lineage_report
        String vibecheck_lineage                = vibecheck.top_lineage
        Float vibecheck_confidence              = vibecheck.confidence
        String vibecheck_classification_notes   = vibecheck.classification_notes
        String vibecheck_pass_qc                = vibecheck.pass_qc
        String vibecheck_qc_notes               = vibecheck.qc_notes
        String vibecheck_version                = vibecheck.version
    }
}

task vibecheck {
    meta {
        description: "Rapidly assigns O1 Vibrio cholerae genome sequences to transmission lineages using phylogenetic placement."
    }
    input {
        File query_fasta
        File? usher_tree
        Float? max_ambiquity=0.3
        String docker="watronfire/vibecheck:latest"
    }
    Int disk_size = 16
    command <<<
        set -ex

        vibecheck -v | tee VERSION

        vibecheck "~{query_fasta}" \
            --outdir . \
            ~{"--usher-tree " + usher_tree} \
            ~{"--max-ambiquity " + max_ambiquity}

        python3 <<CODE
        import csv
        with open( "lineage_report.csv", "rt" ) as csv_file:
            reader = csv.DictReader( csv_file )
            line = next( reader )
            for key in ["lineage", "confidence", "classification_notes", "qc_status", "qc_notes"]:
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
        String pass_qc              = read_string("QC_STATUS")
        String qc_notes             = read_string("QC_NOTES")
        String version              = read_string("VERSION")
    }
}
