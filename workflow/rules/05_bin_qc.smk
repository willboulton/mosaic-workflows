rule vacuum_blobology:
rule flye_reassembly:
rule tidy_bins_with_assembly_graph:
rule bandage_view_assembly:
rule eukcc:
rule busco:
rule checkm2:
    input:
    output:
    conda:
        "checkm2"
    shell:
        """
        checkm2 predict \\
            --threads {threads} \\
            --input {input} \\
            --output-directory {output} \\
            --database_path $HOME/scratch/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd \\
            --force -x {params.extension}
        """
rule checkv:
rule gvclass:
