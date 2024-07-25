rule vacuum_blobology:
rule flye_reassembly:
rule tidy_bins_with_assembly_graph:
rule bandage_view_assembly:
rule identify_euks:
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
            --database_path {software.checkm_db} \\
            --force -x {params.extension}
        """
