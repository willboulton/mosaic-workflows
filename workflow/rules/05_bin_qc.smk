def get_extension_by_binner(wc):
    extn_map = {
        "comebin":"fna",
        "metabat":"fna",
        "vacuum":"fna",
        "vamb":"fna"
    }
    return extn_map.get(wc.binner)

rule checkm2:
    input:
        f"{binning_dir}/{{sample}}/{{binner}}"
    output:
        directory(f"{binning_dir}/{{sample}}/{{binner}}_checkm")
    params:
        extension = get_extension_by_binner
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

# rule vacuum_blobology:
# rule identify_euks:
# rule eukcc:
# rule busco:

# rule flye_reassembly:
# rule tidy_bins_with_assembly_graph:
# rule bandage_view_assembly: