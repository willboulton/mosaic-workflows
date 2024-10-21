def get_extension_by_binner(wc):
    extn_map = {
        "comebin":"fa",
        "metabat":"fa",
        "vacuum":"fna",
        "vamb":"fna", 
        "semibin":"fa"
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


# run aggregate_taxonomy per bin
# disregard checkm2 bins
# filter by size

#rule identify_superkingdom:

#rule sort_bins:


rule eukcc:
    input:
        f"{binning_dir}/{{batch}}/{{binner}}_{{coa|multi}}"
    output:
        f"{binning_dir}/{{batch}}/{{binner}}_{{coa|multi}}_eukcc"
    