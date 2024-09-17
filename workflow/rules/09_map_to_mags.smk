rule strobealign_index:
    input:
        fna=f""
    output:
        f"{input.fna}.r150.sti"
    shell:
    """
    strobealign --create-index -t {threads} {input.fna} -r 150 
    """


rule strobealign_map_reads_to_mags:
    input:
        ref=f"",
        fq=f"", 
        idx=f"{input.ref}.r150.sti"
    output:
        f"{mapping_dir}/{{sample}}.strobealign_aemb.tsv"
    conda: "anvio"
    shell:
        """
        strobealign --use-index -t {threads} --aemb {input.ref} --interleaved {input.fq} > {output}
        """
