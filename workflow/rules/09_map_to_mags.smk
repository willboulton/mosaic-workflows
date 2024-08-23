rule strobealign_index:
    input:
        fna=f""
    output:
        f".sti"
    shell:
    """
    strobealign --create-index -t {threads} {input.fna} -r 150 
    """


rule strobealign_map_reads_to_mags:
    input:
        ref=f"",
        fq=f"", 
        idx=f".sti"
    output:
    shell:
        """
        strobealign --use-index {input.idx} -t {threads} --aemb {input.ref} --interleaved {input.fq} > {output}
        """

