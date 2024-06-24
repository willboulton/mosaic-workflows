rule kallisto_map_reads_to_mags:
    input:
        pass 

rule kallisto_index:
    input:
        pass
    output:
        pass
    shell:
    """
    {software.kallisto}/kallisto index --index 3300060026.kindex --distinguish \\
    <(zcat 3300060026.tar.gz | tar xfO - 3300060026/3300060026.a.fna | head) 
    """

rule kallisto_quant:
    input:
        pass
    output:
        pass
    shell:
    """
    {software.kallisto}/kallisto quant \\
    --index {input.index} --output-dir {} \\
    --threads {threads} \\
    {input.fq1} {input.fq2}
    """
