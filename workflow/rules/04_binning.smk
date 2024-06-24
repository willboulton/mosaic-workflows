#include: "00_utils.smk"

rule vamb_single_assembly:
rule vamb_coassembly:
rule vamb_multiassembly:
rule metabat_single_assembly:
    input:
        img_file = f"{img_folder}" + "/{sample}.tar.gz"
    output:
        directory(f"{binning_dir}" + "/{sample}/metabat")
    params:
        fna = "{sample}/{sample}.a.fna",
        depth = "{sample}/{sample}.a.depth.txt"
    conda:
        "anvio"
    shell:
        """
        cd {binning_dir}
        tar xzf {input.img_file} {params.fna}
        tar xzf {input.img_file} {params.depth}
        sort -k 2 -t '_' -n {params.depth} > {params.depth}.sorted
        metabat2 \\
            -i {params.fna} \\
            -a {params.depth}.sorted --cvExt \\
            -o {wildcards.sample}/metabat/mbat \\
            -t {threads}
        rm -f {wildcards.sample}/{wildcards.sample}.*
        """

rule metabat_coassembly:
rule vacuum_single_assembly:
rule vacuum_coassembly:
rule vacuum_multiassembly:
rule semibin_single_assembly:
rule semibin_coassembly:
rule semibin_multiassembly:
rule metawrap:
rule euk_bin_size_filtering:
rule virus_sorting:
rule virsorter2:
