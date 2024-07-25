rule extract_files:
    input:
        f"{binning_dir}/{{sample}}/{{sample}}.tar.gz"
    output:
        fna = "{sample}/{sample}.a.fna",
        depth = "{sample}/{sample}.a.depth.txt",
        phylodist = "{sample}/{sample}.a.phylodist.txt"
    shell:
        """
        tar -xzf {input} > {output.fna}
        tar -xzf {input} > {output.phylodist}
        tar -xzf {input} > {output.depth}
        """

rule vamb_single_assembly:
    input:
        fna=f"{binning_dir}/" + "{sample}/{sample}.a.fna",
        bam="/gpfs/home/qcf21tku/scratch/{sample}/bam/{sample}.sorted.bam"
    output:
        directory(f"{binning_dir}/" + "{sample}/vamb")
    conda:
        "anvio"
    envmodules:
        "cuda/11.8"
    shell:
        """
        vamb \
            --minfasta 50000 \\
            --norefcheck \\
            -p {threads} \
            -m 2000 \\
            --cuda \\
            --bamfiles {input.bam} \\
            --fasta {input.fna} 
        """


rule vamb_coassembly:
rule vamb_multiassembly:
rule metabat_single_assembly:
    input:
        img_file = f"{img_folder}/" + "{sample}.tar.gz"
    output:
        directory(f"{binning_dir}/" + "{sample}/metabat")
    params:
        fna = "{sample}/{sample}.a.fna",
        depth = "{sample}/{sample}.a.depth.txt"
    conda:
        "anvio"
    shell:
        """
        cd {binning_dir}

        metabat2 \\
            -i {params.fna} \\
            -a {params.depth} --cvExt \\
            -m 2000 \\
            -o {wildcards.sample}/metabat/mbat \\
            -t {threads} -v
        """

rule metabat_coassembly:
rule vacuum_single_assembly:
rule vacuum_coassembly:
rule vacuum_multiassembly:

rule comebin_single_assembly:
    input:
        tar = "",
        map = "",
        bam = "",
    output:
        directory(f"{binning_dir}/" + "{sample_id}/comebin")
    conda:
        "comebin_env"
    shell:
        """
        cd /gpfs/data/mock_lab/will/mosaic_tmp/binning_dir/{sample}
        python {comebin_dir}/scripts/Filter_tooshort.py {sample}.a.renamed.fna 1000
        cd {comebin_dir}
        scripts/contig_to_scaffold_ids.py {input.fna} {input.map} ${{input.fna}$.}.renamed.fna
        comebin \\
            -a ${{input.fna}$.}.renamed_1000.fa \\
            -p {input.bam} \\
            -o {output} \\
            -t {threads} 
        """

rule comebin_coassembly:
rule comebin_multiassembly:
rule metawrap:
rule euk_bin_size_filtering:
rule virus_sorting:
rule virsorter2:
