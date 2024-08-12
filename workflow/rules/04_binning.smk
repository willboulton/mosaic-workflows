module comebin:
    snakefile: f"{software.comebin_src}/Snakefile"

use rule * from comebin


rule vamb_single_assembly:
    input:
        fna = f"{binning_dir}/{{sample}}/{{sample}}.a.fna",
        bam = f"{bam_tmpdir}/{{sample}}/bam/{{sample}}.sorted.bam"
    output:
        bins = directory(f"{binning_dir}/{{sample}}/vamb"), 
        embeddings = f"{binning_dir}/{{sample}}/vamb/embeddings.npz"
    conda:
        "anvio"
    envmodules:
        "cuda/11.8"
    shell:
        """
        vamb \\
            --minfasta 50000 \\
            --norefcheck \\
            -p {threads} \
            -m 2000 \\
            --cuda \\
            --bamfiles {input.bam} \\
            --fasta {input.fna} \\
            --outdir {output.bins}
        """


rule vamb_coassembly:
rule vamb_multiassembly:
rule metabat_single_assembly:
    input:
        fna = f"{binning_dir}/{{sample}}/{{sample}}.a.fna",
        depth = f"{binning_dir}/{{sample}}/{{sample}}.a.depth.txt"
    output:
        directory(f"{binning_dir}/{{sample}}/metabat")
    conda:
        "anvio"
    shell:
        """
        cd {binning_dir}
        metabat2 \\
            -i {input.fna} \\
            -a {input.depth} --cvExt \\
            -m 2000 \\
            -o {wildcards.sample}/metabat/mbat \\
            -t {threads} -v
        """

rule metabat_coassembly:
rule vacuum_single_assembly:
    input:
        fna = f"{binning_dir}/{{sample}}/{{sample}}.a.fna",
        embeddings = f"{binning_dir}/{{sample}}/vamb/embeddings.npz"
    output:
        bins = directory(f"{binning_dir}/{{sample}}/vacuum"),
        clusters = f"{binning_dir}/{{sample}}/vacuum/mcl_clusters_11.tsv",
        coords = f"{binning_dir}/{{sample}}/vacuum/umap_coords.tsv"
    conda:
        "blobology"
    envmodules:
        "mcl/14-137"
    shell:
        """
        {software.workflow}/scripts/vacuum_blobology.py \\
            --embeddings {input.embeddings} \\
            --embeddings_format VAMB \\
            --inflation 1.1 \\
            --fasta_in {input.fna} \\
            --fasta_out {output.bins} \\
            --umap_out {output.coords} \\
            --mcl_out {output.clusters} 
        """

rule vacuum_coassembly:
rule vacuum_multiassembly:

rule comebin_single_assembly:
    input:
        fna = f"{binning_dir}/{{sample}}/{{sample}}.a.fna",
        map = f"{binning_dir}/{{sample}}/{{sample}}.a.map.txt",
        bam = f"{bam_tmpdir}/{{sample}}/bam/{{sample}}.sorted.bam",
    output:
        directory(f"{binning_dir}/{{sample}}/comebin")
    conda:
        "comebin_env"
    shell:
        """
        mkfifo {wildcards.sample}_kv && \\
        awk '{print $2 "\t" $1}' {input.map} > {wildcards.sample}_kv &
        seqkit replace -p '(.+)$' -r '{kv}' -k {wildcards.samples}_kv {input.fna} | \\
        seqkit seq -m 1001 > {binning_dir}/{wildcards.sample}/{wildcards.sample}.a.renamed_1000.fa
        rm -f {wildcards.samples}_kv
        
        cd {software.comebin_src}/COMEBin
        bash run_comebin.sh \\
            -a {binning_dir}/{wildcards.sample}/{wildcards.sample}.a.renamed_1000.fa \\
            -p {input.bam} \\
            -o {output} \\
            -t {threads} 
        """

rule comebin_coassembly:
rule comebin_multiassembly:
