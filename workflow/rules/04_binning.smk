#module comebin:
#    snakefile: f"{software.comebin_src}/Snakefile"

#use rule comebin from comebin as comebin_single_assembly with:
#    input:
#        fna = f"{binning_dir}/{{sample}}/{{sample}}.a.fna",
#        map = f"{binning_dir}/{{sample}}/{{sample}}.a.map.txt",
#        bam_dir = f"{bam_tmpdir}/{{sample}}/bam"
#    output:
#        directory(f"{binning_dir}/{{sample}}/comebin/comebin_res")    


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

rule comebin_coassembly:
rule comebin_multiassembly:
