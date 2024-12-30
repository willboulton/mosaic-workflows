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
    input:
        fna = f"{assembly_dir}/{{batch}}/megahit/final.contigs.fa",
        depth = f"{binning_dir}/{{batch}}/{{batch}}.euks_coa.jgi_cov.tsv"
    output:
        directory(f"{binning_dir}/{{batch}}/metabat_coa")
    conda:
        "anvio"
    shell:
        """
        cd {binning_dir}
        metabat2 \\
            -i {input.fna} \\
            -a {input.depth} \\
            -m 2000 \\
            -o {wildcards.batch}/metabat_coa/mbat_coa \\
            -t {threads} -v
        """

rule metabat_multiassembly:
    input:
        fna = f"{binning_dir}/{{batch}}/{{batch}}.euks.fna",
        depth = f"{binning_dir}/{{batch}}/{{batch}}.euks_multi.jgi_cov.tsv"
    output:
        directory(f"{binning_dir}/{{batch}}/metabat_multi")
    conda:
        "anvio"
    shell:
        """
        cd {binning_dir}
        metabat2 \\
            -i {input.fna} \\
            -a {input.depth} \\
            -m 2000 \\
            -o {wildcards.batch}/metabat_multi/mbat_multi \\
            -t {threads} -v
        """

rule semibin_single_assembly:
    input:
        fna = f"{binning_dir}/{{sample}}/{{sample}}.a.fna",
        bam = f"{bam_tmpdir}/{{sample}}/bam/{{sample}}.sorted.bam"
    output:
        directory(f"{binning_dir}/{{sample}}/semibin")
    conda:
        "SemiBin"
    envmodules:
        "hmmer/3.3", "samtools/1.16.1", "bedtools/2.29.2"
    shell:
        """
        SemiBin2 single_easy_bin \\
            -i {input.fna} \\
            -o {output} \\
            -b {input.bam} -t {threads} -m 4000 --environment ocean
        """



rule vamb_single_assembly:
    input:
        fna = f"{binning_dir}/{{sample}}/{{sample}}.a.fna",
        bam = f"{bam_tmpdir}/{{sample}}/bam/{{sample}}.sorted.bam"
    output:
        directory(f"{binning_dir}/{{sample}}/vamb")
    conda:
        "anvio"
    envmodules:
        "cuda/11.8"
    shell:
        """
        vamb \\
            --minfasta 50000 \\
            --norefcheck \\
            -p {threads} \\
            -m 2000 \\
            --bamfiles {input.bam} \\
            --fasta {input.fna} \\
            --outdir {output}
        """



rule vamb_coassembly:
    input:
        fna = f"{assembly_dir}/{{batch}}/megahit/final.contigs.fa",
        cov = f"{binning_dir}/{{batch}}/{{batch}}.euks_coa.jgi_cov.npz"
    output:
        bins = directory(f"{binning_dir}/{{batch}}/vamb_coa"), 
        embeddings = f"{binning_dir}/{{batch}}/vamb_coa/latent.npz",
        contignames = f"{binning_dir}/{{batch}}/vamb_coa/contignames"
    conda:
        "anvio"
    envmodules:
        "cuda/11.8"
    shell:
        """
        cd {binning_dir}/{wildcards.batch}
        rm -rf vamb_coa && vamb \\
            --minfasta 200000 \\
            --norefcheck \\
            -p {threads} \\
            -m 2000 \\
            --rpkm {input.cov} \\
            --fasta {input.fna} \\
            --outdir vamb_coa
        """

rule vamb_multiassembly:
    input:
        fna = f"{binning_dir}/{{batch}}/{{batch}}.euks.fna",
        cov = f"{binning_dir}/{{batch}}/{{batch}}.euks_multi.jgi_cov.npz"
    output:
        bins = directory(f"{binning_dir}/{{batch}}/vamb_multi"), 
        embeddings = f"{binning_dir}/{{batch}}/vamb_multi/latent.npz",
        contignames = f"{binning_dir}/{{batch}}/vamb_multi/contignames"
    conda:
        "anvio"
    envmodules:
        "cuda/11.8"
    shell:
        """
        cd {binning_dir}/{wildcards.batch}
        rm -rf vamb_multi && vamb \\
            --minfasta 200000 \\
            --norefcheck \\
            -p {threads} \\
            -m 2000 \\
            --rpkm {input.cov} \\
            --fasta {input.fna} \\
            --outdir vamb_multi
        """


rule vacuum_coassembly:
    input:
        embeddings = f"{binning_dir}/{{batch}}/vamb_coa/latent.npz",
        contignames = f"{binning_dir}/{{batch}}/vamb_coa/contignames",
        fna = f"{assembly_dir}/{{batch}}/megahit/final.contigs.fa"
    output:
        fna_out = directory(f"{binning_dir}/{{batch}}/vamb_coa/vacuum_bins"),
        umap_out = f"{binning_dir}/{{batch}}/vamb_coa/umap_vamb.tsv", 
        mcl_out = f"{binning_dir}/{{batch}}/vamb_coa/vacuum_bins.tsv"
    envmodules:
        "mcl/14-137"
    conda: 
        "blobology"
    shell:
        """
        mkdir -p {output.fna_out}
        {software.workflow}/scripts/vacuum_blobology.py \\
            --embeddings {input.embeddings} \\
            --inflation 1.1 --embeddings_format VAMB \\
            --fasta_in {input.fna} --contignames {input.contignames} \\
            --umap_out {output.umap_out} --fasta_out {output.fna_out} \\
            --mcl_out {output.mcl_out}
        """


rule vacuum_multiassembly:
    input:
        embeddings = f"{binning_dir}/{{batch}}/vamb_multi/latent.npz",
        contignames = f"{binning_dir}/{{batch}}/vamb_multi/contignames",
        fna = f"{binning_dir}/{{batch}}/{{batch}}.euks.fna"
    output:
        fna_out = directory(f"{binning_dir}/{{batch}}/vamb_multi/vacuum_bins"),
        umap_out = f"{binning_dir}/{{batch}}/vamb_multi/umap_vamb.tsv", 
        mcl_out = f"{binning_dir}/{{batch}}/vamb_multi/vacuum_bins.tsv"
    envmodules:
        "mcl/14-137"
    conda: 
        "blobology"
    shell:
        """
        mkdir -p {output.fna_out}
        {software.workflow}/scripts/vacuum_blobology.py \\
            --embeddings {input.embeddings} \\
            --inflation 1.1 --embeddings_format VAMB \\
            --fasta_in {input.fna} --contignames {input.contignames} \\
            --umap_out {output.umap_out} --fasta_out {output.fna_out} \\
            --mcl_out {output.mcl_out}
        """


rule split_bins:
    input:
        f""
    output:
        directory("")
    conda:
        "blobology"
    shell:
        """
        seqkit grep -r -p {input} > {}
        """
