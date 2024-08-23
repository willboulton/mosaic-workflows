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
        directory(f"{bam_tmpdir}/{{batch}}/metabat_coa")
    conda:
        "anvio"
    shell:
        """
        cd {bam_tmpdir}
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
        directory(f"{bam_tmpdir}/{{batch}}/metabat_multi")
    conda:
        "anvio"
    shell:
        """
        cd {bam_tmpdir}
        metabat2 \\
            -i {input.fna} \\
            -a {input.depth} \\
            -m 2000 \\
            -o {wildcards.batch}/metabat_multi/mbat_multi \\
            -t {threads} -v
        """

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
    input:
        fna = f"{assembly_dir}/{{batch}}/megahit/final.contigs.fa",
        cov = f"{binning_dir}/{{batch}}/{{batch}}.euks_coa.jgi_cov.npz"
    output:
        bins = directory(f"{bam_tmpdir}/{{batch}}/vamb_coa"), 
        embeddings = f"{bam_tmpdir}/{{batch}}/vamb_coa/embeddings.npz"
    conda:
        "anvio"
    envmodules:
        "cuda/11.8"
    shell:
        """
        mkdir -p {output.bins}
        vamb \\
            --minfasta 200000 \\
            --norefcheck \\
            -p {threads} \
            -m 2000 \\
            --cuda \\
            --rpkm {input.cov} \\
            --fasta {input.fna} \\
            --outdir {output.bins}
        """

rule vamb_multiassembly:
    input:
        fna = f"{binning_dir}/{{batch}}/{{batch}}.euks.fna",
        cov = f"{binning_dir}/{{batch}}/{{batch}}.euks_multi.jgi_cov.npz"
    output:
        bins = directory(f"{bam_tmpdir}/{{batch}}/vamb_multi"), 
        embeddings = f"{bam_tmpdir}/{{batch}}/vamb_multi/embeddings.npz"
    conda:
        "anvio"
    envmodules:
        "cuda/11.8"
    shell:
        """
        mkdir -p {output.bins}
        vamb \\
            --minfasta 200000 \\
            --norefcheck \\
            -p {threads} \
            -m 2000 \\
            --cuda \\
            --rpkm {input.cov} \\
            --fasta {input.fna} \\
            --outdir {output.bins}
        """


rule vacuum_coassembly:
rule vacuum_multiassembly:

rule comebin_coassembly:
rule comebin_multiassembly:

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