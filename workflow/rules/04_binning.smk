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
            --cuda \\
            --bamfiles {input.bam} \\
            --fasta {input.fna} \\
            --outdir {output}
        """



rule vamb_coassembly:
    input:
        fna = f"{assembly_dir}/{{batch}}/megahit/final.contigs.fa",
        cov = f"{binning_dir}/{{batch}}/{{batch}}.euks_coa.jgi_cov.npz"
    output:
        bins = directory(f"{bam_tmpdir}/{{batch}}/vamb_coa"), 
        embeddings = f"{bam_tmpdir}/{{batch}}/vamb_coa/latent.npz"
    conda:
        "anvio"
    envmodules:
        "cuda/11.8"
    shell:
        """
        cd {bam_tmpdir}/{wildcards.batch}
        rm -rf vamb_coa && vamb \\
            --minfasta 200000 \\
            --norefcheck \\
            -p {threads} \\
            -m 2000 \\
            --cuda \\
            --rpkm {input.cov} \\
            --fasta {input.fna} \\
            --outdir vamb_coa
        """

rule vamb_multiassembly:
    input:
        fna = f"{binning_dir}/{{batch}}/{{batch}}.euks.fna",
        cov = f"{binning_dir}/{{batch}}/{{batch}}.euks_multi.jgi_cov.npz"
    output:
        bins = directory(f"{bam_tmpdir}/{{batch}}/vamb_multi"), 
        embeddings = f"{bam_tmpdir}/{{batch}}/vamb_multi/latent.npz"
    conda:
        "anvio"
    envmodules:
        "cuda/11.8"
    shell:
        """
        cd {bam_tmpdir}/{wildcards.batch}
        rm -rf vamb_multi && vamb \\
            --minfasta 200000 \\
            --norefcheck \\
            -p {threads} \\
            -m 2000 \\
            --cuda \\
            --rpkm {input.cov} \\
            --fasta {input.fna} \\
            --outdir vamb_multi
        """


rule vacuum_coassembly:
    input:
        embeddings = f"{bam_tmpdir}/{{batch}}/vamb_multi/latent.npz",
        fna = f"{binning_dir}/{{batch}}/{{batch}}.euks.fna"
    output:
        
    conda: "blobology"
    shell:
        """
        {software.workflow}/scripts/vacuum_blobology.py \\
            
        """


rule vacuum_multiassembly:
    input:
    output:
    conda: "blobology"
    shell:
        """
        
        """


rule comebin_coassembly:
    input:
    output:
    conda:
        "comebin_env"
    shell:
        """
        """

rule comebin_multiassembly:
    input:
    output:
    conda:
        "comebin_env"
    shell:
        """
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