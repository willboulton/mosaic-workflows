rule comebin_single_assembly:
    input:
        fna = f"{binning_dir}/{{sample}}/{{sample}}.a.fna",
        map = f"{binning_dir}/{{sample}}/{{sample}}.a.map.txt",
        bam = f"{bam_tmpdir}/{{sample}}/bam/{{sample}}.sorted.bam",
        bam_dir = f"{bam_tmpdir}/{{sample}}/bam"
    output:
        directory(f"{binning_dir}/{{sample}}/comebin_old")
    conda:
        "comebin_env"
    shell:
        """
        mkfifo {wildcards.sample}_kv && \\
        awk '{{print $2 "\\t" $1}}' {input.map} > {wildcards.sample}_kv &
        seqkit replace -p '(.+)$' -r '{{kv}}' -k {wildcards.sample}_kv {input.fna} | \\
        seqkit seq -m 1001 > {binning_dir}/{wildcards.sample}/{wildcards.sample}.a.renamed_1000.fa
        rm -f {wildcards.sample}_kv

        cd {software.comebin_src}/COMEBin
        bash run_comebin.sh \\
            -a {binning_dir}/{wildcards.sample}/{wildcards.sample}.a.renamed_1000.fa \\
            -p {input.bam_dir} \\
            -o {output} \\
            -t {threads} 
        """    

rule cbin_00_preprocessing:
    input:
        map = f"{binning_dir}/{{sample}}/{{sample}}.a.map.txt",
        fna = f"{binning_dir}/{{sample}}/{{sample}}.a.fna"
    output:
        f"{binning_dir}/{{sample}}/{{sample}}.a.renamed_1000.fa"
    conda: "comebin_env"
    shell:
        """
        mkfifo {wildcards.sample}_kv && \\
        awk '{print $2 "\t" $1}' {input.map} > {wildcards.sample}_kv & \\
        seqkit replace -p '(.+)$' -r '{kv}' -k {wildcards.sample}_kv {input.fna} | \\
        seqkit seq -m 1001 > {output}
        rm -f {wildcards.sample}_kv
        """

rule cbin_01_augmentation_data:
    input:
        contigs_file = f"{binning_dir}/{{sample}}/{{sample}}.a.renamed_1000.fa", 
        bam_dir = f"{bam_tmpdir}/{{sample}}/bam"
    output:
        directory(f"{binning_dir}/{{sample}}/comebin/augmentation_data")
    conda: "comebin_env"
    shell:
        """
        cd ~/software/COMEBin/COMEBin
        bash run_comebin_01.sh \\
            -a {input.contigs_file} \\
            -p {input.bam_dir} \\
            -o {binning_dir}/{wildcards.sample}/comebin \\
            -t {threads} 
        """


rule cbin_02_training:
    input:
        contigs_file = f"{binning_dir}/{{sample}}/{{sample}}.a.renamed_1000.fa", 
        bam_dir = f"{bam_tmpdir}/{{sample}}/bam",
        augmentation_data = f"{binning_dir}/{{sample}}/comebin/augmentation_data"
    output:
        directory(f"{binning_dir}/{{sample}}/comebin/comebin_res")
    conda: "comebin_env"
    shell:
        """
        cd ~/software/COMEBin/COMEBin
        bash run_comebin_02.sh \\
            -a {input.contigs_file} \\
            -p {input.bam_dir} \\
            -o {binning_dir}/{wildcards.sample}/comebin \\
            -t {threads} 
        """


rule cbin_03_leiden_clustering:
    input:
        contigs_file = f"{binning_dir}/{{sample}}/{{sample}}.a.renamed_1000.fa", 
        bam_dir = f"{bam_tmpdir}/{{sample}}/bam",
        res = f"{binning_dir}/{{sample}}/comebin/comebin_res"
    output:
        directory(f"{binning_dir}/{{sample}}/comebin")
    conda: "comebin_env"
    shell:
        """
        cd ~/software/COMEBin/COMEBin
        bash run_comebin_03.sh \\
            -a {input.contigs_file} \\
            -p {input.bam_dir} \\
            -o {binning_dir}/{wildcards.sample}/comebin \\
            -t {threads} 
        mv {output}/comebin_res_bins/*.fa {binning_dir}/{wildcards.sample}/comebin
        rm -f {contigs_file}.*
        """


rule cbin_04_bin_qc:
    input:
        contigs_file = f"{binning_dir}/{{sample}}/{{sample}}.a.renamed_1000.fa", 
        bam_dir = f"{bam_tmpdir}/{{sample}}/bam",
        clust = f"{binning_dir}/{{sample}}/comebin/comebin_clusters"
    output:
        directory(f"{binning_dir}/{{sample}}/comebin/comebin_out")
    conda: "comebin_env"
    shell:
        """
        cd ~/software/COMEBin/COMEBin
        bash run_comebin_04.sh \\
            -a {input.contigs_file} \\
            -p {input.bam_dir} \\
            -o {binning_dir}/{wildcards.sample}/comebin \\
            -t {threads} 
        """

