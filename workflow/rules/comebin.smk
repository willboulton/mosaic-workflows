rule cbin_00_preprocessing:
    input:
        f"{binning_dir}/{{sample}}/{{sample}}.a.fna"
    output:
        f"{binning_dir}/{{sample}}/{{sample}}.a.renamed_1000.fa"
    conda: "comebin_env"
    shell:
        """
        seqkit seq -m 1001 {input} > {output}
        """

use rule cbin_00_preprocessing as cbin_00_multiassembly with:
    input:
        f"{binning_dir}/{{batch}}/{{batch}}.euks.fna"
    output:
        f"{binning_dir}/{{batch}}/{{batch}}.euks.renamed_1000.fna"

use rule cbin_00_preprocessing as cbin_00_coassembly with:
    input:
        f"{assembly_dir}/{{batch}}/megahit/final.contigs.fa"
    output:
        f"{assembly_dir}/{{batch}}/megahit/final.contigs.renamed_1000.fa"


rule cbin_01_data_augmentation:
    input:
        contigs_file = f"{binning_dir}/{{sample}}/{{sample}}.a.renamed_1000.fa", 
        bam_dir = f"{bam_tmpdir}/{{sample}}/bam"
    output:
        directory(f"{binning_dir}/{{sample}}/comebin/data_augmentation")
    conda: "comebin_env"
    shell:
        """
        cd ~/software/COMEBin/COMEBin
        out=$(dirname {output})
        bash run_comebin_01.sh \\
            -a {input.contigs_file} \\
            -p {input.bam_dir} \\
            -o $out \\
            -t {threads} 
        """

use rule cbin_01_data_augmentation as cbin_01_multiassembly with:
    input:
        contigs_file = f"{binning_dir}/{{batch}}/{{batch}}.euks.renamed_1000.fna", 
        bam_dir = f"{bam_tmpdir}/{{batch}}/bam_multi"
    output:
        directory(f"{binning_dir}/{{batch}}/comebin_multi/data_augmentation")

use rule cbin_01_data_augmentation as cbin_01_coassembly with:
    input:
        contigs_file = f"{assembly_dir}/{{batch}}/megahit/final.contigs.renamed_1000.fa", 
        bam_dir = f"{bam_tmpdir}/{{batch}}/bam_coa"
    output:
        directory(f"{binning_dir}/{{batch}}/comebin_coa/data_augmentation")

rule cbin_02_training:
    input:
        contigs_file = f"{binning_dir}/{{sample}}/{{sample}}.a.renamed_1000.fa", 
        bam_dir = f"{bam_tmpdir}/{{sample}}/bam",
        data_augmentation = f"{binning_dir}/{{sample}}/comebin/data_augmentation"
    output:
        directory(f"{binning_dir}/{{sample}}/comebin/comebin_res")
    conda: "comebin_env"
    shell:
        """
        cd ~/software/COMEBin/COMEBin
        out=$(dirname {output})
        bash run_comebin_02.sh \\
            -a {input.contigs_file} \\
            -p {input.bam_dir} \\
            -o $out \\
            -t {threads} 
        """

use rule cbin_02_training as cbin_02_multiassembly with:
    input:
        contigs_file = f"{binning_dir}/{{batch}}/{{batch}}.euks.renamed_1000.fna", 
        bam_dir = f"{bam_tmpdir}/{{batch}}/bam_multi",
        data_augmentation = f"{binning_dir}/{{batch}}/comebin_multi/data_augmentation"
    output:
        directory(f"{binning_dir}/{{batch}}/comebin_multi/comebin_res")

use rule cbin_02_training as cbin_02_coassembly with:
    input:
        contigs_file = f"{assembly_dir}/{{batch}}/megahit/final.contigs.renamed_1000.fa", 
        bam_dir = f"{bam_tmpdir}/{{batch}}/bam_coa",
        data_augmentation = f"{binning_dir}/{{batch}}/comebin_coa/data_augmentation"
    output:
        directory(f"{binning_dir}/{{batch}}/comebin_coa/comebin_res")

rule cbin_03_leiden_clustering:
    input:
        contigs_file = f"{binning_dir}/{{sample}}/{{sample}}.a.renamed_1000.fa", 
        bam_dir = f"{bam_tmpdir}/{{sample}}/bam",
        res = f"{binning_dir}/{{sample}}/comebin/comebin_res"
    output:
        directory(f"{binning_dir}/{{sample}}/comebin/comebin_clusters")
    conda: "comebin_env"
    shell:
        """
        cd ~/software/COMEBin/COMEBin
        out=$(dirname {output})
        bash run_comebin_03.sh \\
            -a {input.contigs_file} \\
            -p {input.bam_dir} \\
            -o $out \\
            -t {threads} 
        """

use rule cbin_03_leiden_clustering as cbin_03_multiassembly with:
    input:
        contigs_file = f"{binning_dir}/{{batch}}/{{batch}}.euks.renamed_1000.fna", 
        bam_dir = f"{bam_tmpdir}/{{batch}}/bam_multi",
        res = f"{binning_dir}/{{batch}}/comebin_multi/comebin_res"
    output:
        directory(f"{binning_dir}/{{batch}}/comebin_multi/comebin_clusters")

use rule cbin_03_leiden_clustering as cbin_03_coassembly with:
    input:
        contigs_file = f"{assembly_dir}/{{batch}}/megahit/final.contigs.renamed_1000.fa", 
        bam_dir = f"{bam_tmpdir}/{{batch}}/bam_coa",
        res = f"{binning_dir}/{{batch}}/comebin_coa/comebin_res"
    output:
        directory(f"{binning_dir}/{{batch}}/comebin_coa/comebin_clusters")

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
        out=$(dirname {output})
        bash run_comebin_04.sh \\
            -a {input.contigs_file} \\
            -p {input.bam_dir} \\
            -o $out \\
            -t {threads} 
        """

use rule cbin_04_bin_qc as cbin_04_multiassembly with:
    input:
        contigs_file = f"{binning_dir}/{{batch}}/{{batch}}.euks.renamed_1000.fna", 
        bam_dir = f"{bam_tmpdir}/{{batch}}/bam_multi",
        clust = f"{binning_dir}/{{batch}}/comebin_multi/comebin_clusters"
    output:
        directory(f"{binning_dir}/{{batch}}/comebin_multi/comebin_out")

use rule cbin_04_bin_qc as cbin_04_coassembly with:
    input:
        contigs_file = f"{assembly_dir}/{{batch}}/megahit/final.contigs.renamed_1000.fa", 
        bam_dir = f"{bam_tmpdir}/{{batch}}/bam_coa",
        clust = f"{binning_dir}/{{batch}}/comebin_coa/comebin_clusters"
    output:
        directory(f"{binning_dir}/{{batch}}/comebin_coa/comebin_out")