rule collate_bins:
    input:
        metabat=f"{binning_dir}/{{sample}}/metabat",
        semibin=f"{binning_dir}/{{sample}}/semibin",
        vamb=f"{binning_dir}/{{sample}}/vamb"
    conda:
        "blobology"
    output:
        directory(f"{binning_dir}/{{sample}}/bins"),
        directory(f"{binning_dir}/{{sample}}/ebins"),
        directory(f"{binning_dir}/{{sample}}/stats")
    shell:
        """
        {software.workflow}/scripts/move_bins.sh \\
            {binning_dir}/{wildcards.sample} \\
            {taxo_tmpdir}/{wildcards.sample}/{wildcards.sample}.a.tiara.txt
        """

rule das_tool:
    input:
        f"{binning_dir}/{{sample}}"
    output:
        directory("f{binning_dir}/{{sample}}/bins_das")
    envmodules:
        "R/4.2.1", "ruby/3.3.5"
    conda:
        "bio"
    shell:
        """
        {software.workflow}/scripts/run_dastool.sh {input}
        """


rule checkm2:
    input:
        f"{binning_dir}/{{sample}}/bins"
    output:
        directory(f"{binning_dir}/{{sample}}/bins_checkm")
    conda:
        "checkm2"
    shell:
        """
        checkm2 predict \\
            --threads {threads} \\
            --input {input} \\
            --output-directory {output} \\
            --database_path {software.checkm_db} \\
            --force -x fna
        """

rule eukcc_single_sample:
    input:
        f"{binning_dir}/{{sample}}/ebins"
    output:
        directory(f"{binning_dir}/{{sample}}/ebins_eukcc")
    envmodules:
        "singularity/3.6.0" 
    shell:
        """
        [ -z "$(ls {input})" ] && mkdir -p {output} && touch {output}/eukcc.csv && exit 0
        singularity exec {software.eukcc}/eukcc_latest.sif eukcc \\
            folder -t {threads} --db {software.eukcc}/eukcc2_db_ver_1.2 \\
            --suffix .fna --out {output} {input}
        """

# run aggregate_taxonomy per bin
# filter by size

rule eukcc:
    input:
        f"{binning_dir}/{{batch}}/{{binner}}_{{coa|multi}}"
    output:
        f"{binning_dir}/{{batch}}/{{binner}}_{{coa|multi}}_eukcc"
    envmodules:
        "singularity/3.6.0" 
    shell:
        """
        singularity exec {software.eukcc}/eukcc_latest.sif eukcc \\
            folder -t {threads} --db {software.eukcc}/eukcc2_db_ver_1.2 \\
            --suffix .fna --out {output} {input}
        """