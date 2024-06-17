rule globus:
    output:
        f"{fastq_tmpdir}" +  "/{sample}-filter-METAG.fastq.gz",
    conda:
        "bio"
    params:
        nersc_fastq_path = set_param('nersc_fastq_path')
    resources:
        globus_workers = 1
    shell:
        """
        TASK_ID=`globus transfer --skip-activation-check \
        {globus.nersc_uuid}:{params.nersc_fastq_path} \
        {globus.ada_uuid}:{params.img_taxon_id}-filter-METAG.fastq.gz`
        T=`echo $TASK_ID | cut -f 3 -d ':'`
        FINISHED=""
        while [ ! "$FINISHED" ]; do
            $FINISHED=`globus task show $T | grep Status | grep SUCCEEDED`
            sleep 60
        done
        mv {home}/{params.img_taxon_id}-filter-METAG.fastq.gz {output}
        """

rule get_img_file:

rule get_coverage_file:

rule append_coverage_file:
