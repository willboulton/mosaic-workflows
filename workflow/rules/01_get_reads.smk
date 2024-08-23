def weighted_sum(counter_aggregation, node):
    
    # remove root node if it's not the only thing there
    # as this only drags estimates to be less specific
    
    if 1 in counter_aggregation:
        root_count = counter_aggregation.pop(1)
    
    if len(counter_aggregation) == 0:
        return root_count * node.get_distance(taxid_dict[1])
    
    return sum((w * node.get_distance(taxid_dict[n]) for n, w in counter_aggregation.items()))

def minimise_weighted_sum(counter_aggregation):
    
    ans = {}
    
    def agg_weights(counter_aggregation, node):

        ans[node.name] = weighted_sum(counter_aggregation, node)

        if node.up:
            agg_weights(counter_aggregation, node.up)
    
    for n in counter_aggregation.keys():
        
        ans = {}
        agg_weights(counter_aggregation, taxid_dict[n])
    
    #print(ans)
    min_dist = min([d for n,d in ans.items()])
    #print("min dist is", min_dist)
    
    candidate_nodes = [taxid_dict[n] for n, d in ans.items() if d == min_dist]
    #print(candidate_nodes)
    
    if len(candidate_nodes) == 1:
        return candidate_nodes[0]
    
    lineages = set()
    for c in candidate_nodes:
        x = set(lineage(c))
        x.remove(c)
        lineages = lineages.union(x)

    #print("lineages:", lineages)

    cn = []
    for c in candidate_nodes:
        if not c in lineages:
            cn.append(c)

    #print(cn)

    if len(cn) == 1:
        return cn[0]

    return root.get_common_ancestor(cn)

def load_ncbi_taxonomy(taxdump_dir):
    taxid_to_parent = {}
    taxid_to_rank = {}
    
    # Read the nodes.dmp file
    with open(f"{taxdump_dir}/nodes.dmp") as nodes_file:
        reader = csv.reader(nodes_file, delimiter='|', quoting=csv.QUOTE_MINIMAL)
        for row in reader:
            taxid = int(row[0].strip())
            parent_taxid = int(row[1].strip())
            rank = row[2].strip()
            taxid_to_parent[taxid] = parent_taxid
            taxid_to_rank[taxid] = rank
    
    return taxid_to_parent, taxid_to_rank


def get_lineage(taxid, taxid_to_parent):
    lineage = set()
    current_taxid = taxid
    while current_taxid != 1 and current_taxid in taxid_to_parent:
        lineage.add(current_taxid)
        current_taxid = taxid_to_parent[current_taxid]
    return lineage



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

rule bwa_mem:
    input:
        tar = f"{img_folder}/{{sample}}.tar.gz",
        fq = f"{fastq_tmpdir}/{{batch}}/{{sample}}-filter-METAG.fastq.gz"
    output:
        f"{bam_tmpdir}/{{batch}}/{{sample}}.sorted.bam"
    envmodules:
        "bwa/0.7.17",
        "samtools/1.16.1"
    shell:
        """
        cd {bam_tmpdir}/{wildcards.batch}
        zcat {input.tar} | \\
            tar xfO - {wildcards.sample}/{wildcards.sample}.a.fna > {wildcards.sample}.a.fna
        bwa index -a bwtsw {wildcards.sample}.a.fna
        bwa mem -M -t {threads} \\
            {wildcards.sample}.a.fna {input.fq} | \\
        samtools sort --threads {threads} -o {output} 
        """


rule split_fastq:
    input:
        fna = f"{binning_dir}/{{sample}}/{{sample}}.a.fna",
        bam = f"{bam_tmpdir}/{{sample}}/bam/{{sample}}.sorted.bam",
        unmapped = f"{taxo_tmpdir}/{{sample}}/{{sample}}.euks.cf",
        mapped = f"{taxo_tmpdir}/{{sample}}/{{sample}}.a.euks.txt"
    output:
        mapped = f"{fastq_tmpdir}/{{batch}}/{{sample}}-EVU-mapped.fastq.gz",
        unmapped = f"{fastq_tmpdir}/{{batch}}/{{sample}}-EVU-unmapped.fastq.gz"
    envmodules:
        "bbmap/37.28", 
        "samtools/1.16.1"
    conda:
        "blobology"
    shell:
        """
        samtools fastq -@ {threads} -f 12 {input.bam} | \\
        filterbyname.sh \\
            include=true \\
            in=stdin.fq.gz \\
            names={input.unmapped} \\
            out={output.unmapped}

        seqkit grep -n {input.fna} > {taxo_tmpdir}/{wildcards.sample}/{wildcards.sample}.a.euks.fna
        reformat.sh \\
            in={input.bam} \\
            ref={taxo_tmpdir}/{wildcards.sample}/{wildcards.sample}.a.euks.fna \\
            out={output.mapped} \\
            filterbits=12
        """

rule get_unmapped_reads:
    input:
        f"{bam_tmpdir}/{{sample}}/bam/{{sample}}.sorted.bam"
    output:
        fq1=f"{fastq_tmpdir}/{{sample}}-unmapped_1.fastq.gz",
        fq2=f"{fastq_tmpdir}/{{sample}}-unmapped_2.fastq.gz"
    envmodules:
        "samtools/1.16.1"
    shell:
        # -f 12 indicates both the read and the mate are unmapped. 
        """
        samtools fastq -@ {threads} -f 12 \\
            -1 {output.fq1} \\
            -2 {output.fq2} \\
            -0 /dev/null -s /dev/null -n {input}
        """


rule centrifuge_on_reads:
    input:
        fq1=f"{fastq_tmpdir}/{{sample}}-unmapped_1.fastq.gz",
        fq2=f"{fastq_tmpdir}/{{sample}}-unmapped_2.fastq.gz",
    output:
        all = f"{taxo_tmpdir}/{{sample}}/{{sample}}.unmapped.fastq.cf",
        evu = f"{taxo_tmpdir}/{{sample}}/{{sample}}.euks.cf"
    envmodules:
        "centrifuge/1.0.4"
    conda:
        "bio"
    shell:
        """
        centrifuge -x {software.centrifuge_idx} \\
            -1 {input.fq1} -2 {input.fq2} \\
            -q -p {threads} -S {output.all}

        {software.workflow}/scripts/centrifuge_reads_filter_EVU.py \\
            {output.all} \\
            {software.taxdump} \\
            {output.evu}
        """

rule split_fastq:
    input:
        fna = f"{binning_dir}/{{sample}}/{{sample}}.a.fna",
        bam = f"{bam_tmpdir}/{{sample}}/bam/{{sample}}.sorted.bam",
        mapped = f"{taxo_tmpdir}/{{sample}}/{{sample}}.a.euks.txt"
    output:
        mapped = f"{fastq_tmpdir}/{{batch}}/{{sample}}-EVU-mapped.fastq.gz"
    envmodules:
        "bbmap/37.28", 
        "samtools/1.16.1"
    conda:
        "blobology"
    shell:
        """
        seqkit grep -f {input.mapped} \\
            -o {taxo_tmpdir}/{wildcards.sample}/{wildcards.sample}.a.euks.fna \\
            {input.fna} 
            
        reformat.sh \\
            in={input.bam} \\
            ref={taxo_tmpdir}/{wildcards.sample}/{wildcards.sample}.a.euks.fna \\
            out={output.mapped} \\
            filterbits=12
        """

rule concat_fasta:
    input: unpack(collate_fastas)
    output: f"{binning_dir}/{{batch}}/{{batch}}.euks.fna"
    conda: "blobology"
    shell:
        """
        touch {output}
        paste \\
            <(echo "{input.euks}" | tr ' ' '\\n') \\
            <(echo "{input.fna}" | tr ' ' '\\n') \\
            <(echo "{input.map}" | tr ' ' '\\n') > {wildcards.batch}_concat_fastq_tmp
        
        while read -r euk fna map; do
            sample=$(basename $(dirname $euk))
            seqkit grep -f <(grep -xf $euk $map | cut -f 2) $fna | seqkit replace -p ^ -r ${{sample}}_ >> {output}
        done < {wildcards.batch}_concat_fastq_tmp
        
        rm -f {wildcards.batch}_concat_fastq_tmp
        if [ ! -s {output} ]; then
            rm -f {output}
            exit 1
        fi
        """


rule get_img_file:

rule get_coverage_file:

rule append_coverage_file:

rule comebin_single_assembly:
    # old preprocessing 
    shell:
        """
        {software.workflow}/scripts/contig_to_scaffold_ids.py \\
        {input.fna} \\
        {input.map} \\
        {binning_dir}/{wildcards.sample}/{wildcards.sample}.a.renamed.fna
        
        cd {binning_dir}/{wildcards.sample}
        python {software.comebin_src}/COMEBin/scripts/Filter_tooshort.py {wildcards.sample}.a.renamed.fna 1000
        """


rule kallisto_index:
    input:
        f"{binning_dir}/{{batch}}/{{batch}}.euks.fna"
    output:
        f"{binning_dir}/{{batch}}/{{batch}}.euks.fna.kindex"
    shell:
        """
        {software.kallisto}/kallisto index --distinguish -i {output} -t {threads} {input}
        """

rule deinterleave:
    input:
        f"{fastq_tmpdir}/{{sample}}/{{sample}}-filter-METAG.fastq.gz",
    output:
        fq1=f"{fastq_tmpdir}/{{sample}}/{{sample}}-filter-METAG_1.fastq.gz",
        fq2=f"{fastq_tmpdir}/{{sample}}/{{sample}}-filter-METAG_2.fastq.gz",
    conda:
        "bio"
    shell:
        """
        cd {fastq_tmpdir}/{wildcards.batch}
        {software.utils}/deinterleave_parallel {input}
        """


rule kallisto_x_sample:
    input:
        kindex = f"{binning_dir}/{{batch}}/{{batch}}.euks.fna.kindex", 
        fq1 = f"{fastq_tmpdir}/{{batch}}/{{sample}}-filter-METAG_1.fastq.gz",
        fq2 = f"{fastq_tmpdir}/{{batch}}/{{sample}}-filter-METAG_2.fastq.gz"
    output:
        directory(f"{binning_dir}/{{batch}}/{{sample}}_x_{{batch}}_euks_kallisto")
    shell:
        """
        {software.kallisto}/kallisto quant \\
            -i {input.kindex} \\
            -o {output} \\
            -t {threads} \\
            --plaintext {input.fq1} {input.fq2}
        """

# strobealign -t {threads} --aemb {input.ref} --interleaved {input.fq} > {output}

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


rule comebin_single_assembly:
    input:
        fna = f"{binning_dir}/{{sample}}/{{sample}}.a.renamed_1000.fna",
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
        awk '{{print $2 "\\t" $1}}' {input.map} > {wildcards.sample}_kv &&
        seqkit replace -p '(.+)$' -r '{{kv}}' -k {wildcards.sample}_kv {input.fna} | \\
        seqkit seq -m 1001 > {binning_dir}/{wildcards.sample}/{wildcards.sample}.a.renamed_1000.fa
        rm -f {wildcards.sample}_kv

        cd {software.comebin_src}/COMEBin
        bash run_comebin.sh \\
            -a {input.fna} \\
            -p {input.bam_dir} \\
            -o {output} \\
            -t {threads} 
        """    