from functools import lru_cache

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

@lru_cache(maxsize = 64000)
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

rule deinterleave:
    input:
        f"{fastq_tmpdir}/{{batch}}/{{sample}}-filter-METAG.fastq.gz",
    output:
        fq1=f"{fastq_tmpdir}/{{batch}}/{{sample}}-filter-METAG_1.fastq.gz",
        fq2=f"{fastq_tmpdir}/{{batch}}/{{sample}}-filter-METAG_2.fastq.gz",
    conda:
        "bio"
    shell:
        """
        cd {fastq_tmpdir}/{wildcards.batch}
        {software.utils}/deinterleave_parallel {input}
        """

rule get_img_file:

rule get_coverage_file:

rule append_coverage_file:
