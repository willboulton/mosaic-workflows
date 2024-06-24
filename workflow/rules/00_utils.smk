configfile: "config/config.yaml"

from collections import namedtuple
import pandas as pd

def initialise_config(config):

    """
    """
    software = namedtuple("software", (k for k in config["software"].keys()))(**config["software"])

    home = config['home']
    fastq_tmpdir = config['fastq_tmpdir']
    taxo_tmpdir = config['taxo_tmpdir']
    bam_tmpdir = config['bam_tmpdir']
    binning_dir = config['binning_dir']

    globus = namedtuple("globus", (k for k in config["globus"].keys()))(**config["globus"])

    samples = to_samplelist(pd.read_csv(config['samples'], sep = '\t'))
    run_batch = config['run_batch']
    img_folder = config['img_folder']
    return (
        software, 
        home, 
        fastq_tmpdir, 
        taxo_tmpdir, 
        bam_tmpdir, 
        binning_dir,
        globus, 
        samples, 
        run_batch, 
        img_folder
    ) 


def to_samplelist(df):
    """
    """
    STuple = namedtuple("STuple", df.columns)
    return [STuple(**d) for d in df.to_dict('records')]

def sample_dict(wc):
    result = [s for s in samples if str(s.img_taxon_id) == wc.sample][0]._asdict()
    return result

def set_param(varname):
    return (lambda wc: sample_dict(wc)[varname])

(
    software, 
    home, 
    fastq_tmpdir, 
    taxo_tmpdir,
    bam_tmpdir,
    binning_dir,
    globus,
    samples, 
    run_batch, 
    img_folder
) = initialise_config(config)
