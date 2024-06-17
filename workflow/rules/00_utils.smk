from collections import namedtuple
import pandas as pd

def initialise_config(config):

    """
    """
    software = namedtuple("software", (k for k in config["software"].keys()))(**config["software"])

    home = config['home']
    fastq_tmpdir = config['fastq_tmpdir']

    globus = namedtuple("globus", (k for k in config["globus"].keys()))(**config["globus"])

    samples = to_samplelist(pd.read_csv(config['samples'], sep = '\t'))
    run_batch = config['run_batch']
    img_path = config['img_path']
    return (software, home, fastq_tmpdir, globus, samples, run_batch, img_path) 


def to_samplelist(df):
    """
    """
    STuple = namedtuple("STuple", df.columns)
    return [STuple(**d) for d in df.to_dict('records')]
