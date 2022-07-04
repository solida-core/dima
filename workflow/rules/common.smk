#######################
import errno
import pandas as pd
import os
import multiprocessing
import psutil
from snakemake.utils import validate


report: "../report/workflow.rst"

validate(config, schema="../schemas/config.schema.yaml")


samples = pd.read_table(config.get("samples"), index_col="sample")
units = pd.read_table(config.get("units"), index_col=["unit"], dtype=str)

def resolve_results_filepath(basepath, outname):
    return os.path.join(basepath, outname)

def expand_filepath(filepath):
    filepath = os.path.expandvars(os.path.expanduser(filepath))
    if not os.path.isabs(filepath):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT)+" (path must be absolute)", filepath)
    return filepath

def get_fastq(wildcards,units):
    if units.loc[wildcards.unit,["fq2"]].isna().all():
        return expand_filepath(units.loc[wildcards.unit,["fq1"]].dropna()[0])
    else:
        return expand_filepath(units.loc[wildcards.unit,["fq1"]].dropna()[0]),expand_filepath(units.loc[wildcards.unit,["fq2"]].dropna()[0])

def get_a_fastq(wildcards,units, fq="fq1"):
        return expand_filepath(units.loc[wildcards.unit,[fq]].dropna()[0])

def get_r2_fastq(wildcards,units):
    if units.loc[wildcards.unit,["fq2"]].isna().all():
        return expand_filepath(units.loc[wildcards.unit,["fq1"]].dropna()[0])
    else:
        return expand_filepath(units.loc[wildcards.unit,["fq1"]].dropna()[0]),expand_filepath(units.loc[wildcards.unit,["fq2"]].dropna()[0])




def get_trimmed_reads(wildcards,units):
    if units.loc[wildcards.unit,["fq2"]].isna().all():
        return rules.post_rename_fastq_se.output.r1
    else:
        return rules.post_rename_fastq_pe.output.r1,rules.post_rename_fastq_pe.output.r2

def total_physical_mem_size():
    mem = psutil.virtual_memory()
    return mem.total


def cpu_count():
    return multiprocessing.cpu_count()


def conservative_cpu_count(reserve_cores=1, max_cores=5):
    cores = max_cores if cpu_count() > max_cores else cpu_count()
    return max(cores - reserve_cores, 1)





def expand_filepath(filepath):
    filepath = os.path.expandvars(os.path.expanduser(filepath))
    if not os.path.isabs(filepath):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT)+" (path must be absolute)", filepath)
    return filepath


def resolve_single_filepath(basepath, filename):
    return os.path.join(basepath, filename)



def resolve_multi_filepath(basepath, dictionary):
    for k, v in dictionary.items():
        dictionary[k] = os.path.join(basepath, v)
    return dictionary


def tmp_path(path=''):
    """
    if does not exists, create path and return it. If any errors, return
    default path
    :param path: path
    :return: path
    """
    default_path = os.getenv('TMPDIR', config.get("paths").get("tmp_dir"))
    if path:
        try:
            os.makedirs(path)
        except OSError as e:
            if e.errno != errno.EEXIST:
                return default_path
        return path
    return default_path


def java_params(tmp_dir='', percentage_to_preserve=20, stock_mem=1024 ** 3,
                stock_cpu=2, multiply_by=1):
    """
    Set Java params
    :param tmp_dir: path to tmpdir
    :param percentage_to_preserve: percentage of resources to preserve
    :param stock_mem: min memory to preserve
    :param stock_cpu: min cpu to preserve
    :param multiply_by: multiply base resource by this param
    :return: string to return to configure java environments
    """

    def bytes2human(n):
        # http://code.activestate.com/recipes/578019
        # >>> bytes2human(10000)
        # '9.8K'
        # >>> bytes2human(100001221)
        # '95.4M'
        symbols = ('K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y')
        prefix = {}
        for i, s in enumerate(symbols):
            prefix[s] = 1 << (i + 1) * 10
        for s in reversed(symbols):
            if n >= prefix[s]:
                value = float(n) / prefix[s]
                return '%.0f%s' % (value, s)
        return "%sB" % n

    def preserve(resource, percentage, stock):
        preserved = resource - max(resource * percentage // 100, stock)
        return preserved if preserved != 0 else stock

    # def preserve(resource, percentage, stock):
    #     return resource - max(resource * percentage // 100, stock)

    params_template = "'-Xms{} -Xmx{} -XX:ParallelGCThreads={} " \
                      "-Djava.io.tmpdir={}'"

    mem_min = 1024 ** 3 * 2  # 2GB

    mem_size = preserve(total_physical_mem_size(), percentage_to_preserve,
                        stock_mem)
    cpu_nums = preserve(cpu_count(), percentage_to_preserve, stock_cpu)
    tmpdir = tmp_path(tmp_dir)

    return params_template.format(bytes2human(mem_min).lower(),
                                  bytes2human(max(mem_size//cpu_nums*multiply_by,
                                                  mem_min)).lower(),
                                  min(cpu_nums, multiply_by),
                                  tmpdir)


def get_units_by_sample(wildcards, samples, label='units', prefix=resolve_results_filepath(config.get("paths").get("results_dir"),'reads/sorted/'),
                        suffix='_sorted.cram'):
    return [prefix+i+suffix for i in samples.loc[wildcards.sample,
                                                  [label]][0].split(',')]

def get_odp(wildcards,samples,optical_dup='odp'):
    return "OPTICAL_DUPLICATE_PIXEL_DISTANCE={}".format(samples.loc[wildcards.sample, [optical_dup]].dropna()[0])

def get_known_sites(known_sites=['dbsnp','mills','ph1_indel']):
    known_variants = config.get("resources").get("known_variants")
    ks = []
    if len(known_sites) == 0:
        known_sites = known_variants.keys()
    for k, v in known_variants.items():
        if k in known_sites:
            ks.append("--known-sites {} ".format(v))
    return "".join(ks)