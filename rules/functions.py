import errno
import multiprocessing
import os.path
import psutil


def total_physical_mem_size():
    mem = psutil.virtual_memory()
    return mem.total


def cpu_count():
    return multiprocessing.cpu_count()


def tmp_path(path=''):
    """
    if does not exists, create path and return it. If any errors, return
    default path
    :param path: path
    :return: path
    """
    default_path = os.getenv('TMPDIR', '/tmp')
    if path:
        try:
            os.makedirs(path)
        except OSError as e:
            if e.errno != errno.EEXIST:
                return default_path
    return path


def conservative_cpu_count(reserve_cores=1, max_cores=10):
    cores = max_cores if cpu_count() > max_cores else cpu_count()
    return cores - reserve_cores


def java_params(tmp_dir='', stock_mem=1024 ** 2, stock_cpu=2, fraction_for=1):
    """
    Set Java params
    :param tmp_dir: path to tmpdir
    :param fraction_for: divide resource for this param
    :param stock_mem: default memory to stock
    :param stock_cpu: default cpu to reserve
    :return: string to be used to configure java
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

    params_template = "-Xms2g -Xmx{} -XX:ParallelGCThreads={} -Djava.io.tmpdir={}"

    reserve_mem = stock_mem//fraction_for
    mem_size = total_physical_mem_size() if reserve_mem > total_physical_mem_size() \
        else total_physical_mem_size() - reserve_mem

    reserve_cpu = stock_cpu//fraction_for
    cpu_nums = cpu_count() if reserve_cpu > cpu_count() \
        else cpu_count() - reserve_cpu

    tmpdir = tmp_path(tmp_dir)
    return params_template.format(bytes2human(mem_size//fraction_for).lower(),
                                  cpu_nums//fraction_for,
                                  tmpdir)


def references_abs_path():
    references = config.get('references')
    basepath = references['basepath']
    provider = references['provider']
    genome = references['genome_release']

    return [os.path.join(basepath, provider, genome)]


def resolve_single_filepath(basepath, filename):
    return [os.path.join(basepath, filename)]


def resolve_multi_filepath(basepath, dictionary):
    for k, v in dictionary.items():
        dictionary[k] = os.path.join(basepath, v)
    return dictionary
