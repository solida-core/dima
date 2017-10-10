from subprocess import run

if len(snakemake.input) > 1:
    subprocess.run(['samtools', 'merge', snakemake.output, snakemake.input])
else:
     subprocess.run(['cp', snakemake.input, snakemake.output])
     subprocess.run(['touch', '-h', snakemake.output])
