from subprocess import run

if len(snakemake.input) > 1:
    run(['samtools', 'merge', snakemake.output, snakemake.input])
else:
     run(['cp', snakemake.input, snakemake.output])
     run(['touch', '-h', snakemake.output])
