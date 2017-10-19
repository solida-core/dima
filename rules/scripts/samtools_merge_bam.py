from subprocess import run

if len(snakemake.input) > 1:
    run(['samtools', 'merge', snakemake.output[0], snakemake.input[0], snakemake.input[1]])
else:
     run(['cp', snakemake.input[0], snakemake.output[0]])
     run(['touch', '-h', snakemake.output[0]])
