from subprocess import run

if len(snakemake.input) > 1:
    cmd = [snakemake.params['cmd'], 'merge',
           '--threads', str(snakemake.threads),
           '-O', snakemake.params['output_fmt']]
    if 'genome' in snakemake.params:
        cmd.append('--reference')
        cmd.append(snakemake.params['genome'][0])
    cmd.append(snakemake.output[0])
    for i in snakemake.input:
        cmd.append(i)
    run(cmd)
else:
     run(['cp', snakemake.input[0], snakemake.output[0]])
     run(['touch', '-h', snakemake.output[0]])
