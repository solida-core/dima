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
    if snakemake.params['output_fmt'] == 'CRAM':
        _opt_ = '-C'
    else:
        _opt_ = '-b'
    cmd = [snakemake.params['cmd'], 'view', _opt_,
           '--threads', str(snakemake.threads),
           '-O', snakemake.params['output_fmt']]
    if snakemake.params['output_fmt'] == 'CRAM':
        cmd.append('--reference')
        cmd.append(snakemake.params['genome'][0])
    cmd.append('-o')
    cmd.append(snakemake.output[0])
    cmd.append(snakemake.input[0])
    run(cmd)
    run(['touch', '-h', snakemake.output[0]])
