@click.command()
@click.option('--artifact_id', help='artifact id', required=True)
@click.option('--out_dir', help='output directory', required=True)
@click.option('--job_id', help='job id', required=True)
def runner(artifact_id, out_dir, job_id):
    qclient = client_connect('https://qiita.ucsd.edu')
    njobs = generate_sample_list(qclient, artifact_id, out_dir)
    generate_templates(out_dir, job_id, njobs)
    print(qclient.artifact_and_preparation_files(artifact_id))

    step0_job = run(['sbatch', f'{out_dir}/step-0/step-0.slurm'], stdout=PIPE)
    jid0 = step0_job.stdout.decode('utf8').split()[-1]

    step1_job = run(['sbatch', '-d', f'{out_dir}/step-1/step-1.slurm'], stdout=PIPE)
    jid1 = step1_job.stdout.decode('utf8').split()[-1]

    step2_job = run(['sbatch', '-d', f'aftercorr:{jid1}', f'{out_dir}/step-2/step-2.slurm'], stdout=PIPE)
    jid2 = step2_job.stdout.decode('utf8').split()[-1]

    step3_job = run(['sbatch', '-d', f'aftercorr:{jid2}', f'{out_dir}/step-3/step-3.slurm'], stdout=PIPE)
    jid3 = step3_job.stdout.decode('utf8').split()[-1]

    step4_job = run(['sbatch', '-d', f'aftercorr:{jid3}', f'{out_dir}/step-4/step-4.slurm'], stdout=PIPE)
    jid4 = step4_job.stdout.decode('utf8').split()[-1]

    step5_job = run(['sbatch', '-d', f'aftercorr:{jid4}', f'{out_dir}/step-5/step-5.slurm'], stdout=PIPE)
    jid5 = step5_job.stdout.decode('utf8').split()[-1]

    step6_job = run(['sbatch', '-d', f'aftercorr:{jid5}', f'{out_dir}/step-6/step-6.slurm'], stdout=PIPE)
    jid6 = step6_job.stdout.decode('utf8').split()[-1]

    step7_job = run(['sbatch', '-d', f'aftercorr:{jid6}', f'{out_dir}/step-7/step-7.slurm'], stdout=PIPE)
    jid7 = step7_job.stdout.decode('utf8').split()[-1]

    print(f'Submitted jobs: {jid0}, {jid1}, {jid2}, {jid3}, {jid4}, {jid5}, {jid6}')

if __name__ == '__main__':
        runner()
