"""Get fastqc reports"""

def get_fastq(wildcards):
    return 'data/' + samples.loc[wildcards.sample, ["fq1", "fq2"]].dropna()

rule fastqc:
	"""Create fastqc report"""
	input: 
		R1 = 'data/{sample}_R1.fastq.gz',
		R2 = 'data/{sample}_R2.fastq.gz'
	output: 'logs/{sample}_R1_fastqc.html'
	conda: '../envs/fastqc.yaml'
	params:
		sge_opts="-cwd -j y -pe shared 2 -l h_data=6G,h_rt=2:00:00 -M "+email+" -m bea -N "+job_label+"{sample} -o cluster_logs -e cluster_logs",
		FASTQC = config['SOFTWARE']['FASTQC']
	shell:
		"""{params.FASTQC} {input.R1} {input.R2} -t 2 -o logs --extract
		sleep 500
		"""

rule multiqc_fastqc:
	input: expand('logs/{sample}_R1_fastqc.html', sample=samples.index)
	output:
		html = 'reports/fastqc.html',
		txt = 'reports/fastqc_data/multiqc_general_stats.txt'
	params:
		multiqc = config["SOFTWARE"]["MULTIQC"],
		sge_opts="-cwd -j y -l h_data=8G,h_rt=2:00:00 -M "+email+" -m bea -N "+job_label+" -o cluster_logs -e cluster_logs"
	shell:
		"""{params.multiqc} logs/ -m fastqc -o reports/ -n fastqc.html -f
		sleep 500
		"""