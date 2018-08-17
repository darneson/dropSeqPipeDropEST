"""Filter data"""

rule fastq_to_sam:
	"""Create an empty bam file linking cell/UMI barcodes to reads"""
	input:
		R1 = 'data/{sample}_R1.fastq.gz',
		R2 = 'data/{sample}_R2.fastq.gz'
	output:
		temp('data/{sample}_unaligned.bam')
	params:
		TMPDIR = config['LOCAL']['TMPDIR'],
		memory = config['LOCAL']['MEMORY'],
		picard = config['SOFTWARE']['PICARD'],
		java = config['SOFTWARE']['JAVA'],
		sge_opts = "-cwd -j y -l h_data=12G,h_rt=4:00:00 -M "+email+" -m bea -N "+job_label+"{sample} -o cluster_logs -e cluster_logs"
	shell:
		"""{params.java} -Djava.io.tmpdir={params.TMPDIR} -Xmx{params.memory} -jar {params.picard} FastqToSam\
		F1={input[0]}\
		F2={input[1]}\
		SM=DS O={output}
		sleep 500
		"""

rule quality_tags:
	"""Add tags to .bam file for cell and UMI barcodes quality scores which can be used for dropEST"""
	input:
		'data/{sample}_unaligned.bam'
	output:
		temp('data/{sample}_unaligned_qualityTagged.bam')
	params:
		CB_start = config['FILTER']['Cell_barcode']['start'],
		CB_end = config['FILTER']['Cell_barcode']['end'],
		UMI_start = config['FILTER']['UMI']['start'],
		UMI_end = config['FILTER']['UMI']['end'],
		Quality_Tag=config['SOFTWARE']['QUALITY_TAG'],
		Python=config['SOFTWARE']['PYTHON_PATH'],
		sge_opts = "-cwd -j y -l h_data=8G,h_rt=6:00:00 -M "+email+" -m bea -N "+job_label+"{sample} -o cluster_logs -e cluster_logs"
	shell:
		"""{params.Python} {params.Quality_Tag}\
		-CBBR {params.CB_start}-{params.CB_end}\
		-UMIBR {params.UMI_start}-{params.UMI_end}\
		-I {input}\
		-O {output}\
		-CBTN CQ\
		-UMITN MQ\
		-BCD 1
		sleep 500
		"""

rule BC_tags:
	input:
		'data/{sample}_unaligned_qualityTagged.bam'
	output: 
		bam = temp('data/{sample}_BC_tagged_unmapped.bam'),
		BC_summary = 'logs/{sample}_CELL_barcode.txt'
	params:
		BC_start = config['FILTER']['Cell_barcode']['start'],
		BC_end = config['FILTER']['Cell_barcode']['end'],
		BC_min_quality = config['FILTER']['Cell_barcode']['min_quality'],
		BC_min_quality_num = config['FILTER']['Cell_barcode']['num_below_quality']+1,
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory = config['LOCAL']['MEMORY'],
		TMPDIR = config['LOCAL']['TMPDIR'],
		sge_opts = "-cwd -j y -l h_data=12G,h_rt=2:00:00 -M "+email+" -m bea -N "+job_label+"{sample} -o cluster_logs -e cluster_logs"
	shell:
		"""{params.DROPSEQ_wrapper} -m {params.memory} -t {params.TMPDIR} -p TagBamWithReadSequenceExtended\
		SUMMARY={output.BC_summary}\
		BASE_RANGE={params.BC_start}-{params.BC_end}\
		BASE_QUALITY={params.BC_min_quality}\
		BARCODED_READ=1\
		DISCARD_READ=false\
		TAG_NAME=XC\
		NUM_BASES_BELOW_QUALITY={params.BC_min_quality_num}\
		INPUT={input}\
		OUTPUT={output.bam}
		sleep 500
		"""

rule UMI_tags:
	input: 'data/{sample}_BC_tagged_unmapped.bam'
	output:
		bam = temp('data/{sample}_BC_UMI_tagged_unmapped.bam'),
		UMI_summary = 'logs/{sample}_UMI_barcode.txt'
	params:
		UMI_start = config['FILTER']['UMI']['start'],
		UMI_end = config['FILTER']['UMI']['end'],
		UMI_min_quality = config['FILTER']['UMI']['min_quality'],
		UMI_min_quality_num = config['FILTER']['UMI']['num_below_quality']+1,
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory = config['LOCAL']['MEMORY'],
		TMPDIR = config['LOCAL']['TMPDIR'],
		sge_opts = "-cwd -j y -l h_data=8G,h_rt=2:00:00 -M "+email+" -m bea -N "+job_label+"{sample} -o cluster_logs -e cluster_logs"
	shell:
		"""{params.DROPSEQ_wrapper} -m {params.memory} -t {params.TMPDIR} -p TagBamWithReadSequenceExtended\
		SUMMARY={output.UMI_summary}\
		BASE_RANGE={params.UMI_start}-{params.UMI_end}\
		BASE_QUALITY={params.UMI_min_quality}\
		BARCODED_READ=1\
		DISCARD_READ=true\
		TAG_NAME=XM\
		NUM_BASES_BELOW_QUALITY={params.UMI_min_quality_num}\
		INPUT={input}\
		OUTPUT={output.bam}
		sleep 500
		"""

rule filter_tags:
	input: 'data/{sample}_BC_UMI_tagged_unmapped.bam'
	output: 
		temp('data/{sample}_tags_filtered_unmapped.bam')
	params:
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory = config['LOCAL']['MEMORY'],
		TMPDIR = config['LOCAL']['TMPDIR'],
		sge_opts = "-cwd -j y -l h_data=8G,h_rt=2:00:00 -M "+email+" -m bea -N "+job_label+"{sample} -o cluster_logs -e cluster_logs"
	shell:
		"""{params.DROPSEQ_wrapper} -m {params.memory} -t {params.TMPDIR} -p FilterBAM\
		TAG_REJECT=XQ\
		INPUT={input}\
		OUTPUT={output}
		sleep 500
		"""

rule count_reads:
	input: 'data/{sample}_tags_filtered_unmapped.bam'
	output: 'logs/{sample}_reads_left.txt'
	params:
		SAMTOOLS=config['SOFTWARE']['SAMTOOLS'],
		sge_opts="-cwd -j y -l h_data=4G,h_rt=2:00:00 -M "+email+" -m bea -N "+job_label+"{sample} -o cluster_logs -e cluster_logs"
	shell:
		"""{params.SAMTOOLS} view {input} | wc -l > {output}
		sleep 500
		"""

rule start_trim:
	input: 'data/{sample}_tags_filtered_unmapped.bam'
	output:
		bam = temp('data/{sample}_tags_start_filtered_unmapped.bam'),
		trim_summary = 'logs/{sample}_start_trim.txt'
	params:
		SmartAdapter = config['FILTER']['5PrimeSmartAdapter'],
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory = config['LOCAL']['MEMORY'],
		TMPDIR = config['LOCAL']['TMPDIR'],
		sge_opts="-cwd -j y -l h_data=8G,h_rt=2:00:00 -M "+email+" -m bea -N "+job_label+"{sample} -o cluster_logs -e cluster_logs"
	shell:
		"""{params.DROPSEQ_wrapper} -m {params.memory} -t {params.TMPDIR} -p TrimStartingSequence\
		OUTPUT_SUMMARY={output.trim_summary}\
		SEQUENCE={params.SmartAdapter}\
		MISMATCHES=1\
		NUM_BASES=6\
		INPUT={input}\
		OUTPUT={output.bam}
		sleep 500
		"""

rule polya_trim:
	input: 'data/{sample}_tags_start_filtered_unmapped.bam'
	output:
		bam =temp('data/{sample}_trimmed_unmapped.bam'),
		trim_summary = 'logs/{sample}_polyA_trim.txt'
	params:
		DROPSEQ_wrapper=config['LOCAL']['DROPSEQ-wrapper'],
		memory = config['LOCAL']['MEMORY'],
		TMPDIR = config['LOCAL']['TMPDIR'],
		sge_opts="-cwd -j y -l h_data=8G,h_rt=2:00:00 -M "+email+" -m bea -N "+job_label+"{sample} -o cluster_logs -e cluster_logs"
	shell:
		"""{params.DROPSEQ_wrapper} -m {params.memory} -t {params.TMPDIR} -p PolyATrimmer\
		OUTPUT_SUMMARY={output.trim_summary}\
		MISMATCHES=0\
		NUM_BASES=5\
		INPUT={input}\
		OUTPUT={output.bam}
		sleep 500
		"""

rule sam_to_fastq:
	input: 'data/{sample}_trimmed_unmapped.bam'
	output: temp('data/{sample}_trimmed_unmapped.fastq.gz')
	params:
		TMPDIR = config['LOCAL']['TMPDIR'],
		memory = config['LOCAL']['MEMORY'],
		picard = config['SOFTWARE']['PICARD'],
		java = config['SOFTWARE']['JAVA'],
		sge_opts="-cwd -j y -l h_data=8G,h_rt=2:00:00 -M "+email+" -m bea -N "+job_label+"{sample} -o cluster_logs -e cluster_logs"
	shell:
		"""{params.java} -Xmx{params.memory} -jar -Djava.io.tmpdir={params.TMPDIR}	{params.picard} SamToFastq\
		INPUT={input}\
		FASTQ=/dev/stdout COMPRESSION_LEVEL=0|\
		gzip > {output}
		sleep 500
		"""

rule trim_single:
	input: 'data/{sample}_trimmed_unmapped.fastq.gz'
	output:
		data = temp('data/{sample}_filtered.fastq.gz'),
		log = 'logs/{sample}_trimlog.txt'
	params:
		ILLUMINACLIP = config['FILTER']['SeqAdapters'],
		TRIMMOMATIC = config['SOFTWARE']['TRIMMOMATIC'],
		java = config['SOFTWARE']['JAVA'],
		sge_opts="-cwd -j y -pe shared 2 -l h_data=6G,h_rt=2:00:00 -M "+email+" -m bea -N "+job_label+"{sample}"
		# sge_opts="-cwd -j y -pe shared 2 -l h_data=6G,h_rt=2:00:00 -M "+email+" -m bea -N "+job_label+"{sample} -o cluster_logs -e cluster_logs"
	shell:
		"""{params.java} -jar {params.TRIMMOMATIC} \
		SE {input} {output.data} \
		-threads 2 \
		ILLUMINACLIP:{params.ILLUMINACLIP}:2:30:10 \
		LEADING:3 \
		TRAILING:3 \
		SLIDINGWINDOW:4:20 \
		MINLEN:15 > {output.log} 2>&1
		sleep 500
		"""

rule plot_polyA_trim:
	input: 'logs/{sample}_polyA_trim.txt'
	output:
		pdf = 'plots/{sample}_polya_trimmed.pdf',
		png = 'plots/png/{sample}_polya_trimmed.png'
	params:
		sge_opts="-cwd -j y -l h_data=8G,h_rt=2:00:00 -M "+email+" -m bea -N "+job_label+"{sample} -o cluster_logs -e cluster_logs"
	script:
		'../scripts/plot_polyA_trim.R'

rule plot_barcode_start_trim:
	input: 'logs/{sample}_start_trim.txt'
	output:
		pdf = 'plots/{sample}_start_trim.pdf',
		png = 'plots/png/{sample}_start_trim.png'
	params:
		sge_opts="-cwd -j y -l h_data=8G,h_rt=2:00:00 -M "+email+" -m bea -N "+job_label+"{sample} -o cluster_logs -e cluster_logs"
	script:
		'../scripts/plot_start_trim.R'


rule plot_UMI_filtering:
	input: 'logs/{sample}_UMI_barcode.txt'
	output:
		pdf = 'plots/{sample}_UMI_dropped.pdf',
		png = 'plots/png/{sample}_UMI_dropped.png'
	params: 
		min_quality = config['FILTER']['UMI']['min_quality'],
		num_below_quality = config['FILTER']['UMI']['num_below_quality'],
		sge_opts="-cwd -j y -l h_data=8G,h_rt=2:00:00 -M "+email+" -m bea -N "+job_label+"{sample} -o cluster_logs -e cluster_logs"
	script:
		'../scripts/plot_umi_drop.R'

rule plot_CELL_filtering:
	input: 'logs/{sample}_CELL_barcode.txt'
	output:
		pdf = 'plots/{sample}_CELL_dropped.pdf',
		png = 'plots/png/{sample}_CELL_dropped.png'
	params:
		min_quality = config['FILTER']['Cell_barcode']['min_quality'],
		num_below_quality = config['FILTER']['Cell_barcode']['num_below_quality'],
		sge_opts="-cwd -j y -l h_data=8G,h_rt=2:00:00 -M "+email+" -m bea -N "+job_label+"{sample} -o cluster_logs -e cluster_logs"
	script:
		'../scripts/plot_cell_drop.R'

rule plot_BC_drop:
	input:
		BC_tagged = expand('logs/{sample}_CELL_barcode.txt', sample=samples.index),
		UMI_tagged = expand('logs/{sample}_UMI_barcode.txt', sample=samples.index),
		reads_left = expand('logs/{sample}_reads_left.txt', sample=samples.index)
	output:
		pdf = 'plots/BC_drop.pdf',
		png = 'plots/png/BC_drop.png'
	params:
		BC_length = config['FILTER']['Cell_barcode']['end'] - config['FILTER']['Cell_barcode']['start']+1,
		UMI_length = config['FILTER']['UMI']['end'] - config['FILTER']['UMI']['start']+1,
		min_num_below_BC = config['FILTER']['Cell_barcode']['num_below_quality'],
		min_num_below_UMI = config['FILTER']['UMI']['num_below_quality'],
		min_BC_quality = config['FILTER']['Cell_barcode']['min_quality'],
		min_UMI_quality = config['FILTER']['UMI']['min_quality'],
		sample_names = lambda wildcards: samples.index,
		sge_opts="-cwd -j y -l h_data=8G,h_rt=2:00:00 -M "+email+" -m bea -N "+job_label+" -o cluster_logs -e cluster_logs"
	script:
		'../scripts/plot_BC_drop.R'

rule multiqc_trimmomatic:
	input: expand('logs/{sample}_trimlog.txt', sample=samples.index)
	output: 'reports/filter.html'
	params:
		multiqc = config["SOFTWARE"]["MULTIQC"],
		sge_opts="-cwd -j y -l h_data=8G,h_rt=2:00:00 -M "+email+" -m bea -N "+job_label+" -o cluster_logs -e cluster_logs"
	shell:
		"""{params.multiqc} logs/ -m trimmomatic -o reports/ -n filter.html -f
		sleep 500
		"""