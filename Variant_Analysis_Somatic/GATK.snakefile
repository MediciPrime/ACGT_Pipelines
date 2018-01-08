configfile:  "config.yaml"

targets = []
for sample in config['samples']:
    targets.extend(
        expand(
            'data/bam_files/{reference}/{sample}/{sample}_recalb.bam',
            sample=sample,
            reference=config['reference']
        )
    )
    targets.extend(
        expand(
            'data/variants/{reference}/normal_v_tumor.table',
            reference=config['reference']
        )
    )

rule all:
    input:
        targets

rule cutadapt:
    input:
        r1 = 'data/raw_reads/{sample}_R1.fastq.gz',
        r2 = 'data/raw_reads/{sample}_R2.fastq.gz'
    output:
        r1 = 'data/trimmed_reads/{sample}_R1_trimmed.fq',
        r2 = 'data/trimmed_reads/{sample}_R2_trimmed.fq'
    params:
        five_prime = 'GTTCAGAGTTCTACAGTCCGACGATC',
        three_prime = 'TGGAATTCTCGGGTGCCAAGG'
    shell:
        'cutadapt -a {params.five_prime} -A {params.three_prime} '
        '-q 20 -m 25 -o {output.r1} -p {output.r2} '
        '{input.r1} {input.r2}'

rule bwa_index:
    input:
        ref_fasta = 'data/reference/{reference}.fasta'
    output:
        'data/reference/{reference}.fasta.bwt'
    shell:
        'bwa index {input.ref_fasta}'

rule bwa_align:
    input:
        ref_index = 'data/reference/{reference}.fasta.bwt',
        ref_fasta = 'data/reference/{reference}.fasta',
        r1 = 'data/trimmed_reads/{sample}_R1_trimmed.fq',
        r2 = 'data/trimmed_reads/{sample}_R2_trimmed.fq'
    output:
        'data/sam_files/{reference}/{sample}/{sample}.sam'
    params:
        reference = config['reference']
    threads: 10
    shell:
        'bwa mem {input.ref_fasta} {input.r1} {input.r2} > {output}'

def get_bar(wc):
    barcode = config['samples'][wc.sample]['barcode']
    return barcode
        
rule bwa_picard_rg:
    input:
        sam_file = 'data/sam_files/{reference}/{sample}/{sample}.sam'
    output:
        'data/bam_files/{reference}/{sample}/{sample}_rg_sorted.bam'
    params:
        barcode=get_bar
    shell:
        'picard AddOrReplaceReadGroups I={input.sam_file} O={output} '
        'SORT_ORDER=coordinate RGLB={sample} RGPL=illumina RGPU={params.barcode} '
        'RGSM={sample}'
        
rule bwa_picard_md:
    input:
        bam_sorted = 'data/bam_files/{reference}/{sample}/{sample}_rg_sorted.bam'
    output:
        bam_dedup = 'data/bam_files/{reference}/{sample}/{sample}_rg_dedupped.bam',
        output_met = 'data/bam_files/{reference}/{sample}/{sample}_output.metrics'
    shell:
        'picard MarkDuplicates I={input.bam_sorted} O={output.bam_dedup} '
        'CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M={output.output_met}'

rule fasta_index:
    input:
        ref = 'data/reference/{reference}.fasta'
    output:
        ref_index = 'data/reference/{reference}.fasta.fai'
    shell:
        'samtools faidx {input.ref}'

rule fasta_dictionary:
    input:
        ref = 'data/reference/{reference}.fasta'
    output:
        ref_dict = 'data/reference/{reference}.dict'
    shell:
        'picard CreateSequenceDictionary R={input.ref}'      
        
rule gatk_baseREcalab:
    input:
        bam_dedup = 'data/bam_files/{reference}/{sample}/{sample}_rg_dedupped.bam',
        ref_fasta = 'data/reference/{reference}.fasta',
        ref_index = 'data/reference/{reference}.fasta.fai',
        ref_dict = 'data/reference/{reference}.dict',
        gold_snps = 'data/reference/gold_snps.vcf',
        gold_indels = 'data/reference/gold_indels.vcf'
    output:
        covar_table = 'data/base_recalab/{reference}/{sample}/recal_data.table'
    threads: 30
    shell:
        'java -jar GenomeAnalysisTK.jar -T BaseRecalibrator '
        '-nct {threads} -R {input.ref_fasta} -I {input.bam_dedup} '
        '-knownSites {input.gold_snps} -knownSites {input.gold_indels} '
        '-o {output.covar_table} '

rule gatk_secondREcalab:
    input:
        bam_dedup = 'data/bam_files/{reference}/{sample}/{sample}_rg_dedupped.bam',
        ref = 'data/reference/{reference}.fasta',
        gold_snps = 'data/reference/gold_snps.vcf',
        gold_indels = 'data/reference/gold_indels.vcf',
        covar_table = 'data/base_recalab/{reference}/{sample}/recal_data.table',
    output:
        covar_recal = 'data/base_recalab/{reference}/{sample}/post_recal_data.table'
    threads: 30
    shell:
        'java -jar GenomeAnalysisTK.jar -T BaseRecalibrator '
        '-nct {threads} -R {input.ref} -I {input.bam_dedup} '
        '-knownSites {input.gold_snps} -knownSites {input.gold_indels} '
        '-BQSR {input.covar_table} -o {output.covar_recal}'

rule gatk_plots:
    input:
        ref = 'data/reference/{reference}.fasta',
        covar_table = 'data/base_recalab/{reference}/{sample}/recal_data.table',
        covar_recal = 'data/base_recalab/{reference}/{sample}/post_recal_data.table',
    output:
        plots = 'data/base_recalab/{reference}/{sample}/recalibration_plots.pdf'
    threads: 30
    shell:
        'java -jar GenomeAnalysisTK.jar -T AnalyzeCovariates '
        '-nct {threads} -R {input.ref} -before {input.covar_table} '
        '-after {input.covar_recal} -plots {output.plots}'

rule gatk_applyREcalab:
    input:
        ref = 'data/reference/{reference}.fasta',
        bam_dedup = 'data/bam_files/{reference}/{sample}/{sample}_rg_dedupped.bam',
        covar_table = 'data/base_recalab/{reference}/{sample}/recal_data.table'
    output:
        bam_recalb = 'data/bam_files/{reference}/{sample}/{sample}_recalb.bam'
    threads: 30
    shell:
        'java -jar GenomeAnalysisTK.jar -T PrintReads -nct {threads} -R {input.ref} '
        '-I {input.bam_dedup} -BQSR {input.covar_table} '
        '-o {output.bam_recalb}'
        
rule gatk_mutect2:
    input:
        ref = 'data/reference/{reference}.fasta',
        tumor_bam = 'data/bam_files/{reference}/tumor_1/tumor_1_recalb.bam',
        normal_bam = 'data/bam_files/{reference}/normal_1/normal_1_recalb.bam'
    output:
        variants = 'data/variants/{reference}/normal_v_tumor.vcf'
    threads: 60
    shell:
        'java -jar GenomeAnalysisTK.jar -T MuTect2 -nct {threads} -R {input.ref} '
        '-I:tumor {input.tumor_bam} -I:normal {input.normal_bam} '
        '--max_alt_alleles_in_normal_count 10000000 '
        '--max_alt_allele_in_normal_fraction 0.10 -o {output.variants}'

rule filtered_vcf:
    input:
        variants = 'data/variants/{reference}/normal_v_tumor.vcf'
    output:
        filtered_vcf = 'data/variants/{reference}/normal_v_tumor_filtered.vcf'
    shell:
        "cat {input.variants} | SnpSift filter ' ( FILTER = 'PASS')' > {output}"

rule variants_table:
    input:
        ref = 'data/reference/{reference}.fasta',
        filtered_vcf = 'data/variants/{reference}/normal_v_tumor_filtered.vcf'
    output:
        variant_table = 'data/variants/{reference}/normal_v_tumor.table'
    shell:
        'java -jar GenomeAnalysisTK.jar -T VariantsToTable -R {input.ref} '
        '-V {input.filtered_vcf} -F CHROM -F POS -F REF -F ALT -GF AD -GF AF '
        '-EMD -o {output.variant_table} '
