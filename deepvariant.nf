#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process FORMAT {
   
   publishDir "${params.output}/", mode: 'copy'

    input:
    path bam
    
    output:
    path "sample.txt", emit: sample 

    script:
    """
    ls bam | sed 's/\\..*//' |  sort --unique > sample.txt 

    """
   } 

process DEEPVARIANT {

	publishDir "${params.output}/vcf", mode: 'copy'

	container "google/deepvariant:1.4.0"
	cpus 8
	memory "16 GB"

	input:
	
	val(sample)
	path(bam)
	// path(bai)
	

	path(ref2)
	path(ref_fai2)
	path(ref_dict)
	path(ref_gzi)
	

	output:

	tuple val(sample), \
     path("${sample}.vcf.gz"), emit: vcf
	
	tuple val(sample), \
	 path("${sample}.vcf.gz.tbi"), emit: vcf_tbi

	tuple val(sample), \
	 path("${sample}.gvcf.gz"), emit: gvcf
	
	tuple val(sample), \
	 path("${sample}.gvcf.gz.tbi"), emit: gvcf_tbi

	path("${sample}.visual_report.html"), emit: report



	script:

	sample = sample[0]
		
		"""
		
    	/opt/deepvariant/bin/run_deepvariant \
    	--model_type=WES \
    	--ref="${ref2}" \
    	--reads="${bam}/${sample}.bam" \
    	--output_vcf="${sample}.vcf.gz" \
    	--output_gvcf="${sample}.gvcf.gz" \
    	--num_shards=8
 

		"""

}



process FILTER_VCF {

	publishDir "${params.output}/vcf_filter", mode: 'copy'
	container = 'ensemblorg/ensembl-vep'
	cpus = 4
	memory = 8.GB


	input:
		tuple val(sample), path(vcf)
		
	output:
		tuple val(sample), \
 		 path("${sample}.final.vcf"), emit: vcf

	script:

		"""
		filter_vep \
		-i ${vcf} --format vcf -o ${sample}.final.vcf \
		--filter "(FILTER = PASS) and \
		(CHROM in chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY)" \
		--force_overwrite

		"""
}

process CLINVAR {

	publishDir "${params.output}/clinvar", mode: 'copy'
	container "ensemblorg/ensembl-vep"
	errorStrategy 'ignore'

	cpus 8

	input:
		path(clinvar_old)
		path(clinvar_old_i)
		path(clinvar_new)
		path(clinvar_new_i)
		path(gNOMADg)
		path(gNOMADg_tbi)
		val vep_threads
		path vep_cache
		path vep_fasta
		path vep_fai
		path vep_gzi
		val vep_assembly
		tuple val(sample), path(vcf)
		
		

	output:
 		tuple val(sample), \
		 path("${sample}.tsv"), emit: prueba_vep

	script:
		

		"""

		vep \\
		--cache --offline --dir_cache ${vep_cache} \\
		--refseq --species homo_sapiens --assembly ${vep_assembly} \\
		--force_overwrite --use_transcript_ref --use_given_ref \\
		--verbose --fork ${vep_threads} --tab --format vcf --no_stats \\
		--fasta ${vep_fasta} \\
		--input_file ${vcf} \\
		--output_file ${sample}.vep.tsv \\
		--check_existing --canonical --numbers --hgvs --biotype --regulatory --symbol --protein \\
		--sift p --polyphen p --allele_number --variant_class --pubmed \\
		--custom ${clinvar_old},ClinVar_old,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \\
		--custom ${clinvar_new},ClinVar_new,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \\
		--custom ${gNOMADg},gnomADg,vcf,exact,0,AF,AC,AN,nhomalt,popmax,AF_popmax,AC_popmax,AF_nfe,AC_nfe,filt 

		
		grep -v "##" ${sample}.vep.tsv > ${sample}.tsv

		"""

}



workflow {

	//sample = Channel.fromPath(params.sample).splitCsv()

	
	FORMAT (
		params.bam
		)


	DEEPVARIANT (
		FORMAT.out.sample.splitCsv(),
		params.bam, 
		// params.bai,
		params.ref2,
		params.ref_fai2,
		params.ref_dict,
		params.ref_gzi
		)


	FILTER_VCF (		
		DEEPVARIANT.out.vcf 
			
	)

	CLINVAR (
		params.clinvar_old,
		params.clinvar_old_i,
		params.clinvar_new,
		params.clinvar_new_i,
		params.gNOMADg,
		params.gNOMADg_tbi,
		params.vep_threads,
		params.vep_cache,
		params.vep_fasta,
		params.vep_fai,
		params.vep_gzi,
		params.vep_assembly,
		FILTER_VCF.out.vcf
	
	)

}





