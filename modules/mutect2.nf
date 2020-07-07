nextflow.preview.dsl=2

include BedToIntervalList from './gatkCNV.nf'
include PreprocessIntervals from './gatkCNV.nf'
include samtoolsIndex from './gatkCNV.nf'
include Funcotator from './single_varfilter.nf'

/*
params.min-base-quality-score =  20 
params.pcr-indel-model = 'AGGRESSIVE'
params.callable-depth = 14 
params.minimum-allele-fraction = 0.2 
params.base-quality-score-threshold = 20 
*/
process DownloadData {
container 'broadinstitute/gatk:4.1.7.0'
errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
maxRetries 100
label 'small'                
	output:
        path funcotator_dataSource

                """
                gatk FuncotatorDataSourceDownloader --somatic --validate-integrity --extract-after-download
                rm funcotator_dataSource*gz
                mv funcotator_dataSource* funcotator_dataSource
                """

}


process samtoolsRemoveSecondary {
        container "fredhutch/bwa:0.7.17"
//        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
  //      maxRetries 100
	label 'medium'

	input:
	tuple val(sampleID), val(kitID), val(type), val(patientID), file(bam_file)

	output:
	tuple val("${sampleID}"), val("${kitID}"), val("${type}"), val("${patientID}"), file("${sampleID}_primary.bam")


	"""
	samtools view -bh -f 0 -F 256 ${bam_file} > ${sampleID}_primary.bam
	"""
}
//       tuple val("${kitID}"), val("${patientID}"), val("${type}"), val("${sampleID}"), file("${bam}"), file("${bam}.bai")

process mutect2_normal_only {
        container 'broadinstitute/gatk:4.1.7.0'
        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 10
	label 'medium'
	
	input:
	tuple val(sampleID), val(kitID), val(type), val(patientID), file(bam_file), file(bam_index)
	file reference
        file reference_dict
        file reference_index

	output:
	tuple val("${sampleID}"), val("${kitID}"), val("${type}"), file("${sampleID}.vcf.gz")
	
	"""
	gatk Mutect2 --java-options "-Xmx30G" \
	-R ${reference} \
	-I ${bam_file} \
	-O ${sampleID}.vcf.gz

	"""

}




process GenomicsDBImport {
	container 'broadinstitute/gatk:4.1.7.0'
        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 100
	label 'medium'

	input:
	tuple val(kitID), file(interval), file(normal_list)
	file reference
	file reference_dict
	file reference_index

	output:
	tuple val("${kitID}"), path("${kitID}_db")	

"""
#!/usr/bin/python

import os
os.system('for F in *.vcf.gz ; do   tabix -f -p vcf \${F}  ; done')


files = '${normal_list}'
names = files.replace('.vcf.gz', '\t')
name_list = names.split(' ')
file_list = files.split(' ')
res = [i + j for i, j in zip(name_list, file_list)]
outF = open("map.txt", "w")

for line in res:
    outF.write(line)
    outF.write("\\n")

outF.close()

files = files.replace(' ', ' -V ')
     
script1 = 'gatk GenomicsDBImport  -R ${reference} -L ${interval} --genomicsdb-workspace-path ${kitID}_db'
       
script_f = script1 + ' --sample-name-map map.txt --merge-input-intervals'

os.system(script_f)
	   
"""
}


process CreateSomaticPanelOfNormals {
        container 'broadinstitute/gatk:4.1.7.0'
	errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
        maxRetries 100
	label 'medium'

	input:
	tuple val(kitID), path(db)
        file reference
        file reference_dict
        file reference_index

	output:
	tuple val("${kitID}"), file("${kitID}.vcf.gz"), file("${kitID}.vcf.gz.tbi")
	
	"""
	gatk CreateSomaticPanelOfNormals -R ${reference} -V gendb://${db} -O ${kitID}.vcf.gz
	"""
} 


process	mutect2_tumor_only {
	container 'broadinstitute/gatk:4.1.7.0'
	label 'medium'



	input:
	tuple val(kitID), path(normals), path(normal_index), val(sampleID), file(bam_file), file(bam_index)	
	file reference
	file reference_dict
	file reference_index
	file common_variants
	file common_variants_index


	output:
	tuple val("${sampleID}"), file("${sampleID}.vcf.gz"),  file("${sampleID}.vcf.gz.tbi"), file("${sampleID}.vcf.gz.stats")
	"""
  	gatk  --java-options "-Xmx30G" Mutect2 \
  	     -R ${reference} \
  	     -I ${bam_file} \
	     --panel-of-normals ${normals} \
	     --germline-resource ${common_variants} \
	     --min-base-quality-score ${params.min_base_quality_score} \
	     --pcr-indel-model ${params.pcr_indel_model} \
	     --callable-depth ${params.callable_depth} \
	     --minimum-allele-fraction ${params.minimum_allele_fraction} \
	     --base-quality-score-threshold ${params.base_quality_score_threshold} \
    	     -O ${sampleID}.vcf.gz

	"""
}

process mutect2_matched_normal {
        container 'broadinstitute/gatk:4.1.7.0'
//errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
//maxRetries 100
	label 'medium'

        input:
        tuple val(kitID), path(normals), path(normal_index), val(sampleID), file(bam_file), file(bam_index), val(normal_sampleID), file(normal_bam_file), file(normal_index)
        file reference
        file reference_dict
        file reference_index
        file common_variants
        file common_variants_index


        output:
        tuple val("${sampleID}"), file("${sampleID}.vcf.gz"),  file("${sampleID}.vcf.gz.tbi"), file("${sampleID}.vcf.gz.stats")
        """

	gatk AddOrReplaceReadGroups \
       	-I ${bam_file} \
       	-O ${sampleID}.bam \
       	--RGID none \
       	--RGLB none \
       	--RGPL none \
       	--RGPU none \
       	--RGSM ${sampleID}

        gatk AddOrReplaceReadGroups \
        -I ${normal_bam_file} \
        -O ${normal_sampleID}.bam \
        --RGID none \
        --RGLB none \
        --RGPL none \
        --RGPU none \
        --RGSM ${normal_sampleID}
	
	echo *
       samtools index ${sampleID}.bam
       samtools index ${normal_sampleID}.bam

        gatk  --java-options "-Xmx30G" Mutect2 \
             -R ${reference} \
             -I ${sampleID}.bam \
	     -I ${normal_sampleID}.bam \
	     -normal ${normal_sampleID} \
             --panel-of-normals ${normals} \
             --germline-resource ${common_variants} \
             --min-base-quality-score ${params.min_base_quality_score} \
             --pcr-indel-model ${params.pcr_indel_model} \
             --callable-depth ${params.callable_depth} \
             --minimum-allele-fraction ${params.minimum_allele_fraction} \
             --base-quality-score-threshold ${params.base_quality_score_threshold} \
             -O ${sampleID}.vcf.gz

        """
}


process FilterMutectCalls {
        container 'broadinstitute/gatk:4.1.7.0'
        errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
	maxRetries 100
	label 'medium'

        input:
        tuple val(sampleID), file(vcf),  file(vcf_index), file(stats)
        file reference
        file reference_dict
        file reference_index

        output:
	tuple val("${sampleID}"), file("${sampleID}_filtered.vcf.gz")

	"""
	gatk FilterMutectCalls \
	     -R ${reference} \
	     -V ${vcf} \
	     -O ${sampleID}_filtered.vcf.gz
	"""
}

process annotateVariants {
        container 'broadinstitute/gatk:4.1.7.0'
	label 'medium'

        input:
        tuple val(sampleID), file(vcf)
        path reference
        path reference_dict
        path reference_index
	path clinvar

        output:
        tuple val("${sampleID}"), file("${sampleID}_clinvar.vcf.gz")

        """
	 gatk IndexFeatureFile -I ${clinvar}
	 gatk IndexFeatureFile -I ${vcf}
	 gatk VariantAnnotator \
   	 		  -R ${reference} \
   			  -V ${vcf} \
   			  --output ${sampleID}_clinvar.vcf.gz \
			  --resource:ClinVar ${clinvar} \
			   -E ClinVar.AF_ESP \
			   -E ClinVar.AF_EXAC \
			   -E ClinVar.AF_TGP \
			   -E ClinVar.ALLELEID \
			   -E ClinVar.CLNDNINCL \
			   -E ClinVar.CLNDISDBINCL \
			   -E ClinVar.CLNHGVS \
			   -E ClinVar.CLNSIGINCL \
			   -E ClinVar.CLNVI \
			   -E ClinVar.CLNVC \
			   -E ClinVar.MC \
			   -E ClinVar.ORIGIN \
			   -E ClinVar.AF_ESP \
			   -E ClinVar.AF_EXAC \
			   -E ClinVar.AF_TGP \
			   -E ClinVar.ALLELEID \
			   -E ClinVar.CLNDNINCL \
			   -E ClinVar.CLNDISDBINCL \
			   -E ClinVar.CLNHGVS \
			   -E ClinVar.CLNSIGINCL \
			   -E ClinVar.CLNVI \
			   -E ClinVar.CLNVC \
			   -E ClinVar.MC \
			   -E ClinVar.ORIGIN \

        """
}

workflow mutect2_normals_wf {
	 take:
         bams_ch
         beds_ch
         reference
         reference_dict
         reference_index

	 main:

	 mutect2_normal_only(bams_ch.map{r ->[r[3], r[0], r[2], r[1], r[4], r[5]] }, reference, reference_dict, reference_index)
	 BedToIntervalList(beds_ch, reference, reference_dict)
	 normal_vcfs = mutect2_normal_only.out.map{ r -> [r[1], r[3]] }.groupTuple()
	 grouped_normals = BedToIntervalList.out.cross(normal_vcfs).map{ r -> [r[0][0], r[0][1], r[1][1]]}
	 GenomicsDBImport(grouped_normals, reference, reference_dict, reference_index)
	 CreateSomaticPanelOfNormals(GenomicsDBImport.out, reference, reference_dict, reference_index)
//       tuple val("${kitID}"), file("${kitID}.vcf.gz"), file("${kitID}.vcf.gz.tbi")	 
	 
	 emit:
	 PON = CreateSomaticPanelOfNormals.out

}


workflow mutect2_wf {
	 take:
	 bams_ch
	 beds_ch
	 reference
	 reference_dict
	 reference_index
	 common_variants
	 common_variants_index
	 clinvar


	 main:
	 samtoolsRemoveSecondary(bams_ch)
	 samtoolsIndex(samtoolsRemoveSecondary.out)


	 samtoolsIndex.out.branch {
                Normal : it[2] == 'Normal'
                Tumor : it[2] == 'Tumor'
               }.set { bams_branched }


	 mutect2_normals_wf(bams_branched.Normal, beds_ch, reference, reference_dict, reference_index)
	 PON = mutect2_normals_wf.out.PON
	 
	 tumor_bams = bams_branched.Tumor.map { r -> [r[1],r[0],r[3], r[4], r[5]] } // remove type
	 
	 normal_bams = bams_branched.Normal.map { r -> [r[1],r[0],r[3], r[4], r[5]] }
	 
	 unmatched_bams = normal_bams.cross(tumor_bams).map{ r -> [r[0][0], r[0][1]]}.join(tumor_bams, remainder : true).filter{it[1] == null }.map{r -> [r[2], r[3], r[4], r[5]]}

	 unmatched_pon = PON.cross(unmatched_bams).map{r -> [r[0][0], r[0][1], r[0][2], r[1][1], r[1][2], r[1][3]]}

	 mutect2_tumor_only(unmatched_pon, reference, reference_dict, reference_index, common_variants, common_variants_index)

 	 matched_bams = normal_bams.cross(tumor_bams).map{ r -> [r[0][1], r[1][2], r[1][3], r[1][4], r[0][2], r[0][3], r[0][4]] }

	 matched_pon = PON.cross(matched_bams).map{ r -> [r[0][0], r[0][1], r[0][2], r[1][1], r[1][2], r[1][3], r[1][4], r[1][5], r[1][6]] }

	 mutect2_matched_normal(matched_pon, reference, reference_dict, reference_index, common_variants, common_variants_index)
	 
	 mutect2_calls = mutect2_matched_normal.out.mix(mutect2_tumor_only.out)

	 FilterMutectCalls(mutect2_calls, reference, reference_dict, reference_index)
	 annotateVariants(FilterMutectCalls.out, reference, reference_dict, reference_index, clinvar)
	 DownloadData()
	 Funcotator(annotateVariants.out, reference, reference_index, reference_dict, DownloadData.out)

	 emit:
	 rawVCF = mutect2_calls
	 filteredVCF = FilterMutectCalls.out
	 annotated = annotateVariants.out
	 funcotated =  Funcotator.out 	 

}



