nextflow.preview.dsl=2

include gatkCNV_wf from './modules/gatkCNV.nf'
include mutect2_wf from './modules/mutect2.nf'




workflow {

	 reference = file(params.reference)
	 reference_index = file(params.reference_index)
	 reference_dict = file(params.reference_dict)

	 input_csv = file(params.input_csv)
	 output_folder = file(params.output_folder)

	 contig_dict = file(params.contig_dict)
	 common_variants = file(params.common_variants)
	 common_variants_index = file(params.common_variants_index)
	 clinvar = file(params.clinvar)


	bams_ch = Channel
            .fromPath(params.input_csv)
            .splitCsv(header:true)
            .map{ row-> tuple(row.sampleID, row.kitID, row.type, row.patientID, file(row.bam)) }


	beds_ch = Channel
       	     .fromPath(params.input_beds)
             .splitCsv(header:true)
             .map{ row-> tuple(row.kitID, file(row.capture_bed)) }


	main:
//        gatkCNV_wf(beds_ch, reference, reference_dict, reference_index, contig_dict, bams_ch)
	mutect2_wf(bams_ch, beds_ch, reference, reference_dict, reference_index, common_variants, common_variants_index, clinvar)

	publish:
	mutect2_wf.out.rawVCF to: "${params.output_folder}/rawVCF"
	mutect2_wf.out.filteredVCF to: "${params.output_folder}/filteredVCF"
	mutect2_wf.out.annotated to: "${params.output_folder}/annotatedVCF"
	mutect2_wf.out.funcotated to: "${params.output_folder}/funcotatedVCF"

}


