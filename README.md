# nf-fh-pcp-wes-mutect2
Variant calling pipeline with mutect2


## input format
   
   ### input_csv:
   ```sampleID,kitID,type,patientID,bam```
   
   
   example:
   sample_1,v3+UTR,Tumor,001,sample_1.R1,sample_1.R2

   ### input_beds
   There is a parameter to supply bed files exome capture kits. These must be the same genome version as the alignment. It is important that the kitID in the 
   This csv file is supplied to the pipeline as a separate parameter ``--input_beds ./local_capture.csv ``
   It is merely a two column file, one column with kit ids (user supplied), this MUST match the kit id in the ```input_csv``` because this id is used to map each kit to the sample it was used for. During the pipeline a panel of normals is created based on ALL normal samples with an identical kit id and is used to subtract potential germline mutect2 calls. Therefore if there is a normal you do not want as part of the panel of normals, do not include it in the ```input_csv```.

   Example of this file type is as follows:
   ```
   kitID,capture_bed
   v2,/fh/scratch/delete90/nelson_p/james/references/capture_kits/V2.bed
   v3,/fh/scratch/delete90/nelson_p/james/references/capture_kits/V3.bed
   v3+UTR,/fh/scratch/delete90/nelson_p/james/references/capture_kits/V3_UTR.bed
   ```

## Details:
   
   This pipeline runs the basic Mutect2 detailed here: https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2

   The pipeline gathers ALL normals for each provided kitID (will produce a different panel for each capture kit provided). It also pairs matched normals to the corresponding tumors/PDX samples by joining on the provided ```patientID``` column. Please make sure that no samples are provided that were sequenced with different capture kits from the matched normal. This is not allowed due to biases in the exome capture kits. 

   The pipeline will subtract as many normal sources as provided. Ideally we will have a matched normal, PON, and a germline source such as gnomad.vcf.gz (required).


   Additionally this pipeline will Filter the mutect calls with filterMutectCalls https://gatk.broadinstitute.org/hc/en-us/articles/360036730391-FilterMutectCalls,

   Annotate variants with a provided variant file (this was obtained from https://www.ncbi.nlm.nih.gov/variation/docs/ClinVar_vcf_files/). Depending on the alignment genome format the chromosome names may have to be altered (I had to prepend 'chr' to the chromosome column)

   Call funcotator: https://gatk.broadinstitute.org/hc/en-us/articles/360040507011-Funcotator for further functional annotation. This command produces a maf file.

## Running

Can run straight from github with the following run file: 

```
bash run.sh
```

The run file should read as the following:
The 
```
#!/bin/bash                              
set -e

ml Singularity/3.5.3                     
export PATH=$SINGULARITYROOT/bin/:$PATH

BASE_BUCKET="/fh/scratch/delete90/nelson_p/james"

# Load the module                                                                                                                                 
ml nextflow

NXF_VER=20.01.0 nextflow \
    run \
    -resume \
    jamez-eh/nf-fh-pcp-wes-mutect2 \
    -profile local \
    -work-dir /fh/scratch/delete90/nelson_p/james/test_work \
    --input_csv ./test_manifest.txt \
    --output_folder ./out \
    --reference $BASE_BUCKET/references/hg38/Homo_sapiens_assembly38.fa \
    --reference_index $BASE_BUCKET/references/hg38/Homo_sapiens_assembly38.fa.fai \
    --reference_dict $BASE_BUCKET/references/hg38/Homo_sapiens_assembly38.dict \
    --common_variants $BASE_BUCKET/references/hg38/af-only-gnomad.hg38.vcf.gz \
    --common_variants_index $BASE_BUCKET/references/hg38/af-only-gnomad.hg38.vcf.gz.tbi \
    --ref_name hg38 \
    --input_beds local_capture.csv \
    --clinvar  $BASE_BUCKET/references/hg38/clinvar_20200520_chr.vcf \
```


