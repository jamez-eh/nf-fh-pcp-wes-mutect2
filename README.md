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

