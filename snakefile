rulePath=config["general_path"]["RULE_PATH"]

include: rulePath+"/iget_samples_rule" #Temporaire
include: rulePath+"/trimmomatic_rule"
include: rulePath+"/fastqc_rule"
include: rulePath+"/multiqc_rule"
include: rulePath+"/bwa_mem_rule"
include: rulePath+"/samtools_sam_to_bam_rule"
include: rulePath+"/samtools_sort_rule"
include: rulePath+"/picard_markdup_rule"
include: rulePath+"/samtools_index_rule"
include: rulePath+"/pindel_one_file_rule"
include: rulePath+"/pindel2vcf_one_file_rule"
include: rulePath+"/gatk_realigner_target_creator_markdup_rule"
include: rulePath+"/gatk_indel_realigner_rule"
include: rulePath+"/gatk_base_recalibrator_rule"
include: rulePath+"/gatk_print_reads_rule"
include: rulePath+"/samtools_mpileup_one_file_rule"
include: rulePath+"/varscan_snp_rule"
include: rulePath+"/varscan_indel_rule"
include: rulePath+"/vcfcombine_rule"
include: rulePath+"/clintools_checkvariants_one_file_rule"
include: rulePath+"/vardict_rule"
include: rulePath+"/bcftools_index_one_file_rule"
include: rulePath+"/bcftools_call_one_file_rule"
include: rulePath+"/vt_decompose_one_file_rule"
include: rulePath+"/vt_normalize_one_file_rule"
include: rulePath+"/vt_uniq_final_one_file_rule"
include: rulePath+"/combineVCF2Leaves_v4_germline_rule"
include: rulePath+"/snpeff_rule"
include: rulePath+"/convert_rule"
include: rulePath+"/coverage_maker_sorted_rule"
include: rulePath+"/samtools_flagstat_rule"
include: rulePath+"/report_rule"


workdir: config["general_path"]["OUTPUT_PATH"]
input_path = config["general_path"]["INPUT_PATH"]
output_path = config["general_path"]["OUTPUT_PATH"]

import re
import glob

sample_ids = []
sampleName=glob.glob(input_path + "/*_R1_001.fastq.gz")
for name in sampleName:
    path = input_path + "/"
    name = name.replace(path, '')
    a = re.split('/', name)
    name = a[0]
    name = name.replace('_R1_001.fastq.gz', '')
    sample_ids.append(name)

mate_ids = ["R1","R2"]
pindel_types = ["D","BP","INV","TD","LI","SI","RP","CloseEndMapped","INT_final"]
callers = ["ST","VS","PL","VD"]

#iget_samples = expand((output_path+"/{sample_id}/{sample_id}_{mate_id}_001.fastq.gz"), sample_id =sample_ids, mate_id = mate_ids)
#trimmomatic = expand((output_path+"/{sample_id}/{sample_id}_R1_paired.fq.gz", output_path+"/{sample_id}/{sample_id}_R2_paired.fq.gz", output_path+"/{sample_id}/{sample_id}_R1_unpaired.fq.gz", output_path+"/{sample_id}/{sample_id}_R2_unpaired.fq.gz"), sample_id = sample_ids)
#fastqc = expand((output_path+"/multiqc/{sample_id}/{sample_id}_{mate_id}_paired_fastqc.html", output_path+"/multiqc/{sample_id}/{sample_id}_{mate_id}_paired_fastqc.zip"), sample_id = sample_ids, mate_id = mate_ids)
multiqc = (output_path+"/multiqc/results/multiqc_report.html")
#bwa_mem = expand((output_path+"/{sample_id}/{sample_id}.sam"), sample_id =sample_ids)
#sam2bam = expand((output_path+"/{sample_id}/{sample_id}.bam"), sample_id =sample_ids)
#samtools_sort = expand((output_path+"/{sample_id}/{sample_id}_sorted.bam"), sample_id =sample_ids)
#markduplicates = expand((output_path+"/{sample_id}/{sample_id}_marked.bam", output_path+"/{sample_id}/{sample_id}_metrics.txt"), sample_id = sample_ids)
#samtools_index = expand((output_path+"/{sample_id}/{sample_id}.bam.bai"), sample_id=sample_ids)
#realignertargetcreator = expand((output_path+"/{sample_id}/{sample_id}.intervals"), sample_id =sample_ids)
#indelrealigner = expand((output_path+"/{sample_id}/{sample_id}_realign.bam", output_path+"/{sample_id}/{sample_id}_realign.bai"), sample_id =sample_ids)
#gatk_base_recalibrator = expand(output_path+"/{sample_id}/{sample_id}{ext}_recalib.txt", sample_id=sample_ids, ext=["_fixed_sorted"])
#gatk_print_reads = expand((output_path+"/{sample_id}/{sample_id}_cleaned.bam", output_path+"/{sample_id}/{sample_id}_cleaned.bai"), sample_id=sample_ids)
#mpileup = expand(expand(output_path+"/{{sample_id}}/mpileup/{{sample_id}}_{caller}.{ext}", zip, caller=["ST","VS"], ext=["pileup","bcf"]), sample_id=sample_ids)
#varscan_snp = expand((output_path+"/{sample_id}/{sample_id}_snp.vcf"), sample_id =sample_ids)
#varscan_indel = expand((output_path+"/{sample_id}/{sample_id}_indel.vcf"), sample_id =sample_ids)
#vcfcombine = expand((output_path+"/{sample_id}/{sample_id}.vcf"), sample_id =sample_ids)
#clintools_checkvariants = expand(output_path+"/{sample_id}/{sample_id}_{caller}.tsv", sample_id=sample_ids, caller=["CT"])
#bcftools_index = expand(output_path+"/{sample_id}/mpileup/{sample_id}_ST.bcf.tbi", sample_id=sample_ids)
#bcftools_call = expand(output_path+"/{sample_id}/{sample_id}_ST.vcf", sample_id=sample_ids)
#pindel = expand(output_path+"/{sample_id}/pindel/{sample_id}_PL_{pindel_type}", sample_id=sample_ids, pindel_type=pindel_types)
#pindel2vcf = expand(output_path+"/{sample_id}/{sample_id}_PL.vcf", sample_id=sample_ids)
#vardict = expand(output_path+"/{sample_id}/{sample_id}_VD.vcf", sample_id=sample_ids)
#vt_decompose = expand(output_path+"/{sample_id}/{sample_id}_{caller}_VTD.vcf", sample_id=sample_ids, caller=callers)
#vt_normalize = expand(output_path+"/{sample_id}/{sample_id}_{caller}_VTN.vcf", sample_id=sample_ids, caller=callers)
#vt_uniq = expand(output_path+"/{sample_id}/{sample_id}_{caller}_final.vcf", sample_id=sample_ids, caller=callers)
#combineVCF2Leaves = expand(output_path+"/{sample_id}/{sample_id}_final.vcf", sample_id=sample_ids)
#snpeff = expand((output_path+"/{sample_id}/{sample_id}_snpeff.vcf", output_path+"/multiqc/{sample_id}/{sample_id}.csv"), sample_id = sample_ids)
convert = expand((output_path+"/{sample_id}/{sample_id}_converted.vcf"), sample_id =sample_ids)
coverage_maker = expand((output_path+"/{sample_id}/{sample_id}_coverage_maker.txt"), sample_id =sample_ids)
#flagstat = expand((output_path+"/multiqc/{sample_id}/{sample_id}.txt"), sample_id =sample_ids)
report = expand((output_path+"/{sample_id}/{sample_id}.html"), sample_id =sample_ids)

rule all:
    input:
#        iget_samples,
#        trimmomatic,
#        fastqc,
        multiqc,
#        bwa_mem,
#        sam2bam,
#        samtools_sort,
#        markduplicates,
#        index,
#        realignertargetcreator,
#        indelrealigner,
#        gatk_base_recalibrator,
#        gatk_print_reads,
#        mpileup,
#        varscan_snp,
#        varscan_indel,
#        vcfcombine,
#        clintools_checkvariants,
#        bcftools_index,
#        bcftools_call,
#        pindel,
#        pindel2vcf,
#        vardict,
#        vt_decompose,
#        vt_normalize,
#        vt_uniq,
#        combineVCF2Leaves,
#        snpeff,
        convert,
        coverage_maker,
#        flagstat,
        report
    shell:
        "touch " + output_path + "/done"
