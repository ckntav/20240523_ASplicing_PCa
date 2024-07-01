workdir="/Users/chris/Desktop/20240409_NIPBL_project/output/chipseq-pipeline_NIPBL-GRCh38_SE"
antibody="NIPBL"
cutoff=0.05
cutoff_chr="0p05"
Mval=1
Mval_chr=1

time /Users/chris/miniconda3/bin/manorm \
    --verbose \
    --p1 $workdir/peak_call_withWCErep1/A549_DEX_${antibody}_rep1/${antibody}/A549_DEX_${antibody}_rep1.${antibody}_peaks.narrowPeak.stdchr.bed \
    --p2 $workdir/peak_call_withWCErep1/A549_CTRL_${antibody}_rep1/${antibody}/A549_CTRL_${antibody}_rep1.${antibody}_peaks.narrowPeak.stdchr.bed \
    --r1 $workdir/alignment/A549_DEX_${antibody}_rep1/${antibody}/A549_DEX_${antibody}_rep1.${antibody}.sorted.dup.filtered.bam \
    --r2 $workdir/alignment/A549_CTRL_${antibody}_rep1/${antibody}/A549_CTRL_${antibody}_rep1.${antibody}.sorted.dup.filtered.bam \
    --read-format bam \
    -m $Mval \
    -p $cutoff \
    -o $workdir/binding_diff/$antibody/A549_${antibody}_M${Mval_chr}_p${cutoff_chr} > $workdir/binding_diff/$antibody/A549_${antibody}_M${Mval_chr}_p${cutoff_chr}.log 2>&1


workdir="/media/fblab/Dataset/202405ZXXX"
PSIsigma_path="~/Programmes/PSI-Sogma-2.1/dummyai.pl"
gtf_path="Homo_sapiens.GRCh38.111.sorted.gtf"


type_read = 1
nread_min = 10
fmode_val = 3

perl PSIsigma_path 




VCap_RNA_EtOH_rep1.Aligned.sortedByCoord.out.bam
VCap_RNA_EtOH_rep2.Aligned.sortedByCoord.out.bam
VCap_RNA_EtOH_rep3.Aligned.sortedByCoord.out.bam

VCap_RNA_R1881_rep1.Aligned.sortedByCoord.out.bam
VCap_RNA_R1881_rep2.Aligned.sortedByCoord.out.bam
VCap_RNA_R1881_rep3.Aligned.sortedByCoord.out.bam
    
