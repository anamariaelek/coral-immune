#!/bin/bash
#$ -N map_sample
#$ -cwd
#$ -q long-sl7,mem_512_12h,mem_512
#$ -l virtual_free=20G,h_rt=12:00:00,disk=10G
#$ -o map_sample.out
#$ -j y
#$ -pe smp 4


if [ -z "$6" ] ; then echo -e "
Wait! You need to specify:
	1) run ID (user-defined string)
	2) path to gene annotation (BED or GTF)
	3) path to genome assembly
	4) path to reads (1st pairs)
	5) path to reads (2nd pairs)
	6) number of threads
	7) OPTIONAL: output folder (defaults to results_atac_{run ID})
	"
	exit
fi

pid="$1"
gtfo="$2"
faso="$3"
re1o="$4"
re2o="$5"
nth="$6"
timestamp=$(date +%s)
if [ -z $7 ]; then
	dir="results_rna_${pid}/"
else
	dir="$7"
fi
log=${dir}/log_${pid}.${timestamp}.log


# create folders to store data
mkdir -p ${dir}/reference
mkdir -p ${dir}/reads
mkdir -p ${dir}/alignments
mkdir -p ${dir}/qualitycontrol

echo -e "# Input parameters
# run ID       = ${pid}
# GTF          = ${gtfo}
# Genome fasta = ${faso}
# Reads 1      = ${re1o}
# Reads 2      = ${re2o}
# Num threads  = ${nth}
# Timestamp    = ${timestamp}
# Output dir   = ${dir}
# Working dir  = $(pwd)
" |& tee ${log}

# function
check_link_file() {
	local f1=$1
	local f2=$2

    if ! cmp --silent $f1 $f2
	# if there is no f2 or if f2 is different from source f1, relink
	# if f2 is identical to source f1, do nothing
		then
			# echo "# symbolic link from $f1 to $f2"
			rm -f $f2
            ln -s $f1 $f2
        else
        	echo "# $f2 already in place, no need to relink"
	fi
}


# PREPARE INPUTS
# reads
echo "# Reads $re1o --> ${dir}/reads/${pid}_1.fastq.gz" |& tee -a ${log}
check_link_file ${re1o} ${dir}/reads/${pid}_1.fastq.gz  |& tee -a ${log}
echo "# Reads $re2o --> ${dir}/reads/${pid}_2.fastq.gz" |& tee -a ${log}
check_link_file ${re2o} ${dir}/reads/${pid}_2.fastq.gz  |& tee -a ${log}

# reference genome
fas=${dir}/reference/${pid}_gDNA.fasta
echo "# Reference genome ${faso} --> ${fas}" |& tee -a ${log} 
check_link_file ${faso} ${fas}               |& tee -a ${log}

# fasta index 
if [[ ! -f "${faso}.fai" ]]
then
  echo "# No fasta index found, indexing ${fas}..." |& tee -a ${log}
  samtools faidx ${fas}  |& tee -a ${log}
else
  check_link_file ${faso}.fai ${fas}.fai
fi

# annotation
gtf=${dir}/reference/${pid}_annot.gtf
echo "# Reference annotation ${gtfo} --> ${gtf}" |& tee -a ${log}
check_link_file ${gtfo} ${gtf}                   |& tee -a ${log}

# genome index
fadr=$( dirname ${faso} )
fado=${fadr}/STAR_index
if [[ ! -d "${fado}" ]]
then
  echo "# No STAR_index found, indexing ${fas}..." |& tee -a ${log}
  mkdir -p ${dir}/reference/STAR_index
  STAR --runThreadN ${nth} \
	--runMode genomeGenerate \
	--genomeDir ${dir}/reference/STAR_index \
	--genomeFastaFiles ${fas} \
	--sjdbGTFfile ${gtf} \
	--sjdbOverhang 99 \
	|& tee -a ${log}
else
  fad=${dir}/reference/STAR_index
  echo "# Genome index ${fado} --> ${fad}" |& tee -a ${log}     
  check_link_file ${fado} ${fad}           |& tee -a ${log}
fi

# READ QC
# quality control
echo "# FASTQC ${pid}, original reads..."                                                   |& tee -a ${log}
fastqc -t ${nth} ${dir}/reads/${pid}_1.fastq.gz ${dir}/reads/${pid}_2.fastq.gz 2> /dev/null |& tee -a ${log}

# create mock FASTA with common adapters
echo -e ">PrefixNX/1
AGATGTGTATAAGAGACAG
>PrefixNX/2
AGATGTGTATAAGAGACAG
>Trans1
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
>Trans1_rc
CTGTCTCTTATACACATCTGACGCTGCCGACGA
>Trans2
GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
>Trans2_rc
CTGTCTCTTATACACATCTCCGAGCCCACGAGAC" > ${dir}/reads/${pid}_adapters_NexteraPE-PE.fasta

# trim reads
echo "# TRIMMOMATIC trimming ${pid} reads..." |& tee -a ${log}
trimmomatic PE \
        -threads ${nth} \
        -phred33 \
        ${dir}/reads/${pid}_1.fastq.gz \
        ${dir}/reads/${pid}_2.fastq.gz \
        ${dir}/reads/${pid}_1.tp.fastq.gz \
        ${dir}/reads/${pid}_1.tu.fastq.gz \
        ${dir}/reads/${pid}_2.tp.fastq.gz \
        ${dir}/reads/${pid}_2.tu.fastq.gz \
        SLIDINGWINDOW:4:20 \
        ILLUMINACLIP:${dir}/reads/${pid}_adapters_NexteraPE-PE.fasta:2:30:10 |& tee -a ${log}

# recheck with fastqc
echo "# FASTQC ${pid}, trimmed reads..."                                             |& tee -a ${log}
fastqc -t ${nth} ${dir}/reads/${pid}_1.tp.fastq.gz ${dir}/reads/${pid}_2.tp.fastq.gz |& tee -a ${log}


# ALIGNMENT
echo "# STAR mapping ${pid} reads..." |& tee -a ${log}
#/home/anamaria/bin/STAR-2.7.10a/bin/Linux_x86_64/STAR \
STAR \
  	--genomeDir ${dir}/reference/STAR_index \
    --runThreadN ${nth} \
    --readFilesIn ${dir}/reads/${pid}_1.tp.fastq.gz ${dir}/reads/${pid}_2.tp.fastq.gz \
  	--readFilesCommand zcat \
    --outFileNamePrefix ${dir}/alignments/${pid}. \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMattributes Standard \
  	--outFilterMultimapNmax 20 \
	--limitBAMsortRAM 50000000000 \
  	--outFilterScoreMinOverLread 0.33 \
  	--outFilterMatchNminOverLread 0.33 \
  	--outFilterMatchNmin 0 \
  	--outFilterMismatchNmax 5 \
  	--alignIntronMax 5000 \
  	--genomeLoad LoadAndRemove \
  	--readNameSeparator ' ' \
	--quantMode GeneCounts |& tee -a ${log}

echo "# STAR mapping ${pid} reads finished" |& tee -a ${log}

samtools index ${dir}/alignments/${pid}.Aligned.sortedByCoord.out.bam
samtools flagstat -@ ${nth} ${bam} > ${bam}.flagstat

echo "# Counting reads mapped to genes in ${pid}" |& tee -a ${log}
featureCounts -a ${gtf} \
	-o ${pid}/alignments/${pid}.featureCounts.tab \
	-T ${nth} \
	-F "GTF" \
	-t "CDS" \
	-p \
	${bam} |& tee -a ${pid}.featureCounts.out

echo "# Done ${pid}" |& tee -a ${log}