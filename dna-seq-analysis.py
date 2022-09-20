##Command Line Applications for NGS
#Analysis of Illumina reads and ENSEMBL data involves QC assessment, filtering, mapping, variant calling and data visualisation.

source /Users/<username>/opt/miniconda3/bin/activate
conda --version # check

## Move & download data    
# download ENSEMBL data from an instance / AWS-Amazon Web Services
pwd
curl -O ftp://ftp.ensemblgenomes.org/pub/release-37/bacteria/species_EnsemblBacteria.txt
which scp
ls

# or use scp <file to move> <new location>
scp species_EnsemblBacteria.txt #<[user@IP_address_login_credentials]AWS instance>


## Organise and search files
mkdir -p experiments/{scripts,fastqc/fastqc_untrimmed} | mkdir backups
curl -L https://osf.io/vuk5y/download -o experiments/shell_data.zip | \
unzip experiments/shell_data.zip | \
mv -v shell_data/*/* experiments/ | rm -r shell_data | \
cp -R experiments/ backups | \
ls backups # check


# simply check file contents
grep NNNNNNNNNNC "SRR<inputfile>.fastq"
grep -B1 GNATNACCACTTCC "SRR<inputfile>.fastq" # '-B' returns identifiers (name; matches)
grep -A2 NNNNNNNNNNATAT "SRR<inputfile>.fastq" # '-A' returns quality scores associated with reads (plus one (5’) or two (3’) lines flanking each match) plus identifier lines

wc "SRR<inputfile>.fastq" # check counts for words, lines & characters

cd experiments

for i in *.fastq; do \
head -n 2 ${i} >> seq_info.txt; echo "Showing $i"; done


## Scripts for new projects
nano scripts/new_projects.sh
#!/bin/bash
  # 2022_xx_xx
  # Directories created for $i project
  for i in *.fastq; do name=$(basename ${i} .fastq); mkdir -p project_$name/{docs,data,results}; echo "Folders created for project_$name."; \
      # OR mkdir project_$name | mkdir project_$name/docs | mkdir project_$name/data | mkdir project_$name/results | mkdir project_$name/scripts
  ls -l new_projects.sh; \
  chmod +x new_projects.sh; \
  done
  bash new_projects.sh

./new_projects.sh # run script



# script for bad reads
nano scripts/bad_reads.sh
#!/bin/bash
for f in *.fastq; do \
name=$(basename ${i} .fastq); \
cp ${i} project_$name/data/; \
grep -B1 -A2 -h NNNNNNNNNN *.fastq | grep -v '^--' > project_$name/data/bad_reads_$name.txt; echo "Script ${name} finished."; \
chmod +x scripts/bad_reads.sh; \
ls -l; done

./scripts/bad_reads.sh # run script to write .txt file


cat bad_reads_script.sh # prints all contents in directory; reads/ concatenates files

wc -l bad_reads_script.sh # check: 264 characters (264/3=88 integer-file is ok)
ls -l bad_reads_script.sh # lists permissions
chmod +x bad_reads_script.sh # adds executable permission
# ./bad_reads_script.sh # ./ look in this directory for the program




## FastQC quality assessment
ls -l -h
conda install -c bioconda/label/cf201901 fastqc | fastqc -help
conda install -c bioconda seqkit #fasta/q file manipulation package


# create fastQC script
nano scripts/fastQC.sh
#!/bin/bash

echo "QC Analysis running...."
fastqc -o fastqc/fastqc_untrimmed */data/*.fastq
set -e  # exit if error occurs
ls fastqc/fastqc_untrimmed

for z in fastqc/fastqc_untrimmed/*.zip; do \
  echo "Decompressing files..."; \
  unzip ${z} -d fastqc/fastqc_untrimmed/; break;   #unzip qc files
done
set -e

ls -l fastqc/fastqc_untrimmed/
chmod +x fastqc/fastqc_untrimmed/* # make all files executable
ls -l fastqc/fastqc_untrimmed/

cat fastqc/fastqc_untrimmed/*_fastqc/summary.txt > fastqc/fastqc_untrimmed/all_fastqc_summaries.txt # concatenates all summary files

ls -l scripts/fastQC.sh
chmod +x scripts/fastQC.sh
./scripts/fastQC.sh # run script ; also runs script 'bash scripts/fastQC.sh'

open fastqc/fastqc_untrimmed/*.html # opens QC results
less fastqc/fastqc_untrimmed/all_fastqc_summaries.txt # preview; 'q' to exit
grep FAIL fastqc/fastqc_untrimmed/all_fastqc_summaries.txt # searches for failed tests



## Trim & Filter

# trimmomatic
conda install -c bioconda/label/cf201901 trimmomatic | trimmomatic help

# observe any adapter seqs by looking at fastQC adapter content & overrepresented seqs for hints
# open fastqc/fastqc_untrimmed/*.html

# OR download/ copy adapters from trimmomatic package for use in working directory
mkdir adapters/ | cp -R /Users/<username>/opt/miniconda3/pkgs/trimmomatic-0.38-1/share/trimmomatic-0.38-1/adapters/ ./adapters


# create script to trim and re-run fastQC to observe qual scores (esp. per base)
                              # selected adapters dependent on sample library prep
nano scripts/fastQC_trim.sh
#!/bin/bash
for f in */*.fastq; do \
name=$(basename ${f} .fastq); \ #re-define variable
trimmomatic SE -threads 4 -phred33 project_$name/data/$f project_$name/data/$name\_trim.fastq ILLUMINACLIP:adapters/TruSeq3-SE.fa:2:40:15 SLIDINGWINDOW:4:0 MINLEN:25; \
set -e; \
echo "Finished trimming $f file(s)."; \
ls -l -h project_$name/data/$name\_trim.fastq; \
done

echo "QC Analysis running...."
fastqc -o fastqc/ project_$name/data/$name\_trim.fastq
set -e
ls fastqc

for z2 in fastqc/*.zip; do \
  echo "Decompressing files..."; \
  unzip ${z2} -d fastqc/; done
set -e

ls -l scripts/fastQC_trim.sh
chmod +x scripts/fastQC_trim.sh
set -e
ls -l scripts/fastQC_trim.sh
# q, quit


./scripts/fastQC_trim.sh



## Variant Calling & Alignment
# download and decompress ref. file

for i in *.fastq; do \
name=$(basename ${i} .fastq); \
curl -L -o project_$name/data/ecoli_rel606.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz; \
mkdir -pv project_$name/results/{sam,\bam,\bcf,\vcf,\igv_files}; done


# make script align_variant_calling

nano scripts/align_variant_calling.sh
#!/bin/bash
for gz in */data/*.gz; do \
echo "Decompressing $gz files..."; \
  for f in */data/*_trim.fastq; do name=$(basename ${f} _trim.fastq); \  #re-define variable
ref=$(basename ${gz} .gz); \
gunzip -c ${gz} > project_$name/data/${ref}; \
 done 


mv project_$name/data/ecoli_rel606.fasta project_$name/data/ref_ecoli_rel606.fasta #rename reference genome
head -n 1 data/ref_ecoli_rel606.fasta # view genome name via first line

tar xvf data/trimmed_fastq.sub.tar.gz -C data/ #untar/ decompress trimmed data
  
  #---OR TRIM USING TRIMMOATIC!!---
  
trimmomatic #check

seqkit stats data/*f{asta,astq} # observes summary about fasta and fastq files

for R1 in data/*_1.fastq; do \
R2=${R1%_1.fastq}_2.fastq; \    # '%' is a modula; remainder when a is divided by b (a%b); creates new variable by replacing old variable - newvariablename=${oldvariablename//oldtext/newtext}
  trimmomatic PE -threads 4 -phred33 $R1 $R2 "${R1%.*}_trim_paired.fastq.gz" "${R1%.*}_trim_unpaired.fastq.gz" "${R2%.*}_trim_paired.fastq.gz" "${R2%.*}_trim_unpaired.fastq.gz" \  # creates two files of paired and unpaired trim data for each paired-end
  ILLUMINACLIP:"/Users/<username>/opt/miniconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa":2:40:15 MINLEN:25; \  #select adapter
      echo "Finished trimming $R1 $R2 file(s)."; \
done

# index ref. file w BWA, align w BWA-MEM, view, sort then create bam file
    # info from a paired read shouldn't be combined (using samtools sort, then merge) w unpaired reads;
  # paired read = measurements from same fragment measured (sequenced) twice
  # unpaired reads = measure different fragments
bwa index data/ref_ecoli_rel606.fasta #creates faidx/ baidx ext. files
  
  for trim1 in data/*\1_trim_paired.fastq.gz; do \  # for all trim1 data in data dir ending in '\1_trim_paired.fastq.gz'
  trim2=${trim1%\1_trim_paired.fastq.gz}\2_trim_paired.fastq.gz; \  # creates variable for second paired trimmed data, trim2
  trim_name=$(basename ${trim1} \1_trim_paired.fastq.gz)trim_paired; \  # creates new variable trim_name from trimmed file name by replacing "\1_trim_paired.fastq.gz" with "_trim_paired"
  bwa mem data/ref_ecoli_rel606.fasta $trim1 $trim2 > results/sam/$trim_name\_aligned.sam; \  #align both paired reads to indexed reference
  done
  
  echo $trim1 #check saved variable
  echo $trim2 #check
  echo $trim_name #check
  
  conda install -c bioconda/label/cf201901 samtools
  conda install -c bioconda/label/cf201901 bcftools
  
  samtools -h
  # error message, so create environment? see error doc..
  samtools help
  
  # view SAM header
  ls results/sam
  head -n 3 results/sam/$trim_name\_aligned.sam   # check header;    'q' to quit
  
  # convert aligned sam to bam, using -b flag
  cd results/sam
  
  # create for loop and new variable to parse in loop
  # then sort bam by coordinates;   can also be sorted by chr., coordinates, read name etc.
  for samfile in *.sam; do \
  samplename=$(basename ${samfile} _trim_paired_aligned.sam); \   #create new variable from old one by removing suffix from sam file using basename function
  samtools view -bh $samfile > "/Users/<username>/Library/CloudStorage/OneDrive-ImperialCollegeLondon/ImperialBits/PROJECTS (1)/BioInformatics/Command Line/alignment_variant_calling/ecoli_exp/results/bam/${samplename}.bam"; \  # -h flag includes header in bam
  cd ..; \
  samtools sort -o "/Users/<username>/Library/CloudStorage/OneDrive-ImperialCollegeLondon/ImperialBits/PROJECTS (1)/BioInformatics/Command Line/alignment_variant_calling/ecoli_exp/results/bam/${samplename}.sorted.bam" "/Users/ishittu/Library/CloudStorage/OneDrive-ImperialCollegeLondon/ImperialBits/PROJECTS (1)/BioInformatics/Command Line/alignment_variant_calling/ecoli_exp/results/bam/${samplename}.bam"; \
  done
  
  cd ..
  samtools view bam/${samplename}.sorted.bam | head -n 3 # check
  
  # observe info about bam file
  samtools flagstat bam/${samplename}.sorted.bam # flagstat counts the number of alignments for each FLAG type
  
  
  # Variant call by;
  # first calculating read coverage of genomic positions 
      # mpileup - uses mapping scores to evaluate variant calling and create genotype likelihoods (used to be samtools mpileup)
  cd ..
  bcftools mpileup -Ob -o results/bcf/${samplename}\_raw.bcf \
  -f data/ref_ecoli_rel606.fasta results/bam/${samplename}.sorted.bam
  
  
  bcftools head -n 5 results/bcf/${samplename}\_raw.bcf # check, ctrl+c to quit
  bcftools view -h results/bcf/${samplename}\_raw.bcf # check, shows header info
  
  
  # detect SNVs (single nt variants) using call command and converting bcf compressed file to vcf 
  bcftools call --ploidy 1 -m -v -o results/vcf/$samplename\_variants.vcf results/bcf/$samplename\_raw.bcf
   
  # results	- lists the values associated with those metrics in order
  vcfutils.pl varFilter results/vcf/$samplename\_variants.vcf > results/vcf/$samplename\_final_variants.vcf
  vcfutils.pl varFilter -p results/vcf/$samplename\_variants.vcf > results/vcf/$samplename\_final_variants.vcf # -p flags prints filtered variants
  # could also try bcftools filter
  
  # explore vcf file, -S flag incl. sample columns/descriptions
  less -S results/vcf/$samplename\_final_variants.vcf # check, q to quit
  
  
  grep -v -c "#" results/vcf/$samplename\_final_variants.vcf # counts no. of SNPs in vcf file (searches and counts for lines not matching #)
  # -c  counts no. of matched lines;   -v shows lines that don't match
  # grep -v "#" results/vcf/$samplename\_final_variants.vcf | wc -l #also works
  
  bcftools stats results/vcf/$samplename\_final_variants.vcf > results/vcf/$samplename\_variant_stats.vcf #no. of records = no of variants
  less results/vcf/$samplename\_variant_stats.vcf # check
  
  
  samtools index results/bam/${samplename}.sorted.bam # indexes bam file 
  
  # visualise IGV, Integrative Genomics Viewer    -   load ref. genome > load sorted bam file > repeat w/ final vcf
  #OR w/
  samtools tview results/bam/${samplename}.sorted.bam data/ref_ecoli_rel606.fasta # command uses simple alignment viewer based on ncurses library, called tview, for short indels
  # output: navigate with 'g'; type chr. name + position ie 'CP000819.1:4377265 is 4377265th base of chr. CP000819.1 - has a canonical nt A with a G variant
  # first line = genome coordinates of ref.
  # second line = ref.
  # third line = consensus seq from reads
  # . = a match to ref. ∴ sample's consensus matches the in most locations, if not reconsider ref.
  
  
  
  #--------Extra--------
  
  # a second file could be used for comparison uses
  conda install -c "bioconda/label/cf201901" bedtools | bedtools -h
  
  intersectBed -a results/vcf/$samplename\_variant_stats.vcf \
  -b <(string.of).comparisonfiles> -sorted | wc -l # counts no. of SNPs in common between two files
  
  
  
  # NanoPlot
  # Plotting tool for long read sequencing data and alignments creates:
  # a statistical summary;  plots & html summary file from bam
  conda install -c "bioconda/label/cf201901" nanoplot | nanoplot -h
  mkdir -pv results/nanoplot/summary | \
  nanoplot -t 6 --fastq data/SRR2584866_1_trim_paired.fastq.gz data/SRR2584866_2_trim_paired.fastq.gz -o results/nanoplot --maxlength 40000 --plots --legacy dot | \
  
  open results/nanoplot/NanoPlot-report.html
