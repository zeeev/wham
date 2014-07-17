NSAMP=10
comma := ,
space :=
space +=
tbams =
bbams =
TIME_STAMP =
POS = 0
all: generateControl generateDeletion splitDeletion simulateReads align sam2bam sortBam indexBam cleanIbam defineTB timeStamp truth runWHAM

generateControl:
	for n in 1 2 3 4 5 6 7 8 9 10; do wgsim -N 100000 -1 100 -2 100 scaffold612.fa control-$$n.fasta.PE1.fastq control-$$n.fasta.PE2.fastq ; touch control-$$n.fasta ; done
generateDeletion:
	perl  generateSV.pl -n $(NSAMP) -o scaffold612.fa -t DELETION > target_deletion.multi.fa
splitDeletion:
	perl ~/tools/fasta_tool --split target_deletion.multi.fa	
simulateReads:
	find Sample-* | xargs -I {} -P $(NSAMP) wgsim -N 100000 -1 100 -2 100 {} {}.PE1.fastq {}.PE2.fastq
align:
	find *fasta | xargs -I {} -P $(NSAMP) ~/tools/bowtie2-2.0.2/bowtie2 -p 3 -x scaffold612.fa -1 {}.PE1.fastq -2 {}.PE2.fastq -S {}.sam
sam2bam:
	 find *sam | xargs -I {} -P $(NSAMP) ~/tools/samtools/samtools view -Sb {} -o {}.bam
sortBam: 
	find *bam | xargs -I {} -P $(NSAMP) ~/tools/samtools/samtools sort {} {}.sort
indexBam: 
	find *sort.bam | xargs -I {} -P $(NSAMP)  ~/tools/samtools/samtools index {}
cleanIbam:
	rm *sam.bam
cleanFastq:
	rm *fastq
cleanFasta:
	rm *fasta
cleanSam:
	rm *sam
cleanBam:
	rm *bai
	rm *bam

clean: cleanFastq cleanSam cleanFasta cleanBam 

defineTB:
	$(eval tbams = $(subst  bam$(space)S,bam$(comma)S, $(wildcard Sample*bam) ))
	echo $(tbams)	
	$(eval bbams = $(subst  bam$(space)c,bam$(comma)c, $(wildcard control*bam)))
	echo $(bbams)
runWHAM:
	~/tools/wham/wham/bin/WHAM-BAM -t $(tbams) -b $(bbams) > DELETION-$(TIME_STAMP)_$(POS).wham.out
timeStamp:
	$(eval TIME_STAMP=$(shell date +"%s-%d-%m"))
truth:
	$(eval POS=$(shell find *DELETION*bam | awk 'BEGIN{FS ="-"}{print $$4}'  | head -n 1))
	echo $(POS)