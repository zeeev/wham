# setup output folder
mkdir test_output
# run wham
../bin/wham  data/mut1_20-40k.bam data/mut1_20-40k.bam.bai 'vbrc|1485' > test_output/mut1.txt 2> test_output/mut1.err
../bin/wham  data/mut2_20-40k.bam data/mut2_20-40k.bam.bai 'vbrc|1485' > test_output/mut2.txt 2> test_output/mut2.err
../bin/wham  data/ref_20-40k.bam  data/ref_20-40k.bam.bai  'vbrc|1485' > test_output/ref.txt  2> test_output/ref.err

#cleanup
mv test_output test1

