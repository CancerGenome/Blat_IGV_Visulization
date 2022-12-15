## Purpose

Visulization for the locally run Blat Result (PSL Format) with IGV report.

## Requirement

+ Best on Linux Server;
+ Have conda installed, https://docs.conda.io/en/latest/
+ Locally installed igv-report, https://github.com/igvteam/igv-reports
+ Locally installed Blat. https://anaconda.org/bioconda/blat
+ Locally installed tabix/bgzip. https://anaconda.org/bioconda/tabix
+ Locally installed samtools. https://anaconda.org/bioconda/samtools

## Procedures
+ Step 1: Download your own genome reference. Replace the REFERENCE in Run.sh with your reference location. 

+ Step 2: Repalce the REFERENCE and REFERENCE_INDEX in Run.sh with your reference and index file (*.fai). http://www.htslib.org/doc/samtools-index.html

+ Step 3: The annotation file used here is GenCode V41 and MANE V1.0 (db folder). All are downloaded from UCSC genome table browser (hg38). Go to step 5 if you are using hg38 and do not want to update your annotation file. Here, I switched the gene symbol and ENST in GenCode file to display ENST ID first, rather than display gene symbol.

+ Step 4 (optional): If you want to use your own annotation file, please download and update them in the config/track.json. Check this https://github.com/igvteam/igv.js/wiki/Tracks-2.0 if you are not familiar with IGV-report json file. 

+ Step 5: conda activate igvreports

+ Step 6: ./Run.sh YOUR_OWN_FASTA

+ Step 7: Check your result IGV.html

## Example
+ Input: test.fa, randomly choose two exons
+ Output: IGV.html

![image](https://user-images.githubusercontent.com/10337703/206536243-2efbb97f-8689-4836-9c8e-64de6ad2d0d8.png)

