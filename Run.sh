#!/bin/bash
if [ $# -le 0 ]
then echo "
#########################################################################
#      USAGE: Run.sh <FASTA>
#      Description: Given FASTA file,
#                   Output a IGV format html with annotation from GenCode V41 and MANE V1.0
#
#      Caveat: start igvreports first, conda activate igvreports
#
#      Author: Wang Yu
#########################################################################
"
exit
fi

#### Install the software if not 
####----------------------------
#conda create -n igvreports python=3.7.1
#conda activate igvreports
#pip install igv-reports

#### Processing the annotation files if not
####---------------------------------------

## Switch GenCode Gene Symol and ENST, to display ENST only

# less GeneCodeV41.gz | awk 'BEGIN{OFS="\t"}{a=$2;$2=$13;$13=a;print $0}' |~/miniconda3/bin/bgzip > GeneCodeV41.ENST.gz
# tabix -s3 -b5 -e6 GeneCodeV41.ENST.gz

## index MANE file
# tabix -p bed MANE_V1.gz

#conda activate igvreports

#### Step 1 Mapping, Select Rep by Score and Covert to Bam, Bed file
####-------------------------------------------------------
FASTA=$1
echo "Start Blat"
blat REFERENCE $FASTA Blat.psl
sed '1,5d' Blat.psl  | sort -k10V,10 -k1n,1r | awk '!a[$10]++' > Blat.rep.psl
bin/psl2samV2.pl Blat.rep.psl $FASTA | sort -k3V,3 -k4n,4 | samtools view -bt REFERENCE_INDEX - > Blat.rep.bam
samtools index Blat.rep.bam

less Blat.rep.psl |  awk '{print $14"\t"$16-10"\t"$17+10"\t"$10"_Score:"$1}' | sort -k1V,1 -k2n,2 -k3n,3  |grep -v chrUn |grep -v random > Blat.bed

#### Step 2 Run IGV report
####-------------------------------------------------------
echo "Create Report"
create_report Blat.bed --type junctions --genome hg38 --track-config  config/track.json --output IGV.html --ideogram db/cytoBandIdeo.txt

# clean
mv Blat* archive
mv $FASTA archive

#conda deactivate
