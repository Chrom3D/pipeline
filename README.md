# Pipeline to create the input [Model Setup File (GTrack)]() from the Hi-C data
#### In this pipeline we will explain how to use Hi-C and Lamin ChIP-Seq data (optional) to create the Model Setup File (GTrack) as an input to the Chrom3D to model genome in 3 dimensions. 

## Preliminary steps and requirements

## To install Chrom3D and its dependencies
#### Please visit this [link](https://github.com/Chrom3D/Chrom3D#installation-instructions)


## Requirements to run this pipeline:

* OS: Linux or Mac
* Python 2.7, [pybedtools](https://daler.github.io/pybedtools/main.html) and [statsmodels](https://pypi.python.org/pypi/statsmodels)
* [bedtools](http://bedtools.readthedocs.io/en/latest/content/installation.html)
* [NCHG](http://folk.uio.no/jonaspau/hic/NCHG_hic.zip)
* [Chimera](https://www.cgl.ucsf.edu/chimera/download.html)

## To install NCHG
```
curl -O http://folk.uio.no/jonaspau/hic/NCHG_hic.zip
unzip NCHG_hic.zip
cd NCHG_hic
make
export PATH=$PATH:${PWD}
```
## Hi-C data analysis

* We recommend to use the [HiC-Pro](https://github.com/nservant/HiC-Pro) for the initial processes of Hi-C data (mapping and creating Hi-C matrix) as this output (Hi-C matrix in COO format) will serve as the input to this pipeline
* TAD callers such as [Arrowhead](https://github.com/theaidenlab/juicer/wiki/Arrowhead), [TopDom](http://zhoulab.usc.edu/TopDom/) or [Armatus](https://github.com/kingsfordgroup/armatus) can be used to define TADs in this pipeline.

## Assumptions:
We assume that the user already have following data to run the pipeline:

* Raw Hi-C matrix (using HiC-Pro)
	- We recommend 50kb resolution for the intra-chromosomal contact matrices and 1mb resolution for the inter-chromosomal contact matrices 
* Topologically Associating Domains (TADs) BED file (created using any TAD callers)
	- In this pipeline, we assume that the user has defined TADs using the Arrowhead algorithm 
* Necessary [scripts](https://github.com/Chrom3D/preprocess_scripts) (<-- the link to Download)
* Unmappable blacklist BED file (User can find the unmappable_blacklist.bed for the hg19 assembly in the scripts folder)
* Chromosomes size file (hg19.chrom.sizes.sorted in the scripts folder)

**Optional:**

* Lamin associated domains (LAD) BED file

---
---


## Main steps

### Converting HiC-Pro output:

Substitute "sample" with your sample name.

The sample\_50000.matrix and sample\_50000\_abs.bed can be found in the HiC-Pro output folder.

Example: USER\_WD/hic_results/matrix/sample/raw/50000

We recommend to copy these files to the "preprocess_scripts" folder.

`
python conv_hicpro_mat.py sample_50000.matrix sample_50000_abs.bed > sample_50000.intermediate.bedpe
`

The script "conv\_hicpro\_mat.py" converts the HiC-Pro COO matrix output to the intermediate BEDPE format file.  

#### Create intra-chromosomal contact matrices

This step creates intra-chromosomal RAW observed contact matrices for each chromosome using the BEDPE file created above at 50kb resolution.   

```
mkdir intra_chr_RAWobserved
 
# Example (chr1):
awk '{if($1=="chr1" && $4=="chr1") print $2 "\t" $5 "\t" $7}' sample_50000.intermediate.bedpe > intra_chr_RAWobserved/chr1_50kb.RAWobserved

#Run the following automated script
bash make_intrachr_rawObserved.sh chrom.sizes.sorted sample_50000.intermediate.bedpe

DESCRIPTION:
bash make_intrachr_rawObserved.sh [chromSizeFile] [intermediateBedpe]

[chromSizeFile] A text file containing the chromosome sizes (sorted numerically)
[intermediateBedpe] An intermediate bedpe created using "conv_hicpro_mat.py"
```

#### Create inter-chromosomal contact matrices

This step creates inter-chromosomal RAW observed contact matrices for each pair of chromosomes at 1mb resolution.

```
python conv_hicpro_mat.py sample_1000000.matrix sample_1000000_abs.bed > sample_1000000.intermediate.bedpe

mkdir inter_chr_RAWobserved 

# Example for inter-chromosome contact between chr1 and chr2:
awk '{if($1=="chr1" && $4=="chr2") print $2 "\t" $5 "\t" $7}' sample_1000000.intermediate.bedpe > inter_chr_RAWobserved/chr1_2_1mb.RAWobserved

#Run the following automated script
bash make_interchr_rawObserved.sh chrom.sizes.sorted sample_1000000.intermediate.bedpe

DESCRIPTION:
bash make_interchr_rawObserved.sh [chromSizeFile] [intermediateBedpe]

[chromSizeFile] A text file containing the chromosome sizes (sorted numerically)
[intermediateBedpe] An intermediate bedpe created using "conv_hicpro_mat.py"
```
---

### Aggregation of Hi-C contact counts for all pairs of TADs

**(i) To create the BED file specifying the genomic positions of the beads, Arrowhead domains (TADs) need to be merged, and gaps between them need to be filled.**

```
bash arrowhead_to_domains.sh sample_Arrowhead_domainlist.txt chrom.sizes.sorted

DESCRIPTION:
bash arrowhead_to_domains.sh [arrowheadFile] [chromSizeFile]

[arrowheadFIle] A text file listing the Arrowhead domains
[chromSizeFile] A text file containing the chromosome sizes (sorted numerically)
```

**(ii) Concatenate all the .domains to use in a later step**

```
cat *.chr*.domains > sample_Arrowhead_domainlist.domains
```

**(iii) Compute intra-chromosomal interaction counts between TADs**

```
mkdir intrachr_bedpe
bash make_NCHG_input.sh sample_Arrowhead_domainlist.chr19.domains intra_chr_RAWobserved/chr19_50kb.RAWobserved chr19 > intrachr_bedpe/chr19_50kb.domains.RAW.bedpe
```

Run the following script to run the above command for all chromosomes

```
bash intrachr_NCHG_input_auto.sh sample_Arrowhead_domainlist chrom.sizes.sorted 50kb

DESCRIPTION:
bash intrachr_NCHG_input_auto.sh [domainBase] [chromSizeFile] [resolution]

[domainBase] Basename of the domain files
[chromSizeFile] Text file containing the chromosome sizes (sorted numerically)
[resolution] Resolution of the interaction matrix as given in the matrix filename (eg. 50kb or 1mb)

```
**(iv) Concatenate all the bedpe files**

`cat intrachr_bedpe/chr*.bedpe > intrachr_bedpe/sample_50kb.domain.RAW.bedpe`

**(v) Remove domains that contain centromeres from the BEDPE file**

```
curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz" | gunzip -c | grep acen | pairToBed -a intrachr_bedpe/sample_50kb.domain.RAW.bedpe -b stdin -type neither > intrachr_bedpe/sample_50kb.domain.RAW.no_cen.bedpe
```

---

###Identification of significant intra-chromosomal interaction

This step uses a non-central hypergeometric test to calclate P-value and odds ratio for each TAD-TAD interactions. Then, the significant TAD-TAD interactions are filtered using FDR and odds ratio.

**(i) Calculate the P-value and odds ratio for each pair of TADs**

```
NCHG -m 50000 -p intrachr_bedpe/sample_50kb.domain.RAW.no_cen.bedpe > sample_50kb.domain.RAW.no_cen.NCHG.out

DESCRIPTION:
NCHG -m [minDistance] -p [NULL] [inputFile]

-m Minimum genomic distance (in bp) allowed between domains, below which interactions are excluded
-p Print output to stdout
[inputFile] Hi-C contact count data in BEDPE format
```

**(ii) Calculate FDR and filter significant interactions**

```
python NCHG_fdr_oddratio_calc.py sample_50kb.domain.RAW.no_cen.NCHG.out fdr_bh 2 0.01 > sample_50kb.domain.RAW.no_cen.NCHG.sig

DESCRIPTION:
python NCHG_fdr_oddratio_calc.py [inputFile] [testMethod] [cutoff] [threshold]

[inputFile] Input filename (NCHG output; BEDPE format)
[testMethod] Multiple hypothesis testing method (bonferroni, sidak, holm-sidak, holm, simes-hochberg, hommel, fdr_bh, fdr_by, fdr_tsbh, fdr_tsbky)
[cutoff] Odds ratio cutoff value
[threshold] Significance threshold value after multiple testing correction
```

---

### Create the Model Setup File [(GTrack)]()

A Model Setup File (GTrack) is the input file to the Chrom3D which specifies the genomic coordinates, unique id and constraints for the beads. Users can add HiC constraints (TAD-TAD interactions in a chromosome and between chromosomes) using the "edge" column and if the lamin ChIP-Seq data is available, the LAD constraints (beads constrain towards nuclear periphery) using the "periphery" column can be added. A simple illustration of a GTrack file could be found at [the end of the this tutorial](). This step involves creating and adding constraints using HiC data and lamin ChIP-Seq data to a GTrack file. 

**(i) Create GTrack using significant interactions**


```
bash make_gtrack.sh sample_50kb.domain.RAW.no_cen.NCHG.sig sample_Arrowhead_domainlist.domains sample_intra_chromosome.gtrack

DESCRIPTION:
bash make_gtrack.sh [sigFile] [domainFile] [outputFile]

[sigFile] Intra-chromosome significant interactions file
[domainFile] Domain file
[outputFile] Output GTrack file (Model Setup File)

```

**(ii) Add LAD information to the GTrack (OPTIONAL)**

This step will add the periphery column to the GTrack file (specifying beads with periphery constrains). Please refer the [illustration]().    

```
bash make_gtrack_incl_lad.sh sample_intra_chromosome.gtrack sample_LAD.bed sample_intra_chromosome_w_LADs.gtrack

DESCRIPTION:
bash make_gtrack_incl_lad.sh [inputFile] [ladFile] [outputFile]

[inputFile] Input GTrack file (without LAD information)
[ladFile] LAD BED file to be added to GTrack file
[outputFile] Output GTrack file with LAD information added 
```

**(iii) Prepare inter-chromosomal Hi-C interaction counts**

```
bash ./interchr_NCHG_input_auto.sh chrom.sizes.sorted unmappable_blacklist.bed 1mb > sample_1mb_inter.bedpe

DESCRIPTION:
bash interchr_NCHG_input_auto.sh [chromSizeFile] [blackList] [resolution]

[chromSizeFile] Text file containing the chromosome sizes (sorted numerically)
[blackList] BED file containing positions of blacklisted regions
[resolution] Resolution of the interaction matrix as given in the matrix filename (e.g. 50kb or 1mb)
```

**(iv) Call significant inter-chromosomal interactions**

```
NCHG -i -p sample_1mb_inter.bedpe > sample_1mb_inter_chr.NCHG.out

DESCRIPTION:
NCHG -i [NULL] -p [NULL] [inputFile] > [outputFile]

-i Instructs NCHG to use inter-chromosomal interactions 
-p Instructs NCHG to print output to stdout [inputFile] Hi-C contact count data (BEDPE format)
[outputFile] File containing P-values for each TAD pair
```

**(v) Calculate FDR and filter significant interactions**

```
python NCHG_fdr_oddratio_calc.py sample_1mb_inter_chr.NCHG.out fdr_bh 2 0.01 > sample_1mb_inter_chr.NCHG.sig

DESCRIPTION:
python NCHG_fdr_oddratio_calc.py [inputFile] [testMethod] [cutoff] [threshold]

[inputFile] Input filename (NCHG output; BEDPE format)
[testMethod] Multiple hypothesis testing method (bonferroni, sidak, holm-sidak, holm, simes-hochberg, hommel, fdr_bh, fdr_by, fdr_tsbh, fdr_tsbky)
[cutoff] Odds ratio cutoff value
[threshold] Significance threshold value after multiple testing correction

```

**(vi) Add significant inter-chromosomal interaction information to the GTrack**

This step adds siginificant inter-chromosomal interaction to the edge column of the Gtrack file. 

```

bash add_inter_chrom_beads.sh sample_intra_chromosome_w_LADs.gtrack sample_1mb_inter_chr.NCHG.sig sample_inter_intra_chr_w_LADs.gtrack

DESCRIPTION:
bash add_inter_chrom_beads.sh [inputFile] [sigFile] [outFile]

[inputFile] Input GTrack file
[sigFile] Inter-chromosome significant interaction file
[outFile] Output GTrack file containing significant inter-chromosomal interactions
```

**(vii) Modify the Model Setup File to make a diploid model**


```
python make_diploid_gtrack.py sample_inter_intra_chr_w_LADs.gtrack > sample_inter_intra_chr_w_LADs.diploid.gtrack

DESCRIPTION:
python make_diploid_gtrack.py [inputFile] > [outputFile]

[inputFile] Input GTrack file with constraints specified for single chromosomes
[outputFile] Output GTrack file with constraints for each chromosome copy
```

---

### Run Chrom3D to generate 3D genome models

```
Chrom3D -y 0.15 -r 5.0 -n 2000000 -o sample_inter_intra_chr_w_LADs.diploid.cmm sample_inter_intra_chr_w_LADs.diploid.gtrack

DESCRIPTION:
Chrom3D -y [scale] -r [radius] -n [iterations] -o [outputFile] [setupFile]

-y Scale total volume of the model beads relative to the volume of the nucleus
-r Radius of the nucleus in micrometers
-n Number of iterations
-o Output filename
[setupFile] Model Setup File in GTrack format
```
---

### Visualisation of 3D genome models using chimera

`
chimera sample_inter_intra_chr_w_LADs.diploid.cmm
`

**Highlighting LAD-containing beads**


```
python color_beads.py sample_inter_intra_chr_w_LADs.diploid.cmm lad_ad04.ids 0,0,255 override > sample_inter_intra_chr_w_LADs.diploid.visLAD.cmm

DESCRIPTION:
python color_beads.py [inputFile] [beadIdFile] [color] [colorScheme]

[inputFile] Input CMM file
[beadIdFile] File containing selected bead ids
[color] Comma-separated RGB value
[colorScheme] ‘blend’ or ‘override’ color of the beads specified
 
```

**Highlighting HiC constrained beads and visualisation**

* LAD-containing beads - blue
* HiC constrained beads - red
* Beads with both LAD and HiC beads - purple

```
python color_beads.py sample_inter_intra_chr_w_LADs.diploid.visLAD.cmm TAD_interaction.ids 255,0,0 blend > sample_inter_intra_chr_w_LADs.diploid.visLADCons.cmm

chimera sample_inter_intra_chr_w_LADs.diploid.visLADCons.cmm
```


##Simple illustration of a Model Setup File 


![Example GTrack](http://folk.uio.no/tmali/git_ups/gtrack_illust.png)


