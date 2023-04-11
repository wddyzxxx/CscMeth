<a name="DUDaV"></a>
# CscMeth - Guided Analysing Tutorial
<a name="WHhGS"></a>
#### Compiled: Arpil 10, 2023 			By Yuan Zhen
<a name="dHhke"></a>
# Prepare ALLC files for analysis
For this tutorial, we will be analyzing a simple methylation data of two cells, one normal cell and one colorectal cancer cell, i.e., `CRC04_NC_302.sort.rmdup.bam`and`CRC04_PT1_589.sort.rmdup.bam`.<br />We start by transforming them into ALLC format files, which consists of seven columns storing different information:

| index | column name | example | note |
| --- | --- | --- | --- |
| 1 | chromosome | chr12 | The same as genome FASTA |
| 2 | position | 18283342 | 1-based |
| 3 | strand | + | either + or - |
| 4 | sequence context | CGT | can be more than 3 bases, used to determine mC type |
| 5 | mc | 1 | count of reads supporting methylation |
| 6 | cov | 2 | read coverage, cov >= mc |
| 7 | methylated | 1 | indicator of significant methylation (1 if no test is performed) |

We first use `cscmeth bam2allc`to finish the transforming part. For convenience, we can just omit the partient information.
```
cscmeth bam2allc --bam_path ../4.0_40-bam/CRC04_NC_302.sort.rmdup.bam \
--reference_fasta /data/yuanzhen/database/0.Genome/hg19/hg19.fa \
--output_path NC_302

cscmeth bam2allc --bam_path ../4.0_40-bam/CRC04_PT1_589.sort.rmdup.bam \
--reference_fasta /data/yuanzhen/database/0.Genome/hg19/hg19.fa \
--output_path PT_589

###After this, we have four new generated files, two index files and two ALLC files:
PT_589.gz.tbi
PT_589.gz
NC_302.gz.tbi
NC_302.gz
```

<a name="nNgBo"></a>
# Standard pre-processing workflow
<a name="W3RSJ"></a>
## extract
Let's assume we are only interested in CpG sites, i.e., cytosines within NCG contexts. Then we can use `cscmeth extract` to do extraction job.
```
cscmeth extract --allc_path NC_302.gz \
--output_prefix NC_302 \
--mc_contexts CGN \
--chrom_size hg19

extract --allc_path PT_589.gz  \
--output_prefix PT_589 \
--mc_contexts CGN \
--chrom_size hg19

###And we will have two ALLC files merely with cytosines under CGN context. 
PT_589.CGN-Both.allc.tsv.gz.tbi
PT_589.CGN-Both.allc.tsv.gz
NC_302.CGN-Both.allc.tsv.gz.tbi
NC_302.CGN-Both.allc.tsv.gz
```
<a name="xo04X"></a>
## merge
Usually there are a lot of single cells belonging to onetype, so it's necessary to merge them for later DMR analysis with `cscmeth merge`command. But now since we just  take two cells as an simple example, let's skip it for now.
<a name="GhlJb"></a>
## intersect
As there will definitely exist cytosines which only one of them captured, we have to do a intersecting operation to make sure our following analysis are based on common CpG sites.
```
cscmeth intersect-allc --allc_file1 NC_302.CGN-Both.allc.tsv.gzz \
--allc_file2 PT_589.CGN-Both.allc.tsv.gz \
--prefix_file1 NC_302 \
--prefix_file2 PT_589

###And we have these two files:
PT_589_intersected.gz
NC_302_intersected.gz
```
Above are the necessary pre-processing steps to prepare intersected ALLC files for analysis.

<a name="VNvWG"></a>
# DMR analysis (premise for enrichment analysis)
Next, for the following Differentially Methylated Regions and genomic elements enrichment analysis, we take a windowing-based calculation method. In short, we group cytosines into fixed and continous windows and then perform paried t-test on cytosine methylation value within the same window.
<a name="dQ0bQ"></a>
## window
```
cscmeth window --allc_files 302_intersected.gz 589_intersected.gz \
--step 500 \
--genome_size hg19

### This process will create a corresponding directory to store your windowing 
### information. In this example, it is 500_step_splitted_results. And it includes 
### these two files:\
Methylation_info_302_intersected.gz
Methylation_info_589_intersected.gz
```
<a name="b8a7N"></a>
## perform test
```
cscmeth test --input1 Methylation_info_302_intersected.gz \
--input2 Methylation_info_589_intersected.gz \
--prefix NC_302-PT_589

### This is probabaly the most time-consuming step throughout the workflow. 
### After it's done, you will see four files, two bed files for future enrichment 
### analysis, two test results of both the background and customized filtered DMRs:
NC_302-PT_589_DMRs-Input.bed
NC_302-PT_589_DMRs.gz
NC_302-PT_589_Tested-Background.bed
NC_302-PT_589_Tested.gz
```
<a name="gE5lo"></a>
# Enrichment analysis
<a name="uxpaK"></a>
## annotate
To understand which elements the DMRs are enriched in the genome, we can first annotate them and then visualize enrichment results. But NOTE here we use `Tested`file as the input rather than another `DMRs` file because the calculation is based on hypergeometric distrbution.<br />And as for elements, there are ready ones in this package, for example in hg19 version:
```
hg19.Alu_raw_4col_sorted.bed
hg19.ERVL-MaLR_raw_4col_sorted.bed
hg19.LINE_raw_4col_sorted.bed
hg19.LTR_raw_4col_sorted.bed
hg19.MIR_raw_4col_sorted.bed
hg19.promoter_raw_4col_sorted.bed
hg19.SINE_raw_4col_sorted.bed
```
And just type:
```
### If there is only one element after --elements, then quotation mark won't be necessary.

cscmeth annotate --tested_file NC_302-PT_589_Tested.gz \
--elements 'Alu SINE LINE' \
--anno_dir hg19



### This step also creates a special directory 302vs589_Tested_annotation to store its 
### annotaion results for calcaulating p-values and plotting.

### command `ls` results
NC_302-PT_589_Tested.gz_Alu
NC_302-PT_589_Tested.gz_LINE
NC_302-PT_589_Tested.gz_SINE
```
<a name="FWNJq"></a>
## ※Plot elements enrichment result
```
### Simply specify the directory just created.

cscmeth elements --anno_dir 302vs589_Tested_annotation

### And you will have three new files, one enrichment input data, one enrichment result
### and one enrichment plot.

### However since we only use two cells as input, so actully the result shows no 
### significance. But out of the purpose of illustration, here shows another picture:
```
![X43WKIIM@22FYAK{G_LYBU3.png](https://cdn.nlark.com/yuque/0/2023/png/32598292/1681198619279-fc80232a-3f08-4cea-99dc-05312a3187ce.png#averageHue=%23f9f9f9&clientId=u4ba62734-d942-4&from=paste&height=833&id=u64146956&name=X43WKIIM%4022FYAK%7BG_LYBU3.png&originHeight=833&originWidth=1794&originalType=binary&ratio=1&rotation=0&showTitle=false&size=22822&status=done&style=none&taskId=u28c59439-9235-4b69-8296-a5e71614ce6&title=&width=1794)
<a name="A32Tb"></a>
## functional enrichment
To associate regions with particular biological processes, we incoporate rGREAT R package to satisfy this need. In this part we assume users are interested in DMRs and what biological activities they are involved, so we use `302vs589_DMRs-Input.bed` as test regions, and `302vs589_Tested-Background.bed` as the background. But users will have to specify another two arguments, namely, `--gene_set` and`--tss` . For short here, `gene_set` includes GO gene sets and gene sets collections from MSigDB and `tss`supports many formats, please use `cscmeth functional -h` for detail information.
```
### Here is when the two bed files generated during DMR analysis comes in handy.

cscmeth functional --test 302vs589_DMRs-Input.bed \
--background 302vs589_Tested-Background.bed \
--gene_set BP \
--tss hg19

### And three files will be produced:
302vs589_DMRs-Input_GO:BP_txdb:hg19.rds
302vs589_DMRs-Input_GO:BP_txdb:hg19.txt
Volcano_Plot_302vs589_DMRs-Input_GO:BP_txdb:hg19.pdf
```
`302vs589_DMRs-Input_GO:BP_txdb:hg19.rds` is for following step of associating regions to genes;<br />`302vs589_DMRs-Input_GO:BP_txdb:hg19.txt` includes GO:BP results file;<br />`Volcano_Plot_302vs589_DMRs-Input_GO:BP_txdb:hg19.pdf` visualizes GO:BP results.<br />![image.png](https://cdn.nlark.com/yuque/0/2023/png/32598292/1681199839192-b9a848ef-493d-4f62-bfed-b01df3188264.png#averageHue=%23fcfbfb&clientId=u4ba62734-d942-4&from=paste&height=998&id=ue65feb10&name=image.png&originHeight=998&originWidth=1050&originalType=binary&ratio=1&rotation=0&showTitle=false&size=140762&status=done&style=none&taskId=u384372da-407e-4481-bc85-9fac5b3cd21&title=&width=1050)
<a name="tLvmd"></a>
## region2gene
If we are interested in connecting the input regions with nearby genes, that's when `region2gene` should be used:
```
cscmeth region2gene --rds 302vs589_DMRs-Input_GO:BP_txdb:hg19_.rds

### You will get such two files:
RegionGeneAssociation_302vs589_DMRs-Input_GO:BP_txdb:hg19.csv
RegionGeneAssociation_Plot_302vs589_DMRs-Input_GO:BP_txdb:hg19.pdf
```
The table looks like this:

| seqnames | start | end | width | strand | annotated_genes | dist_to_TSS |
| --- | --- | --- | --- | --- | --- | --- |
| chr1 | 567002 | 567500 | 499 | * | OR4F5;SAMD11 | 497911;-293030 |
| chr1 | 1841002 | 1841500 | 499 | * | GNB1;CALML6 | -18476;-4766 |
| chr1 | 2055502 | 2056000 | 499 | * | PRKCZ;FAAP20 | 73593;83172 |
| chr1 | 2160502 | 2161000 | 499 | * | SKI | 368 |
| chr1 | 2371502 | 2372000 | 499 | * | PEX10;PLCH2 | -27492;-35754 |
| chr1 | 2398002 | 2398500 | 499 | * | PEX10;PLCH2 | -53992;-9254 |
| chr1 | 2589502 | 2590000 | 499 | * | MMEL1;TTC34 | -25021;116230 |
| chr1 | 2603002 | 2603500 | 499 | * | MMEL1;TTC34 | -38521;102730 |
| chr1 | 2689502 | 2690000 | 499 | * | MMEL1;TTC34 | -125021;16230 |
| chr1 | 2989002 | 2989500 | 499 | * | PRDM16;ARHGEF16 | 3260;-381647 |
| chr1 | 3021502 | 3022000 | 499 | * | PRDM16;ARHGEF16 | 35760;-349147 |
| chr1 | 3027002 | 3027500 | 499 | * | PRDM16;ARHGEF16 | 41260;-343647 |
| chr1 | 3030502 | 3031000 | 499 | * | PRDM16;ARHGEF16 | 44760;-340147 |
| chr1 | 3117002 | 3117500 | 499 | * | PRDM16;ARHGEF16 | 131260;-253647 |
| chr1 | 3145002 | 3145500 | 499 | * | PRDM16;ARHGEF16 | 159260;-225647 |

![image.png](https://cdn.nlark.com/yuque/0/2023/png/32598292/1681199819026-3ebb7143-a50c-4a98-a94f-8267b2c46e31.png#averageHue=%23f8f7f7&clientId=u4ba62734-d942-4&from=paste&height=704&id=uffc8ebf2&name=image.png&originHeight=704&originWidth=1689&originalType=binary&ratio=1&rotation=0&showTitle=false&size=108790&status=done&style=none&taskId=u4f632fea-2961-4a50-93c6-13e19d2cd84&title=&width=1689)
<a name="MJIgK"></a>
# Pictures
There are currently five plots in CscMeth, with 2 incorporated in fucntions and the other three customized ones:<br />`cscmeth elements`, showed [above](#FWNJq).<br />`cscmeth tanghulu`, suitable for visualizing methylation data in short range.<br />`cscmeth minihg`, suitable for visualizing methylation data in long range.
```
### Please note here what follows --allc_path is a directory where ALLC files to be ploted
### are. And This step may take a while if there are many ALLC files.

cscmeth tanghulu --allc_path ../../allcools/0.plot/ \
--range chr6:566781-570000\
```
![image.png](https://cdn.nlark.com/yuque/0/2023/png/32598292/1681200229379-1531c516-7833-4704-9477-793ad8dcbed0.png#averageHue=%23ededed&clientId=u4ba62734-d942-4&from=paste&height=910&id=u1c63508a&name=image.png&originHeight=910&originWidth=1316&originalType=binary&ratio=1&rotation=0&showTitle=false&size=101108&status=done&style=none&taskId=u226cfa4f-0a19-404a-b745-104cabae874&title=&width=1316)
```
### Please note here what follows --allc_path is a directory where ALLC files to be ploted
### are. And This step may take a while if there are many ALLC files.

cscmeth minihg --allc_path ../../allcools/0.plot/ \
--range chr6:566781-570000\
```
![image.png](https://cdn.nlark.com/yuque/0/2023/png/32598292/1681200267531-d3f7eb83-bcd8-4440-9e9e-64d21eb6eb1e.png#averageHue=%23fafafa&clientId=u4ba62734-d942-4&from=paste&height=1080&id=uc9fe9644&name=image.png&originHeight=1080&originWidth=527&originalType=binary&ratio=1&rotation=0&showTitle=false&size=70315&status=done&style=none&taskId=u77841434-97d7-457f-af2c-b1c7e7a76a2&title=&width=527)<br />And 
