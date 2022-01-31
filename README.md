# What is AlleleSingleCell_RNAseq?
This python script counts maternal and paternal scRNA-seq reads based on maternal and paternal variants.

# Python Libraries
Python 2.7
Python scrict requires the python packages:
pysam
pandas 

To install with anaconda:
conda install -c bioconda pysam
conda install pandas

# Run Python Script
python count_bam_sc_gene.py <BAM File> <VCF> <TSV Barcodes> <out.txt> <Number of Cores>
  
Requires a VCF file indicating the variants that distingish the maternal and paternal DNA.
Requires scRNA-seq alignment file BAM
Requres a barcode file (TSV) with barcodes for individual cells.
Number of Cores: for multiprocessing.   

In our case, we mated C57/B6 and Castaneous mice. We sequced cells from the F1 hybrid offspring mice.
These inbred strains are seperated by about 500,000 years of evolution, enough to have many distingushing variants but still interbeed. 
This allows us to distinguish many RNA-seq reads as maternally or paternally transcribed. See paper for details:
https://www.cell.com/neuron/pdf/S0896-6273(17)30057-0.pdf

# Output
Row names are gene IDs. Each cell has a column of gene level read counts for each allele.
  
