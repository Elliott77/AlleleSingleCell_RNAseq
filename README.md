# AlleleSingleCell_RNAseq
Count maternal and paternal scRNA-seq reads based on parental variants
Requires a VCF file indicating the variants that distingish the maternal and paternal DNA.
Requires scRNA-seq alignment file BAM
Requres a Barcode file (TSV) with barcodes for individual cells.

In our case we mated C57/B6 mice with Castaneous mice. We sequcen tissue from the F1 hybrid mice.
These strains are seperated by about 500,000 years, enough to have many distingushing SNPs.
https://www.cell.com/neuron/pdf/S0896-6273(17)30057-0.pdf
