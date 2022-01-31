#!/usr/bin/env python
## count_bam_fetch_sc_gene.py
## Eliott Ferris
## 12/3/19

## Single Cell RNA-seq. Count Cast and C57 reads for each gene for each cell

useage = """
python count_bam_sc_gene.py <BAM File> <VCF> <TSV Barcodes> <out.txt> <Number of Cores>

"""

import pysam
import csv
import os
import pandas as pd
import re
import sys
import time
import math
import multiprocessing as mp
#from multiprocessing import Process

##methods
def remove_comma(snp):
    p = re.compile(r"\,")
    m = p.search(snp)
    if m:
        snp2 = snp[:m.start()]        
    else:
        snp2 = snp
    return snp2    
    
def getTrio(seq, index):
    return seq[index-1: index + 2]

def getCounts(gene00): ## , vcf_df00, samfile00, barcode_dict00
    vcf_df0000 = vcf_df[vcf_df.GeneIDs == gene00]
    vcf_df00 = vcf_df0000.drop_duplicates()
    samfile = pysam.Samfile(sys.argv[1], "rb")
    #print vcf_df00.shape

    counts_df = pd.DataFrame(columns= col_names_both, index = vcf_df00['ChrPos'], data = 0)
    counts_neither = pd.DataFrame(columns = col_names, index = vcf_df00['ChrPos'], data = 0)
    redandant_reads = 0
    already_written = []
    snp_count = 0
    ii = 0; total_read_count = 0; c57_count = 0; cast_count = 0; neither_count = 0; total_reads = 0
    last_index = ""; not_in_positions = 0; unknown_cell_barcode = 0
    print "g",
    for index00, row in vcf_df00.iterrows():
        if index00 == last_index:
            continue
        ii += 1;
        snp_count += 1
##        print "snp",
        chromosome = str(row['CHROM']) ## "CAST.Mar26_" +
        pos = int(row['POS'])## long?
        chr_pos = "%s_%d" %(chromosome, pos)
####        print chr_pos
        c57_snp = remove_comma(row['REF'])
        cast_snp = remove_comma(row['ALT'])

        c57_snp = (row['REF']).split(",")[0]
        cast_snp = (row['ALT']).split(",")[0]        
##        print 'expected C57:', c57_snp, "expected CAST:", cast_snp
##        print 'vcf_row',
        bam_path = os.path.dirname(sys.argv[1])
        os.chdir(bam_path)
        snp_reads = 0
        for read in samfile.fetch(chromosome, int(pos)-1, int(pos) ):  ## "CAST.Mar26_" +
            if read.is_duplicate == True:
                continue     

            if total_reads > 200000 or snp_reads > 50000:
                print 'total gene reads > 50000 or snp_reads > 10000:', str(total_reads)
                break
            ref_pos = read.get_reference_positions()
            
            if not int(pos) - 1  in ref_pos:
##                print "pos", str(pos - 1), "not found"
##                print " - - - ",
                not_in_positions += 1
                ## this usally means this spans the SNP but does not align to it.
                continue


            read_id = read.query_name + str(read.is_read1)[:1]## or read.qname
            if read_id in already_written:
                redandant_reads += 1
                continue

            cell_barcode = read.get_tag("CB") ##
            if not cell_barcode in barcode_dict00.keys():
                unknown_cell_barcode += 1

                continue
##            else:
####                    print "*** cell barcode found ***", cell_barcode,
            total_read_count += 1
            total_reads += 1
            snp_reads +=1
##            print "r", str(total_reads),
            cluster_id = barcode_dict00[cell_barcode]
    ##        print "cell_barcode:", cell_barcode, "cluster_id:", cluster_id, type(cluster_id)
            #pos_list = read.positions
            ref_pos = read.get_reference_positions()

            index2 = ref_pos.index((pos - 1))
##            print 'index from pos -1', str(index2)
            query_call = read.query_sequence[index2]
            qstring = read.query_alignment_sequence
##            print 'expected C57:', c57_snp, "expected CAST:", cast_snp
##            print "Trio", getTrio(qstring, index2)

            if query_call == c57_snp:
                counts_df.loc[chr_pos, ("c57_" + cluster_id)] += 1
##                print "c57",
                c57_count += 1
                already_written.append(read_id)

            elif query_call == cast_snp:
                counts_df.loc[chr_pos, ("cast_" + cluster_id)] += 1
##                print "cast",
##                print counts_df.loc[chr_pos,:].sum() 
                cast_count += 1
                already_written.append(read_id)

            else:    
                counts_neither.loc[chr_pos, cluster_id] += 1
                neither_count += 1
        last_index = index00         
    samfile.close()

    series_counts = counts_df.sum(axis = 0) ## .iloc[:, 1:]
    series_counts.rename(index = gene00)
    df_counts_out = series_counts.to_frame().T#.set_index(gene00)

    df_counts_out.insert(len(series_counts), 'GeneID', gene00)

    del(already_written)
    del(neither_count)
    del(counts_df)
    return df_counts_out #.rename(column = gene00)

############################################################################
start = time.time()
barcode_dict00 = {}

os.chdir("<barcodes directory>")

f_barcode = open(sys.argv[4], "a")

with open(sys.argv[3]) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter = ',')
    #skip header
    ##header = next(csv_reader)
    i = 0
    for row in csv_reader:
        i += 1
        barcode_dict00[(row[0])] = str(i)
        print >> f_out, str(i)+ "\tc57_" + str(i) + "\tcast_" + str(i) + "\t" + row[0]

print 'len(barcode_dict00)',
print len(barcode_dict00)

print 'len(barcode_dict00)',
print len(barcode_dict00)
col_names = barcode_dict00.values()
col_names_both = []
for col00 in col_names:
    col_names_both.append("c57_" + col00)
    col_names_both.append("cast_" + col00)
    
print 'len(col_names_both):',
print len(col_names_both)
    

#clusters = set(x for x in barcode_dict00.values())
print "CSV loaded"
os.chdir('<BAM directory>')

print "bam loaded"

# open the number of bam files and the same number of clusters, and map the out file handler to the cluster id, write to a bam with wb

## open VCF
dtypes00 = {"CHROM": str, "POS": int, "ID": str, "REF": str, "ALT": str, "QUAL": str, "FILTER": str, "INFO": str, "FORMAT": str, "CAST_EiJ": str}

snpFileName = sys.argv[2]
try:
    vcf_df = pd.read_csv(snpFileName, comment = "#", delimiter = '\t', dtype = dtypes00, names = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "CAST_EiJ"])
    
except:
    print "error opening the file", snpFileName

vcf_df['ChrPos'] = vcf_df.CHROM.astype(str).str.cat(vcf_df.POS.astype(str), sep = "_")

## chop_up by gene
gene_ids = []    
for info in vcf_df['INFO']:
    split_info = info.split("|")
    if len(split_info) > 2 and split_info[1][:3] == "ENS":
        gene_ids.append(split_info[1])
    else:
        gene_ids.append("NoGene")

vcf_df['GeneIDs'] = pd.Series(gene_ids, index = vcf_df.index)    
print vcf_df.head()
genes = list(set(gene_ids))

print "number of genes", str(len(genes))

## paralelize
pool = mp.Pool(sys.argv[5])
list_of_df = (pool.map(getCounts, [gene for gene in genes]))

pool.close()
pool.join()
print 'len(list_of_df):',
print len(list_of_df)
####################################################################################
df_out = pd.concat(list_of_df)
print 'df_out.iloc[:10, :6]',
print df_out.iloc[:10, :6]

print "concat list of df done, time:"
print time.time() - start

df_out.to_csv(sys.argv[4], sep = '\t', index = False)



