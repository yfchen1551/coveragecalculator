#!/home/usr/bin/python

import re
import os
import sys
import glob

from BCBio import GFF
from Bio import SeqIO

if len(sys.argv) != 6:
    print ('Usage: python3 CDS_rRNA_extract.py inpath outpath samtools bedtools data_type')
    sys.exit()
    
inpath, outpath, samtools, bedtools, data_type = sys.argv[1:]

def getRef(gbff, ref):
    from Bio import SeqIO
    input_handle = open(gbff, "r")
    output_handle = open(ref, "w")
    for seq_record in SeqIO.parse(input_handle, "genbank"):
        print ("Dealing with GenBank record %s" % seq_record.id)
        output_handle.write(">%s %s\n%s\n" % (
            seq_record.id,
            seq_record.description, 
            seq_record.seq))
    output_handle.close()
    input_handle.close()

def getIdx(samtools, ref):
    os.system(f'{samtools} faidx {ref}')
    
def getGffInfo(gff, outbed, data_type):
    gff_file = open(gff)
    out_bed = open(outbed, 'w')
    for line in gff_file:
        if re.search('^#', line):
            continue
        line_info = line.strip().split('\t')
        if line_info[2] == data_type:
            #print (line_info[8])
            locus_tag = re.search('locus_tag=(.*?);', line_info[8]).group(1)
            try:
                rna_type = re.search('product=(.*?)[;$]', line_info[8]).group(1)
            except:
                rna_type = re.search('product=(.*)$', line_info[8]).group(1)
            start, end = line_info[3:5]
            start = str(int(start)-1)
            out_bed.write(f'{line_info[0]}\t{start}\t{end}\t{locus_tag} {rna_type}\n')
    gff_file.close()
    out_bed.close()
    
def getfasta(bedtools, bed, ref, out):
    os.system(f'{bedtools} getfasta -fi {ref} -fo {out} -bed {bed} -fullHeader -nameOnly')
    
    
for file in glob.glob(f'{inpath}/*gbff'):
    file_body = '.'.join(file.split('/')[-1].split('.')[:-1])
    in_file = file
    out_file = f'{outpath}/{file_body}.gff'
    in_handle = open(in_file)
    out_handle = open(out_file, "w")
    GFF.write(SeqIO.parse(in_handle, "genbank"), out_handle)
    in_handle.close()
    out_handle.close()
    getRef(file, f'{outpath}/{file_body}.ref.fasta')
    getGffInfo(out_file, f'{outpath}/{file_body}.bed', data_type)
    getIdx(samtools, f'{outpath}/{file_body}.ref.fasta')
    getfasta(bedtools, f'{outpath}/{file_body}.bed', f'{outpath}/{file_body}.ref.fasta', f'{outpath}/{file_body}.extract.{sys.argv[4]}.fasta')
