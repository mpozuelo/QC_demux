#!/usr/bin/env python


import pandas as pd
import argparse
import sys

def parse_args(args=None):
    Description = 'Check samplesheet and add files to bed file'
    Epilog = """Example usage: python modify_samplesheet.py <FILE_IN> <FILE_OUT>"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('FILE_IN', help="Input samplesheet.")
    parser.add_argument('FILE_OUT', help="Output samplesheet.")
    return parser.parse_args(args)


def add_bed_file(FileIn,FileOut):
    #Open input file
    fi = open(FileIn, 'r')

    # Load mosdepth thresholds.bed.gz into a pandas dataframe
    cov = pd.read_csv(fi, delimiter=',', index_col=False, low_memory=False)

    # Open output file
    fo = open(FileOut, 'w')


    basefolder = '/datos/ngs/dato-activo/References/iGenomes/'
    # Dictionary for bed files
    star_gtf = {'hg38': basefolder + 'Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf', 'mm38': basefolder + 'Mus_musculus/Ensembl/GRCm38/Annotation/Genes/genes.gtf'}
    star_index = {'hg38': basefolder + 'Homo_sapiens/NCBI/GRCh38/Annotation/Sequence/STARIndex/', 'mm38': basefolder + 'Mus_musculus/Ensembl/GRCm38/Annotation/Sequence/STARIndex/'}

    # Write header
    #fo.write("%s\n" %('\t'.join(l_th[1:])))

    # Compute percentages
    cov['star_gtf'] = cov['genome'].map(star_gtf)
    cov['star_index'] = cov['genome'].map(star_index)
    cov['fastq1'] = "/datos/ngs/dato-activo/data/04_pfastq/" + cov['platform'] + '/' + cov['machine'] + '/' + cov['run'] + '/' + cov['lane'] + '/' + cov['user'] + '/demux_fastq/' + cov['sampleID'].astype(str) + '_S' + cov['ID'].astype(str) + '_R1_001.fastq.gz'
    cov['fastq2'] = "/datos/ngs/dato-activo/data/04_pfastq/" + cov['platform'] + '/' + cov['machine'] + '/' + cov['run'] + '/' + cov['lane'] + '/' + cov['user'] + '/demux_fastq/' + cov['sampleID'].astype(str) + '_S' + cov['ID'].astype(str) + '_R2_001.fastq.gz'


    cov.to_csv(fo, index = False)
    fi.close()


def main(args=None):
    args = parse_args(args)
    add_bed_file(args.FILE_IN,args.FILE_OUT)


if __name__ == '__main__':
    sys.exit(main())
