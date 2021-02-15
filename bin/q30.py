#!/usr/bin/env python

import os,sys
import fastq
import time

def qual_stat(qstr):
    q30 = 0
    for q in qstr:
        qual = q - 33
        if qual >= 30:
            q30 += 1
    return q30

def stat(filename):
    reader = fastq.Reader(filename)
    total_count = 0
    q20_count = 0
    q30_count = 0
    while True:
        read = reader.nextRead()
        if read == None:
            break
        total_count += len(read[3])
        q30 = qual_stat(read[3])
        q30_count += q30

    print("%.2f%%" % (round(100 * float(q30_count)/float(total_count), 2)))

def main():
    if len(sys.argv) < 2:
        print("usage: python q30.py <fastq_file>")
        sys.exit(1)
    stat(sys.argv[1])

if __name__ == "__main__":
    main()
