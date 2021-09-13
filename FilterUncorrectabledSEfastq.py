"""
author: Jaydee Sereewit
ssereewit at gmail.com

**Modified from:
https://github.com/harvardinformatics/TranscriptomeAssemblyTools/blob/master/FilterUncorrectabledPEfastq.py
to handle single-end fastq files corrected using Rcorrector

This script takes as an input Rcorrector error corrected Illumina single-end reads
in fastq format and:

1. Removes any reads that Rcorrector indentifes as containing an error,
but can't be corrected, typically low complexity sequences. For these,
the header contains 'unfixable'.

2. Strips the ' cor' from headers of reads that Rcorrector fixed, to avoid
issues created by certain header formats for downstream tools.

3. Write a log with counts of (a) total reads processed, (b) total 'unfixable'
reads removed, (c) total reads retained, (d) number of reads corrected by 
RCorrector.

Currently, this script only handles single-end data, and handle either unzipped
or gzipped files on the fly, so long as the gzipped files end with 'gz'.

"""

import sys        
import gzip
from itertools import zip_longest
import argparse
from os.path import basename

def get_input_streams(sefile):
    if sefile[-2:]=='gz':
        sehandle=gzip.open(sefile,'rb')
    else:
        sehandle=open(sefile,'r')
    return sehandle

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)


if __name__=="__main__": 
    parser = argparse.ArgumentParser(description="options for filtering and logging rCorrector single-end fastq outputs")
    parser.add_argument('-i','--single_end_reads',dest='sereads',type=str,help='single end fastq file')
    parser.add_argument('-s','--sample_id',dest='id',type=str,help='sample name to write to log file')
    opts = parser.parse_args()

    seout=open('unfixrm_%s' % basename(opts.sereads).replace('.gz',''),'w')

    se_cor_count=0
    unfix_se_count=0

    se_stream=get_input_streams(opts.sereads)

    with se_stream as f1:
        SE=grouper(f1,4)
        counter=0
        for entry in SE:
            counter+=1
            if counter%100000==0:
                print("%s reads processed" % counter)

            head1,seq1,placeholder1,qual1=[i.strip() for i in entry]

            if 'unfixable' in head1:
                unfix_se_count+=1
            else:
                if 'cor' in head1:
                    se_cor_count+=1

                head1=head1.split('l:')[0][:-1] 
                seout.write('%s\n' % '\n'.join([head1,seq1,placeholder1,qual1]))

    total_retained = counter - unfix_se_count

    unfix_log=open('rmunfixable_%s.log' % opts.id,'w')
    unfix_log.write('total SE reads:%s\nremoved SE unfixable reads:%s\nretained SE reads:%s\nSE reads corrected:%s\n' % (counter,unfix_se_count,total_retained,se_cor_count))

    seout.close()
    unfix_log.close()
