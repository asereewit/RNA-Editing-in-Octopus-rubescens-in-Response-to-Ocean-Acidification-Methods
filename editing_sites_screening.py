import pysam
import numpy as np
import csv
from  statsmodels.stats.multitest import multipletests
from multiprocessing import Process, Queue
from scipy import stats

dna_orf_alignment = pysam.AlignmentFile('octo_dna_orf_alignment_sorted.bam','rb')

octo_bam_file_name = ['octo1_rna_orf_alignment_sorted.bam',
                      'octo2_rna_orf_alignment_sorted.bam',
                      'octo3_rna_orf_alignment_sorted.bam',
                      'octo4_rna_orf_alignment_sorted.bam',
                      'octo5_rna_orf_alignment_sorted.bam',
                      'octo6_rna_orf_alignment_sorted.bam']
num_rna_alignment = len(octo_bam_file_name)

rna_orf_alignment_list = []
for i in octo_bam_file_name:
    rna_orf_alignment_list.append(pysam.AlignmentFile(i,'rb'))

orf_fasta = pysam.FastaFile('swissprotORF.fasta')

codon_dict = {'GCT': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
              'CGT': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
              'AGA': 'Arg', 'AGG': 'Arg', 'AAT': 'Asn', 'AAC': 'Asn',
              'GAT': 'Asp', 'GAC': 'Asp', 'TGT': 'Cys', 'TGC': 'Cys',
              'CAA': 'Gln', 'CAG': 'Gln', 'GAA': 'Glu', 'GAG': 'Glu',
              'GGT': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly',
              'CAT': 'His', 'CAC': 'His', 'ATG': 'Met', 'ATT': 'Ile',
              'ATC': 'Ile', 'ATA': 'Ile', 'CTT': 'Leu', 'CTC': 'Leu',
              'CTA': 'Leu', 'CTG': 'Leu', 'TTA': 'Leu', 'TTG': 'Leu',
              'AAA': 'Lys', 'AAG': 'Lys', 'TTT': 'Phe', 'TTC': 'Phe',
              'CCT': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
              'TCT': 'Ser', 'TCC': 'Ser', 'TCA': 'Ser', 'TCG': 'Ser',
              'AGT': 'Ser', 'AGC': 'Ser', 'ACT': 'Thr', 'ACC': 'Thr',
              'ACA': 'Thr', 'ACG': 'Thr', 'TGG': 'Trp', 'TAT': 'Tyr',
              'TAC': 'Tyr', 'GTT': 'Val', 'GTC': 'Val', 'GTA': 'Val',
              'GTG': 'Val', 'TAA': 'Stop', 'TGA': 'Stop', 'TAG': 'Stop'}


ref_seq = ''
with open('swissprotORF.fasta') as fafile: # no linebreaks for sequences
    data = fafile.read()
    ref_seq = data

ref_seq_index = []
with open('swissprotORF.fasta.fai') as faifile:
    index = csv.reader(faifile, delimiter = '\t')
    for row in index:
        ref_seq_index.append(row)

ref_seq_index_dict = {}
for i in range(len(ref_seq_index)):
    entry = ref_seq_index[i]
    ref_seq_index_dict[entry[0]] = {
        'length': int(entry[1]),
        'offset': int(entry[2]),
        'linebases': int(entry[3]),
        'linewidth': int(entry[4]),
        }

def search_refbase(seqid, locus):
    """Use fai index file to allow random access to reference base."""
    name = str(seqid)
    loc = int(locus)
    if loc >= 0 and loc <= ref_seq_index_dict[name]['length']-1:
        base_offset_loc = ref_seq_index_dict[name]['offset'] + loc
        return ref_seq[base_offset_loc]

"""
# Test search_refbase function with SeqIO from Biopython,
# which can be used as an alternative to find reference base.

ref_fa = []
for record in SeqIO.parse('swissprotORF.fasta','fasta'):
    ref_fa.append([record.id,record.seq])

for i in ref_fa:
    for j in range(len(i[1])):
        x = search_refbase(i[0],j)
        if x != i[1][j]:
            print('id: %s \n locus: %d \n refbase_search: %s \n BioSeqIO: %s \n' %
                (i[0], j, x, i[1][j]))
"""

def strong_editing_site_detection(dna_orf_alignment, ses, ubs):
    """This function finds the bases aligned to each locus
    in the reference. A strong editing site (stored in ses)
    is where the bases are uniform AND different from the
    reference. Loci where there are uniform bases are saved
    in ubs for further weak editing sites analysis.
    """
    mapped_ref_seq = []
    mapstats = dna_orf_alignment.get_index_statistics()
    for i in range(len(mapstats)):
        if mapstats[i][1] != 0: mapped_ref_seq.append(mapstats[i][0])

    counter = 0
    total_seq = len(mapped_ref_seq)

    for i in mapped_ref_seq:
        ref_seq_length = dna_orf_alignment.get_reference_length(i)
        coverage = dna_orf_alignment.count_coverage(i,quality_threshold=30)

        for x in range(ref_seq_length):
            dna_base_profile = ['',0,'',0,'']
            temp = np.array([0,0,0,0,0]) # temp = [#As, #Cs, #Gs, #Ts, total]

            for y in range(4):
                temp[y] = coverage[y][x]
            temp[4] = temp[0] + temp[1] + temp[2] + temp[3]

            if (temp[4] != 0 and (
                temp[0] == temp[4] or
                temp[1] == temp[4] or
                temp[2] == temp[4] or
                temp[3] == temp[4]
                )):
                dna_base_profile[0] = i
                dna_base_profile[1] = x
                dna_base_profile[2] = search_refbase(i,x)
                dna_base_profile[3] = temp[4]

                if temp[0] == temp[4]: dna_base_profile[4] = 'A'
                if temp[1] == temp[4]: dna_base_profile[4] = 'C'
                if temp[2] == temp[4]: dna_base_profile[4] = 'G'
                if temp[3] == temp[4]: dna_base_profile[4] = 'T'
                ubs.append(dna_base_profile)

                if dna_base_profile[2] != dna_base_profile[4]:
                    ses.append(dna_base_profile)

        counter += 1
        if counter % round(total_seq/1000) == 0:
            print("ses progress: %.1f%s" % ((counter/total_seq)*100, '%'))

def weak_editing_site_detection(rna_orf_alignment, q, sample_id):
    """screen for weak editing sites in ubs where the dna bases are uniform"""
    alignment = rna_orf_alignment
    wes = []
    counter = 0
    total = len(ubs)

    for i in ubs:
        coverage = alignment.count_coverage(contig=i[0],start=i[1],stop=i[1]+1,
                                            quality_threshold=30)
        temp = np.array([0,0,0,0,0])

        for j in range(4):
            temp[j] = coverage[j][0]
        temp[4] = temp[0] + temp[1] + temp[2] + temp[3]
        if temp[4] > 1:
            if i[2] == 'A': num_mismatch = temp[4] - temp[0]
            if i[2] == 'C': num_mismatch = temp[4] - temp[1]
            if i[2] == 'G': num_mismatch = temp[4] - temp[2]
            if i[2] == 'T': num_mismatch = temp[4] - temp[3]

            if num_mismatch > 0:
                p_binom = stats.binom_test(num_mismatch, temp[4], p=0.001,
                                           alternative='greater')
                if p_binom < 0.05:
                    wes.append([i[0], i[1], i[2],
                               temp[0], temp[1],
                               temp[2], temp[3],
                               temp[4], p_binom])

        counter += 1
        if counter % round(total/1000) == 0:
            print("wes sample %d progress: %.1f%s" % (sample_id, (counter/total)*100, '%'))

    wes_pval = [i[-1] for i in wes]
    multitest_results = multipletests(wes_pval,alpha=0.1,method='fdr_bh')

    for i, j in zip(multitest_results[1], range(len(wes))):
        wes[j][-1] = i

    for i in wes:
        if i[-1] < 0.05: q.put(i)

def count_bases(rna_orf_alignment, es, q, sample_id):
    """count the rna bases of editing sites in list es
    and put the rna bases profile in multiprocessing queue.
    """
    counter = 0
    total = len(es)
    for i in es:
        coverage = rna_orf_alignment.count_coverage(contig=i[0],start=i[1],
                                                    stop=i[1]+1,quality_threshold=30)
        temp = np.array([0,0,0,0,0])
        for y in range(4):
            temp[y] = coverage[y][0]
        temp[4]=temp[0]+temp[1]+temp[2]+temp[3]
        q.put(temp)

        counter += 1
        if counter % round(total/1000) == 0:
            print("sample %d count bases progress: %.1f%s" % (sample_id, (counter/total)*100, '%'))

def add_edit_info(es):
    """add genomic codon, edited codon, genomic amino acid, edited amino acid,
    edited position within codon, upstream base, and downstream base. 
    """
    for i in range(len(es)):
        orf_id = es[i][0]
        orf_seq = orf_fasta.fetch(orf_id)
        pos = int(es[i][1])
        transcriptome_base = es[i][2]
        dna_base = es[i][4]
        edit_pos = pos % 3
        codon = ''
        if pos != 0:
            upstream_base = orf_seq[pos-1]
        else: upstream_base = ''
        if pos < (len(orf_seq)-1):
            downstream_base = orf_seq[pos+1]
        else: downstream_base = ''

        if edit_pos == 0:
            codon = orf_seq[pos:pos+3]
        if edit_pos == 1:
            codon = orf_seq[pos-1:pos+2]
        elif edit_pos == 2:
            codon = orf_seq[pos-2:pos+1]
        assert(codon[edit_pos] == transcriptome_base)
        assert(codon[edit_pos] == orf_seq[pos])

        Atotal=es[i][5]+es[i][10]+es[i][15]+es[i][20]+es[i][25]+es[i][30]
        Ctotal=es[i][6]+es[i][11]+es[i][16]+es[i][21]+es[i][26]+es[i][31]
        Gtotal=es[i][7]+es[i][12]+es[i][17]+es[i][22]+es[i][27]+es[i][32]
        Ttotal=es[i][8]+es[i][13]+es[i][18]+es[i][23]+es[i][28]+es[i][33]
        edited_base=''
        if dna_base=='A':
            total=[Ctotal,Gtotal,Ttotal]
            max_index = total.index(max(total))
            if max_index == 0: edited_base='C'
            if max_index == 1: edited_base='G'
            if max_index == 2: edited_base='T'
        if dna_base=='C':
            total=[Atotal,Gtotal,Ttotal]
            max_index = total.index(max(total))
            if max_index == 0: edited_base='A'
            if max_index == 1: edited_base='G'
            if max_index == 2: edited_base='T'
        if dna_base=='G':
            total=[Atotal,Ctotal,Ttotal]
            max_index = total.index(max(total))
            if max_index == 0: edited_base='A'
            if max_index == 1: edited_base='C'
            if max_index == 2: edited_base='T'
        if dna_base=='T':
            total=[Atotal,Ctotal,Gtotal]
            max_index = total.index(max(total))
            if max_index == 0: edited_base='A'
            if max_index == 1: edited_base='C'
            if max_index == 2: edited_base='G'

        genomic_codon_list = list(codon)
        genomic_codon_list[edit_codon_pos] = dna_base
        genomic_codon="".join(genomic_codon_list)
        genomic_aminoacid = codon_dict[genomic_codon]

        edited_codon_list = list(codon)
        edited_codon_list[edit_codon_pos] = edited_base
        edited_codon = "".join(edited_codon_list)
        edited_aminoacid = codon_dict[edited_codon]

        for j in [genomic_codon,
                  edited_codon,
                  genomic_aminoacid,
                  edited_aminoacid,
                  edited_codon_pos,
                  upstream_base,
                  downstream_base]:
            es[i].append[j]		


if __name__ == '__main__':

    ses = []
    ubs = []
    strong_editing_site_detection(dna_orf_alignment, ses, ubs)

    wes_queue_list = []
    for i in range(num_rna_alignment):
        wes_queue_list.append(Queue())

    wes_process = []
    for x, y, z in zip(rna_orf_alignment_list,wes_queue_list,range(num_rna_alignment)):
        p = Process(target=weak_editing_site_detection,args=(x,y,z))
        wes_process.append(p)

    for i in wes_process:
        i.start()

    wes_all = []
    for i in wes_queue_list:
        temp = []
        while not i.empty():
            temp.append(i.get())
        wes_all.append(temp)

    for i in wes_process:
        i.join()

    wes_all_id_loc = []
    for i in wes_all:
        wes_id_loc = []
        for j in range(len(i)):
            temp = i[j][0]+','+str(i[j][1])
            wes_id_loc.append(temp)
        wes_all_id_loc.append(wes_id_loc)

    combined_wes = set()
    for i in wes_all_id_loc:
        combined_wes = combined_wes | set(i)
    combined_wes = list(combined_wes)

    ubs_dict={}
    for x in range(len(ubs)):
        entry = ubs[x][0]+','+str(ubs[x][1])
        ubs_dict[entry]=x

    wes = []
    for i in combined_wes:
        temp = ubs[ubs_dict[i]]
        wes.append(temp)

    wes_countbases_queue = []
    for i in range(num_rna_alignment):
        wes_countbases_queue.append(Queue())

    wes_countbases_process = []
    for x,y,z in zip(rna_orf_alignment_list,wes_countbases_queue,range(num_rna_alignment)):
        p = Process(target=count_bases,args=(x,wes[:],y,z))
        wes_countbases_process.append(p)

    for i in wes_countbases_process:
        i.start()

    wes_countbases_list = []
    for i in wes_countbases_queue:
        temp = []
        while not i.empty(): temp.append(i.get())
        wes_countbases_list.append(temp)

    for i in wes_countbases_process:
        i.join()

    #assert(len(wes)
    #       == len(wes_countbases_list[0])
    #       == len(wes_countbases_list[1])
    #       == len(wes_countbases_list[2])
    #       == len(wes_countbases_list[3])
    #       == len(wes_countbases_list[4])
    #       == len(wes_countbases_list[5])
    #       )

    for i in range(len(wes)):
        for j in range(len(wes_countbases_list)):
            wes[i] = wes[i] + wes_countbases_list[j][i]

    add_edit_info(wes)

    ses_countbases_queue = []
    for i in range(num_rna_alignment):
        ses_countbases_queue.append(Queue())

    ses_countbases_process = []
    for x,y,z in zip(rna_orf_alignment_list,ses_countbases_queue,range(num_rna_alignment)):
        p = Process(target=count_bases,args=(x,ses[:],y,z))
        ses_countbases_process.append(p)

    for i in ses_countbases_process:
        i.start()

    ses_countbases_list = []
    for i in ses_countbases_queue:
        temp = []
        while not i.empty(): temp.append(i.get())
        ses_countbases_list.append(temp)

    for i in ses_countbases_process:
        i.join()

    #assert(len(ses)
    #       == len(ses_countbases_list[0])
    #       == len(ses_countbases_list[1])
    #       == len(ses_countbases_list[2])
    #       == len(ses_countbases_list[3])
    #       == len(ses_countbases_list[4])
    #       == len(ses_countbases_list[5])
    #       )

    for i in range(len(ses)):
        for j in range(len(ses_countbases_list)):
            ses[i] = ses[i] + ses_countbases_list[j][i]

    add_edit_info(ses)

    wes_refid_loc = []
    for i in wes:
        temp = i[0] + ',' + str(i[1])
        wes_refid_loc.append(temp)

    ses_refid_loc = []
    for i in ses:
        temp = i[0] + ',' + str(i[1])
        ses_refid_loc.append(temp)

    wes_refid_loc = set(wes_refid_loc)
    ses_refid_loc = set(ses_refid_loc)
    shared_es = wes_refid_loc & ses_refid_loc

    aes = []
    for i in ses:
        temp = i[0] + ',' + str(i[1])
        if temp not in shared_es: aes.append(i)
    aes = aes + wes

    file = open('aes_profile.csv','a+',newline='')
    with file:
        write = csv.writer(file)
        write.writerows(aes)
    file = open('wes_profile.csv','a+',newline='')
    with file:
        write = csv.writer(file)
        write.writerows(wes)
    file = open('ses_profile.csv','a+',newline='')
    with file:
        write = csv.writer(file)
        write.writerows(ses)

