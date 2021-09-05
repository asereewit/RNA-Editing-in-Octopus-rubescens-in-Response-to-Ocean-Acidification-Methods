import csv

swissprot_hits_seqid=[]
swissprot_hits_entries=[]
with open('rubescens_transcriptome_ORF_swissprot_blastx_1bestalignment.csv', newline='') as csvfile:
  csvreader=csv.reader(csvfile)
  for row in csvreader:
    swissprot_hits_seqid.append(row[0])
    swissprot_hits_entries.append(row[1])

hits_seqid_set = set(swissprot_hits_seqid)
hits_accessionid_set = set(swissprot_hits_entries)

with open('swissprot_ORF_seqid.lst', 'w', newline='') as seqidfile:
    csvwriter = csv.writer(seqidfile)
    for row in hits_seqid_set:
        csvwriter.writerow([row])

with open('unique_swissprot_entries.csv', 'w', newline='') as entryfile:
    csvwriter = csv.writer(entryfile)
    for row in hits_accessionid_set:
        csvwriter.writerow([row])
