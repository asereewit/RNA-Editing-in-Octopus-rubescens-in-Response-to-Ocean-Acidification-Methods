# RNA-Editing-in-Octopus-rubescens-in-Response-to-Ocean-Acidification-Methods
Commands and Code for RNA Editing Detection in Octopus Rubescens in Response to Ocean Acidification


* Quality profiling and reads correction on RNA-seq data using Fastp (version 0.20.0)


  ```
  fastp -i Octo1.fastq -o Octo1_trimmed.fastq -h Octo1FastpReport
  ```
  
  
* Kmer-based error correction on RNA-seq data using rCorrector (version 1.0.4). Download rCorrector and run the run_rcorrector.pl inside the folder. 


  ```
  perl run_rcorrector.pl -s Octo1_trimmed.fastq 
  ```


* rRNA removal from RNA-seq data


  * Download O. cyanea 28s rRNA and O. vulgaris 18s rRNA from SILVA database and build a Bowtie2 (version 2.3.4.1) index 


    ```
    bowtie2-build -f Ovulgaris18srRNA.fasta,Ocyanea28srRNA.fasta Ovulgaris18sOcyanea28srRNA  
    ```


  * Map RNA-seq data to the combined Bowtie index. Only reads that are *not* aligned to the rRNA sequences will be used for further analysis


    ```
    bowtie2 --quiet --very-sensitive-local --phred33 -x Ovulgaris18sOcyanea28srRNA -U unfixrm_Octo1_trimmed.cor.fq --threads 20 --met-file blacklist_unpaired_unfixrm_Octo1_trimmed_bowtie2_metrics.txt --al blacklist_unpaired_aligned_unfixrm_Octo1_trimmed.fq --un blacklist_unpaired_unaligned_unfixrm_Octo1_trimmed.fq
    ```
  
  
* Trinity (version 2.8.5) transcriptome assembly


  ```
  Trinity --seqType fq --max_memory 240G --single blacklist_unpaired_unaligned_unfixrm_Octo1_trimmed.fq,blacklist_unpaired_unaligned_unfixrm_Octo2_trimmed.fq,blacklist_unpaired_unaligned_unfixrm_Octo3_trimmed.fq,blacklist_unpaired_unaligned_unfixrm_Octo4_trimmed.fq,blacklist_unpaired_unaligned_unfixrm_Octo5_trimmed.fq,blacklist_unpaired_unaligned_unfixrm_Octo6_trimmed.fq --SS_lib_type R --CPU 20 --output trinity_out_dir > run.log 2>&1
  ```
  
  
* Busco (version 4.1.4) evaluation on transcriptome 


  ```
  busco -i Trinity.fasta -l mollusca_odb10 -o rubescence_tran_busco -m tran -c 20
  ```
  
  
* Identify Open Reading Frames (ORFs) in transcriptome using ORFfinder (version 0.4.3) from NCBI


  ```
  ORFfinder -in Trinity.fasta -out rubescens_transcriptome_ORF_ignore_nested.fasta -n true -outfmt 1 > rubescens_transcriptome_ORFfinder_log.txt
  ```
  
  
* Find ORFs that are significantly matched (Blastx e-value < 1e-6) to SwissProt database 


  * Download the SwissProt database 


    ```
    perl /ncbi-blast/ncbi-blast-2.10.1+/bin/update_blastdb.pl --decompress swissprot
    ```

  
  * Blastx (version 2.10.1) transcriptome ORFs againsts SwissProt database


    ```
    /ncbi-blast/ncbi-blast-2.10.1+/bin/blastx -query rubescens_transcriptome_ORF_ignore_nested.fasta -out rubescens_transcriptome_ORF_swissprotdb_blastx_1bestalignemnt.txt -db /ncbi-blast/blastdb/swissprot -evalue 1e-6 -max_target_seqs 1 -subject_besthit -outfmt "6 qaccver saccver pident bitscore evalue" -num_threads 16
    ```
    
  * Select the transcript names that match to the Swissprot database using the Python snippet below:


    ```python
    import csv
    
    swissprot_hits_seqid=[]
    with open('rubescens_transcriptome_ORF_swissprot_blastx_1bestalignment.txt') as txtfile:
      csvreader=csv.reader(txtfile, delimiter='\t')
      for row in csvreader:
        swissprot_hits_seqid.append(row[0])
      
    hits_seqid_set = set(swissprot_hits_seqid) # remove duplicate transcript names in swissprot blastx hits
    
    with open('seqidname.lst','w') as idlist:
      csvwriter = csv.writer(idlist)
      csvwriter.writerows(hits_seqid_set)
    ```
    
    
  * Subset Swissprot ORFs using ```seqtk``` (version 1.2-r94)


    ```
    seqtk subseq rubescens_transcriptome_ORF_ignore_nested.fasta seqidname.lst > swissprotORF.fasta
    ```
    
    
    
