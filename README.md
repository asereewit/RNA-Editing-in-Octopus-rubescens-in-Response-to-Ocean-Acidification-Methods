# RNA-Editing-in-Octopus-rubescens-in-Response-to-Ocean-Acidification-Methods
Commands and Code for RNA Editing Detection in Octopus Rubescens in Response to Ocean Acidification


* Quality trimming and quality control on DNA-seq data using Trim Galore! (version 0.6.5)


  ```
  trim_galore --fastqc --paired M1_1.fq M1_2.fq M2_1.fq M2_2.fq M3_1.fq M3_2.fq M4_1.fq M4_2.fq M5_1.fq M5_2.fq
   ```


* Combine DNA reads from all samples


  ```
  cat M1_1_val_1.fq M2_1_val_1.fq M3_1_val_1.fq M4_1_val_1.fq M5_1_val_1.fq > pooled_trimmed_reads_1.fq
  ```
  ```
  cat M1_2_val_2.fq M2_2_val_2.fq M3_2_val_2.fq M4_2_val_2.fq M5_2_val_2.fq > pooled_trimmed_reads_2.fq
  ```


* Quality profiling and reads correction on RNA-seq data using Fastp (version 0.20.0) (Repeat for other samples)


  ```
  fastp -i Octo1.fastq -o Octo1_trimmed.fastq -h Octo1FastpReport
  ```
  
  
* Kmer-based error correction on RNA-seq data using rCorrector (version 1.0.4)


  ```
  perl run_rcorrector.pl -s Octo1_trimmed.fastq,Octo2_trimmed.fastq,Octo3_trimmed.fastq,Octo4_trimmed.fastq,Octo5_trimmed.fastq,Octo6_trimmed.fastq
  ```


* rRNA removal from RNA-seq data


  * Download O. cyanea 28s rRNA and O. vulgaris 18s rRNA from SILVA database and build a Bowtie2 (version 2.3.4.1) index 


    ```
    bowtie2-build -f Ovulgaris18srRNA.fasta,Ocyanea28srRNA.fasta Ovulgaris18sOcyanea28srRNA  
    ```


  * Map RNA-seq data to the combined Bowtie index (repeat for other RNA-seq samples) (Only reads that are *not* aligned to the rRNA sequences will be used for further analysis)


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
    
    
  * Subset Swissprot ORFs from transcriptome using ```seqtk``` (version 1.2-r94)


    ```
    seqtk subseq rubescens_transcriptome_ORF_ignore_nested.fasta seqidname.lst > swissprotORF.fasta
    ```
    
    
  * Build a Bowtie2 index of Swissprot ORFs


    ```
    bowtie2-build -f swissprotORF.fasta swissprotORF
    ```
    
    
* Map DNA reads to Swissprot ORFs


  * DNA reads Swissprot ORFs alignment (Samtools version: 1.7) (-F 260 retains only primary mapped reads)


    ```
    bowtie2 --local --threads 16 --quiet -t --met-file pooled_gDNA_orf_alignment_bowtie2_metrics.txt -q -x swissprotORF -1 pooled_trimmed_reads_1.fq -2 pooled_trimmed_reads_2.fq | samtools view -b -F 260 --threads 20 > primary_gDNA_orf_alignment.bam
    ```


  * Sort gDNA-ORF alignment file


    ```
    samtools sort -@ 16 -o primary_gDNA_orf_alignment_sorted.bam primary_gDNA_orf_alignment.bam
    ```
    
    
  * Index gDNA-ORF alignment file


    ```
    samtools index -b -@ 16 primary_gDNA_orf_alignmen_sorted.bam
    ```


* Map RNA reads to Swissprot ORFs (repeat for other samples)

  * RNA reads Swissprot ORFs alignment (Samtools version: 1.7) (-F 260 retains only primary mapped reads)


    ```
    bowtie2 --local --threads 16 --quiet -t --met-file octo1_orf_alignment_bowtie2_metrics.txt -q -x swissprotORF -U blacklist_unpaired_unaligned_unfixrm_Octo1_trimmed.fq | samtools view -b -F 260 --threads 20 > primary_octo1_orf_alignment.bam
    ```


  * Sort RNA-ORF alignment file


    ```
    samtools sort -@ 16 -o primary_octo1_orf_alignment_sorted.bam primary_octo1_orf_alignment.bam
    ```
    
    
  * Index RNA-ORF alignment file


    ```
    samtools index -b -@ 16 primary_octo1_orf_alignment_sorted.bam
    ```
    
  
