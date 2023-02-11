# De-Novo transcriptome assembly

During this lesson we will reconstruct transcriptome assembly for *Heilongjiang brown frog (Rana amurensis)* using RNA-seq data obtained from [https://www.nature.com/articles/s41598-022-24631-6#data-availability](https://www.nature.com/articles/s41598-022-24631-6#data-availability). We already downloaded raw .fastq files and prepared subsets of reads for you. They are available here:
`/home/d.smirnov/homework5/data/`

In the `data` you will find your individual paired-end fastq files, `<your_number>.R1.fq.gz` and `<your_number>.R2.fq.gz`.

Please do not copy files to your directory and do not use more than two computational cores for this tasks!


## Set up conda environment
We will activate `trinity` environment with all required tools preinstalled.

``` bash
/opt/anaconda3/bin/conda init 
conda activate trinity
```

## Read filtering
Befor we start with assembly, it worth to remove adapters and low-quality reads from the raw .fastq files. To do this, we will run `TrimGalore` tool: 

``` bash
trim_galore --cores 2 --paired --gzip --fastqc --fastqc_args "-t 2" --output_dir data_filtered data/1.R1.fq.gz data/1.R2.fq.gz
```
`--fastqc` argument means that we also want to perform quality assesement of filtered reads via `FastQC` using two computational cores (as specified by `--fastqc_args "-t 2"`).


## Run Trinity
This is the most crucial and time consuming step in the analysis.

``` bash
Trinity --seqType fq --SS_lib_type RF --max_memory 2G --left data_filtered/<your_number>.R1_val_1.fq.gz --right data_filtered/<your_number>.R2_val_2.fq.gz --CPU 2
```
where `--SS_lib_type RF` means that the library preparation was done with a first strand protocol (we know that from the original paper)

❓ How many trinity 'genes' and 'transcripts' did you get? (1 point)


## Obtain assembly statistics
Run 
``` bash
TrinityStats.pl  Trinity.fasta
```

## Transcript filtering

``` bash
cd-hit-est -o cdhit -c 0.98 -i trinity_out_dir/Trinity.fasta -p 1 -d 0 -b 3 -T 2 -M 1000
```

After filtering new filtered assembly will be stored in 'cdhit' file. Rename it to avoid confusion: 
``` bash
mv cdhit Trinity.filtered.fasta
```

Run QC script on filtered assembly again:
``` bash
TrinityStats.pl  Trinity.filtered.fasta
```

❓ How many trinity 'transcripts' were filtered out from the assembly? (1 point)

## Assembly completeness with gVolante
To assess completeness of the obtained assembly we will use `gVolante` web tool. Transfer filtered .fasta transcriptome to your local machine and then open https://gvolante.riken.jp/analysis.html. 

Upload your .fasta file, it may take few minutes:

![**Figure 1**. Uploading an assembly to gVolante](gVolante1.png)

Select `Coding/transcribed (nucleotide)` sequence type and `BUSCO v5` as a pipeline to use. Select "Tetrapoda" as an ortholog set for `BUSCO v5`.

![**Figure 2**. Setting up ortholog database to search in](gVolante2.png)


## Transcript quantification
``` bash
align_and_estimate_abundance.pl --SS_lib_type RF --est_method kallisto --transcripts Trinity.filtered.fasta --seqType fq --left data_filtered/<your_number>.R1_val_1.fq.gz --right data_filtered/<your_number>.R2_val_2.fq.gz --output_dir kallisto_output --thread_count 2 --trinity_mode --prep_reference
```
Apart from *kallisto_output* folder it will also generate two additional files:
*Trinity.filtered.fasta.gene_trans_map*  and *Trinity.filtered.fasta.kallisto_idx*


## Transcript annotation
``` bash
Build_Trinotate_Boilerplate_SQLite_db.pl  Trinotate
```

``` bash
makeblastdb -in uniprot_sprot.pep -dbtype prot
```

``` bash
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
```

Run Transdecoder
``` bash
TransDecoder.LongOrfs -m 10 -t Trinity.filtered.fasta
```

``` bash
blastx -query Trinity.filtered.fasta -db uniprot_sprot.pep -num_threads 2 -max_target_seqs 1 -outfmt 6 > blastx.outfmt6
```

``` bash
blastp -query Trinity.filtered.fasta.transdecoder_dir/longest_orfs.pep -db uniprot_sprot.pep -num_threads 2 -max_target_seqs 1 -outfmt 6 > blastp.outfmt6
```


``` bash
hmmscan --cpu 2 --domtblout TrinotatePFAM.out Pfam-A.hmm Trinity.filtered.fasta.transdecoder_dir/longest_orfs.pep > pfam.log
```

Once all three searches are done, we can start uploading information to `Trinotate.sqlite` database:
``` bash
Trinotate Trinotate.sqlite init --gene_trans_map Trinity.filtered.fasta.gene_trans_map \
                                --transcript_fasta Trinity.filtered.fasta \
                                --transdecoder_pep Trinity.filtered.fasta.transdecoder_dir/longest_orfs.pep
```

Load blastp, blastx and pfam results
``` bash
Trinotate Trinotate.sqlite LOAD_swissprot_blastp blastp.outfmt6
Trinotate Trinotate.sqlite LOAD_swissprot_blastx blastx.outfmt6
Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
```

Generate the report:
``` bash
Trinotate Trinotate.sqlite report > trinotate_annotation_report.xls
```

``` bash
Trinotate Trinotate.sqlite report > trinotate_annotation_report.xls
```


``` bash
extract_GO_assignments_from_Trinotate_xls.pl  \
        --Trinotate_xls trinotate_annotation_report.xls \
        -G --include_ancestral_terms > go_annotations.txt
```

## Integrating expression data and annotation
Generate a map of trinity gene IDs and corresponding annotations

``` bash
Trinotate_get_feature_name_encoding_attributes.pl \
                  trinotate_annotation_report.xls  > annot_feature_map.txt
```
❓ How many genes were annotated ? (1 point)


Update expession matrix with gene names:
``` bash
rename_matrix_feature_identifiers.pl Trinity_trans.counts.matrix annot_feature_map.txt > Trinity_trans.counts.wAnnot.matrix
``


## References
[1] https://github.com/trinityrnaseq/trinityrnaseq/wiki

[2] https://github.com/Trinotate/Trinotate/wiki/Software-installation-and-data-required

## Assignment and grading

When you complete all the tasks please upload to Canvas:

* Report with plots, pictures, answers and explanations.

* Your code (txt/pdf/R markdown or something else readable)

How it will be graded:
* The maximum is 10 points.

* For missing step - minus 5 points.

* If the report/code is missing - 0 points.

* For each missing figure/plot - minus 2 points.

* For each missing answer to ❓ - minus 2 points.

* For any other mistake - minus 0.2 point.

