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
Befor we start with assembly, it worth to remove adapters and low-quality reads from the raw .fastq files. To do this, we will run `TrimGalore`: 

``` bash
trim_galore --cores 2 --paired --gzip --fastqc --fastqc_args "-t 2" --output_dir data_filtered data/1.R1.fq.gz data/1.R2.fq.gz
```
`--fastqc` argument means that we also want to perform quality assesement of filtered reads via `FastQC` using two computational cores (as specified by `--fastqc_args "-t 2"`).


## Run Trinity
This is the most crucial and time consuming step in the analysis.

``` bash
Trinity --seqType fq --max_memory 2G --left data_filtered/<your_number>.R1_val_1.fq.gz --right data_filtered/<your_number>.R2_val_2.fq.gz --CPU 2
```


## Obtain assembly statistics

``` bash
TrinityStats.pl  Trinity.fasta
```

## Transcript quantification


## Transcript annotation


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

