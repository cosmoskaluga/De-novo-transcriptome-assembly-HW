# De-Novo transcriptome assembly

During this lesson we will reconstruct transcriptome assembly for Heilongjiang brown frog (Rana amurensis) using data obtained from [https://www.nature.com/articles/s41598-022-24631-6#data-availability](https://www.nature.com/articles/s41598-022-24631-6#data-availability). 


## Setting up conda environment
``` bash
/opt/anaconda3/bin/conda init 
conda activate trinity
```

## Read filtering
``` bash
trim_galore --cores 2 --paired --gzip --output_dir /gss/d.smirnov/Frogs/new_data/2022_03_23/Shekhovtsov_RSF_2021/filtered/ \
        ${sample}R1_001.fastq.gz \
        ${sample}R2_001.fastq.gz
```

## Run fastqc on the data 
``` bash
fastqc -t 2 /path/to/fastq/*.fastq.gz -o /path/to/output/folder/
```

## Run trinity

``` bash
Trinity --seqType fq --max_memory 2G --left reads_1.fq  --right reads_2.fq --CPU 2
```


## Obtain assembly statistics
``` bash
TrinityStats.pl  Trinity.fasta
```

## Transcript quantification


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

