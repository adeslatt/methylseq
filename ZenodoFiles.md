# Using Zenodo

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4625710.svg)](10.5281/zenodo.4625710)

To make the data more accessible and FAIR, the indexed files were transfered to Zenodo using [`zenodo-upload`](https://github.com/jhpoelen/zenodo-upload) from the `Georgetown Lombardi's Comprehensive Cancer Center's Wellstein-Lab` Amazon `S3` bucket.   

The data provided on Zenodo is an indexed version of the human genome hg19. The hg19 fasta file was originally downloaded from UCSC (http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz) and combined with a target portion of the phage lambda genome from 4500 to 6500 base pairs, originally downloaded from NCBI (https://www.ncbi.nlm.nih.gov/nuccore/NC_001416.1.fa.gz)The Bismark program to map bisulfite treated sequencing reads was used to generate the index using the bismark_genome_preparation command. The indexed files then compressed using tar and then pigz to make a small and. Ordered set of the Bisulfite Genome files.References:
Felix Krueger, Simon R. Andrews, Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications, Bioinformatics, Volume 27, Issue 11, 1 June 2011, Pages 1571–1572, https://doi.org/10.1093/bioinformatics/btr167


## Preparing to use [`zenodo-upload`](https://github.com/jhpoelen/zenodo-upload)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4625710.svg)](https://doi.org/10.5281/zenodo.4625710)

Followed required first steps:
## Preparing for `zenodo-upload` on Lifebit's CloudOS

<p align="center"><img src="https://github.com/lifebit-ai/dry-bench-skills-for-researchers/blob/adds-mini-courses/assets/lifebitCloudOS.png"  width="250" align="right" ></p>

Using Lifebit's CloudOS system, there are no steps to be followed to make access to a bucket that is already authorized in the workspace.  If it is not then a two step process needs to be followed.

## Uploading files on Zenodo from Lifebit's CloudOS JupyterLab Notebook


On `Lifebit's CloudOS`, use conda to organize your environment

i. initialize the bash environment

```bash
conda init bash
exec -l bash
```

ii. create and activate a new conda environment `methylseq`.
```bash
conda create -n zenodo
conda activate zenodo
```

Now with the environment ready - begin the pre-requisites

## Clone the [zenodo-upload](https://github.com/jhpoelen/zenodo-upload) repository

```bash
git clone https://github.com/jhpoelen/zenodo-upload.git
```

## Satisfy the Pre-requisites for zenodo upload


i. Install [`jq`](https://stedolan.github.io/jq/).  

```bash
conda install -c conda-forge jq -y
```

ii. Installed `curl`

```bash
conda install -c conda-forge curl -y
```

iii. Bash

already satisfied.


## Transfering data from an S3 bucket

In our case, we are taking files that were stored in an s3 bucket, tar'ing them up, `pigz`-ing them so that they are uploaded in the most efficient manner (and subsequently downloaded) for the purposes of analysis.

If required to copy from another bucket - you can install the appropriate client, the one needed here is `awscli`.

i. install command-line interface to Amazon Web Services `awscli`

```bash
(zenodo) $ conda install -c conda-forge awscli -y
```

ii. configure the aws client

```bash
(zenodo) $ aws configure
AWS Access Key ID [None]: [enter your aws access key]
AWS Secret Access Key [None]:  [enter your secret access key]
```

iii. cp the files

```bash
(zenodo) $ aws s3 cp s3://file-to-be-copied .
```

## tar and pigz

We want to make the files small and tight.   We will tar up loose files to ensure files are the smallest possible to ease in transfer - we will use [pigz](http://zlib.net/pigz/).   Searching [Anaconda](https://anaconda.org/conda-forge/pigz), we find a `conda install`.   Committed to managing our environments with `conda`, we use that.

```bash
(zenodo) $ conda install -c conda-forge pigz
```

Here are the files we wish to transfer.  Note they are not `tar'd` or `zipped`, so we will use `tar` and `pigz`.

## hg19_lambda.tar.gz

There are two different sets we will assemble.   We will make a `hg19_lambda.fa.tar.gz` from the following files.

```bash
(zenodo) $ aws s3 ls s3://file-to-be-copied

2020-11-01 22:32:28 3199908007 hg19_lambda.fa
2020-11-01 22:32:43       8591 hg19_lambda.fa.bis.amb
2020-11-01 22:32:44       4127 hg19_lambda.fa.bis.ann
2020-11-01 22:32:45  784290818 hg19_lambda.fa.bis.pac
2020-11-01 22:33:06 3137163372 hg19_lambda.fa.dau.bwt
2020-11-01 22:33:40 1568581688 hg19_lambda.fa.dau.sa
2020-11-01 22:35:07 3137163372 hg19_lambda.fa.par.bwt
2020-11-01 22:35:40 1568581688 hg19_lambda.fa.par.sa
```

First we `tar` 

```bash
(zenodo) $ tar cvf hg19_lambda.fa.tar hg19_lambda.fa*
```

and then we `pigz`.


```bash
pigz hg19_lambda.fa.tar
```

reducing `hg19_lambda.fa.tar` filesize from `7.1G` to `3.9G`.


## Bisulfite_Genome.tar.gz

We will make a `Bisulfite_Genome.tar.gz` from these files that `Bismark` will expect to find.

At the highest level we find

```bash
(zenodo) $ ls -l Bisulfite_Genome
total 8
drwxr-xr-x 2 jovyan users 4096 Apr  3 18:55 CT_conversion
drwxr-xr-x 2 jovyan users 4096 Apr  3 19:00 GA_conversion
```

Then in each of these two different directories we have:

```bash
(zenodo) $ ls -l Bisulfite_Genome/CT_conversion/
total 7141488
-rw-r--r-- 1 jovyan users  969974376 Nov  1 22:30 BS_CT.1.bt2
-rw-r--r-- 1 jovyan users  724328120 Nov  1 22:30 BS_CT.2.bt2
-rw-r--r-- 1 jovyan users       4823 Nov  1 22:30 BS_CT.3.bt2
-rw-r--r-- 1 jovyan users  724328116 Nov  1 22:30 BS_CT.4.bt2
-rw-r--r-- 1 jovyan users  969974376 Nov  1 22:30 BS_CT.rev.1.bt2
-rw-r--r-- 1 jovyan users  724328120 Nov  1 22:30 BS_CT.rev.2.bt2
-rw-r--r-- 1 jovyan users 3199909184 Nov  1 22:31 genome_mfa.CT_conversion.fa
```
and 

```bash
(zenodo) $ ls -l Bisulfite_Genome/GA_conversion/
total 7141488
-rw-r--r-- 1 jovyan users  969974376 Nov  1 22:31 BS_GA.1.bt2
-rw-r--r-- 1 jovyan users  724328120 Nov  1 22:31 BS_GA.2.bt2
-rw-r--r-- 1 jovyan users       4823 Nov  1 22:31 BS_GA.3.bt2
-rw-r--r-- 1 jovyan users  724328116 Nov  1 22:31 BS_GA.4.bt2
-rw-r--r-- 1 jovyan users  969974376 Nov  1 22:31 BS_GA.rev.1.bt2
-rw-r--r-- 1 jovyan users  724328120 Nov  1 22:32 BS_GA.rev.2.bt2
-rw-r--r-- 1 jovyan users 3199909184 Nov  1 22:32 genome_mfa.GA_conversion.fa
```
First we `tar` 

```bash
(zenodo) $ tar cvf Bisulfite_Genome.tar Bisulfite_Genome/*
```

and then we `pigz`.

```bash
pigz Bisulfite_Genome.tar
```

reducing `Bisulfite_Genome.tar` filesize from `13G` to `7.5G`.


## upload to Zenodo

Uploaded after reserving the `DOI` from Zenodo and getting a personal `zenodo token`, following the instructions [zenodo-upload](https://github.com/jhpoelen/zenodo-upload), I set the ZENODO_TOKEN environment variable.  

```bash
export ZENODO_TOKEN=[`set to your own personal zenodo token`]
```

## First create the repository on Zenodo

First, be sure you have created and saved but not published the repository on Zenodo.   Saving it reserves a DOI.  Do not publish until all desired data files are present.

i. upload the `hg19_lambda.tar.gz`

```bash
 ./zenodo_upload.sh 4625710 ../methylseq/data/hg19_lambda/hg19_lambda.tar.gz
 ```
 

ii.  and
 
 ```bash
  ./zenodo_upload.sh 4625710 ../methylseq/data/hg19_lambda/Bisulfite_Genome.tar.gz
  ```
then it will upload -- this is using cloudOS on Google and we have transfered the files from an aws s3 bucket to make it even more accessible through zenodo.

```bash
(methylseq) jovyan@0b595f0bbcd0:/mnt/shared/gcp-user/session_data/zenodo-upload$ ./zenodo_upload.sh 4625710 ../methylseq/data/hg19_lambda.fa.tar.gz 
++ echo 4625710
++ sed 's+^http[s]*://zenodo.org/deposit/++g'
+ DEPOSITION=4625710
+ FILEPATH=../methylseq/data/hg19_lambda.fa.tar.gz
++ echo ../methylseq/data/hg19_lambda.fa.tar.gz
++ sed 's+.*/++g'
+ FILENAME=hg19_lambda.fa.tar.gz
++ curl -H 'Accept: application/json' -H 'Authorization: Bearer xxxxxxxxxxxxxxxxxxx' https://www.zenodo.org/api/deposit/depositions/4625710
++ jq --raw-output .links.bucket
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100  3067  100  3067    0     0   2248      0  0:00:01  0:00:01 --:--:--  2248
+ BUCKET=https://zenodo.org/api/files/38402d50-855d-4815-9388-f2d231c2952a
+ curl --progress-bar -o /dev/null --upload-file ../methylseq/data/hg19_lambda.fa.tar.gz 'https://zenodo.org/api/files/38402d50-855d-4815-9388-f2d231c2952a/hg19_lambda.fa.tar.gz?access_token=PSpI0gExGdMyWLvaPJAlCAetPJxOkavuRWYReKxoZOvbZakC7S0AJXb0moGE'
######                                                                                                                                                4.5%

