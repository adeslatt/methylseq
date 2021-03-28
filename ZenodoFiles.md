# Using Zenodo

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4625710.svg)](10.5281/zenodo.4625710)

To make the data more accessible and FAIR, the indexed files were transfered to Zenodo using [`zenodo-upload`](https://github.com/jhpoelen/zenodo-upload) from the `Georgetown Lombardi's Comprehensive Cancer Center's Wellstein-Lab` Amazon `S3` bucket.   

The data provided on Zenodo is an indexed version of the human genome hg19. The hg19 fasta file was originally downloaded from UCSC (http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz) and combined with a target portion of the phage lambda genome from 4500 to 6500 base pairs, originally downloaded from NCBI (https://www.ncbi.nlm.nih.gov/nuccore/NC_001416.1.fa.gz)The Bismark program to map bisulfite treated sequencing reads was used to generate the index using the bismark_genome_preparation command. The indexed files then compressed using tar and then pigz to make a small and. Ordered set of the Bisulfite Genome files.References:
Felix Krueger, Simon R. Andrews, Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications, Bioinformatics, Volume 27, Issue 11, 1 June 2011, Pages 1571–1572, https://doi.org/10.1093/bioinformatics/btr167


## Preparing to use [`zenodo-upload`](https://github.com/jhpoelen/zenodo-upload)

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

ii. create and activate a new conda environment `zenodo`.
```bash
conda create -n zenodo
conda activate zenodo
```

Now with the environment ready - begin the pre-requisites

## Clone the `[zenodo-upload](https://github.com/jhpoelen/zenodo-upload)` repository

```bash
git clone https://github.com/jhpoelen/zenodo-upload.git
```

## Satisfy the Pre-requisites for zenodo upload


i. Install [`jq`](https://stedolan.github.io/jq/).  

```bash
conda install -c conda-forge jq
```

ii. Installed `curl`

```bash
conda install -c conda-forge curl
```

iii. Bash

already satisfied.


## Changing MacOS default shell

From the terminal. Opening a new terminal, you may see:

```bash
Last login: Fri Mar 26 17:31:14 on ttys003

The default interactive shell is now zsh.
To update your account to use zsh, please run `chsh -s /bin/zsh`.
For more details, please visit https://support.apple.com/kb/HT208050.
```

Change now to the bash shell with the `chsh` (change shell) command

```bash
chsh -s /bin/bash
```

## Transfering data from an S3 bucket

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
aws s3 cp s3://file-to-be-copied . --recursive

## Upload the files from the MacBook Pro

Two files are uploaded for the publication and ease of access for others

i. hg19_lambda.tar.gz
ii. Bisulfite_Genome.tar.gz

Uploaded after reserving the `DOI` from Zenodo and getting a personal `zenodo token`, following the instructions [zenodo-upload](https://github.com/jhpoelen/zenodo-upload) provides I typed the following two commands:

```bash
 ./zenodo_upload.sh 4625710 ../methylseq/data/hg19_lambda/hg19_lambda.tar.gz
 ```
 
 and
 
 ```bash
  ./zenodo_upload.sh 4625710 ../methylseq/data/hg19_lambda/Bisulfite_Genome.tar.gz
  ```
  

