# ZenodoFiles.md

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4625710.svg)](10.5281/zenodo.4625710)

The data provided on Zenodo is an indexed version of the human genome hg19. The hg19 fasta file was originally downloaded from UCSC (http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz) and combined with a target portion of the phage lambda genome from 4500 to 6500 base pairs, originally downloaded from NCBI (https://www.ncbi.nlm.nih.gov/nuccore/NC_001416.1.fa.gz)The Bismark program to map bisulfite treated sequencing reads was used to generate the index using the bismark_genome_preparation command. The indexed files then compressed using tar and then pigz to make a small and. Ordered set of the Bisulfite Genome files.References:
Felix Krueger, Simon R. Andrews, Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications, Bioinformatics, Volume 27, Issue 11, 1 June 2011, Pages 1571–1572, https://doi.org/10.1093/bioinformatics/btr167


## Uploading files on Zenodo

Using [`zenodo-upload`](https://github.com/jhpoelen/zenodo-upload)

Followed required first steps:

This work was done on a `MacBookPro` running `macOS Big Sur`, `Version 11.2.3`.

Default shell on this version is `zsh`.  Actually there is not an option following the [instructions](https://support.apple.com/kb/HT208050). To change the shell for your account to `bash`.   So did it as follows:

## Changing MacOS default shell

From the terminal. Opening a new terminal, you may see:

```bash
Last login: Fri Mar 26 17:31:14 on ttys003

The default interactive shell is now zsh.
To update your account to use zsh, please run `chsh -s /bin/zsh`.
For more details, please visit https://support.apple.com/kb/HT208050.
```

i. Installed [`jq`](https://stedolan.github.io/jq/).  On the mac, used `brew install jq`

ii. Installed `curl`

iii. Bash

