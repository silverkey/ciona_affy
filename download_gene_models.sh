#!/bin/bash

# --------------- #
#   GENE MODELS   #
# --------------- #

# JGI 1.0
wget ftp://ftp.jgi-psf.org/pub/JGI_data/Ciona/v1.0/ciona.mrna.fasta.gz
gunzip -c ciona.mrna.fasta.gz > jgi.fa

# KYOTOGRAIL2005
wget -U firefox --user=guest --password=welcome http://hoya.zool.kyoto-u.ac.jp/datas/grail2005.zip
unzip -p grail2005.zip > kyotograil.fa

# KH
wget -U firefox --user=guest --password=welcome http://hoya.zool.kyoto-u.ac.jp/datas/KHGene.zip
unzip -p KHGene.zip > kh.fa

# ENSEMBL
wget ftp://ftp.ensembl.org/pub/release-56/fasta/ciona_intestinalis/cdna/Ciona_intestinalis.JGI2.56.cdna.all.fa.gz
gunzip -c Ciona_intestinalis.JGI2.56.cdna.all.fa.gz > ensembl.fa

# --------------- #
#   GENE ANNOTA   #
# --------------- #

# JGI 1.0
wget http://www.aniseed.cnrs.fr/download/JGI_annot.txt.zip
unzip -p JGI_annot.txt.zip > jgi.ann

# KYOTOGRAIL2005
wget http://www.aniseed.cnrs.fr/download/KYOTOGRAIL_annot.txt.zip
unzip -p KYOTOGRAIL_annot.txt.zip > kyotograil.ann

# KH
wget http://www.aniseed.cnrs.fr/download/KH_annot.txt.zip
unzip -p KH_annot.txt.zip > kh.ann

# ENSEMBL
wget http://www.aniseed.cnrs.fr/download/ENSEMBL_annot.txt.zip
unzip -p ENSEMBL_annot.txt.zip > ensembl.ann

# ----------- #
# CLEANING... #
# ----------- #

rm *.gz
rm *.zip
