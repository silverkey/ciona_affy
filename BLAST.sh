#!/usr/bin/bash

makeblastdb -in aniseed_models.fa -dbtype nucl -title aniseed_models.blast -out aniseed_models.blast

blastn -query CINT06a520380F_sif.fa -task megablast -db aniseed_models.blast -out MEGABLAST_DEFAULT.out
