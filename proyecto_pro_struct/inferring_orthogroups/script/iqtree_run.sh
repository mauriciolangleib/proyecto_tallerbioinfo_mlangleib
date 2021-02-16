#!/bin/bash

# sin trimal,sobre modif
iqtree-omp -s ../results/aln/OG0000070.modif.msa -m TEST -ntmax 16 -alrt 1000 -bb 1000
