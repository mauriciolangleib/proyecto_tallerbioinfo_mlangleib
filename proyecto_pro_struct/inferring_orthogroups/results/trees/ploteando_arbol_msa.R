# cargo librerias
library(ape)
library(tidyverse)
library(magrittr)
library(ggtree)
library(Biostrings)
library(ggmsa)

# cargo arbol
arbol_casein_kinasas = ape::read.tree('OG0000070.modif.msa.correct_labels.contree')

# cargo msa
protein_sequences = '../aln/OG0000070.modif.correct_labels.msa'
x = readAAStringSet(protein_sequences)
data = ggmsa::tidy_msa(x)

# ploteo
p = ggtree(arbol_casein_kinasas) + geom_tiplab()
ggtree(arbol_casein_kinasas, branch.length = 'none') + 
 geom_text2(aes(subset = !isTip, label=label), nudge_x = 0.60) +
 geom_tiplab(offset = 0.5) -> p

# guardo plot
plot = facet_plot(p, geom = geom_msa, data = data,  panel = 'MSA',
               font = NULL, color = "Chemistry_AA") +
       xlim_tree(1)
