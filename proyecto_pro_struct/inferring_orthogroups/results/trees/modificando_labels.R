# modificando labels
# cargo librerias
library(ape)
library(magrittr)
library(tidytidbits)
library(tidyverse)
library(bioseq)

# obtengo las labels completas y truncadas para cada secuencia de bicho
correlation_table.tibble = list.files('../../proteomes', full.names = T) %>%
	as.list() %>%
	purrr::map_dfr(., ~{bioseq::read_fasta(.x) %>% as_tibble()}) %>%
	dplyr::select(label) %>%
	dplyr::mutate(tree_lab = label %>% str_split(' ') %>% purrr::map_chr(1))

# cargo arboles
arbol_contree = ape::read.tree('OG0000070.modif.msa.contree')

# filtro a la tabla y armo diccionario
correlation_table.tibble %<>% dplyr::filter(tree_lab %in% arbol_contree$tip.label)

# armo diccionario
seqs_complete.dict = correlation_table.tibble$label
names(seqs_complete.dict) = correlation_table.tibble$tree_lab

# reemplazo labels
arbol_contree$tip.label %<>%
	as.list() %>%
	purrr::map(., ~{.x %>% tidytidbits::lookup_chr(., dict = seqs_complete.dict, default = identity)}) %>%
	unlist()

# exporto arbol
ape::write.tree(phy = arbol_contree, file = 'OG0000070.modif.msa.correct_labels.contree')

# modifico tambien los nombres en el alineamiento
bioseq::read_fasta('../aln/OG0000070.modif.msa', 'AA') %>% 
	as_tibble() %>%
	dplyr::mutate(label = label %>% tidytidbits::lookup_chr(., dict = seqs_complete.dict, default = identity) %>%
					str_replace_all(., ' ', '_')) %>%
	deframe() %>%
	bioseq::write_fasta(., '../aln/OG0000070.modif.correct_labels.msa')

############
# plotting #
############

# cargo librerias
library(ape)
library(tidyverse)
library(magrittr)
library(ggtree)
library(Biostrings)
library(ggmsa)

# cargo arbol
arbol_casein_kinasas = arbol_contree

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

ggtree::ggsave(filename = 'arbol_casein_kinasas_msa.pdf',
               device = "pdf", 
               width = 12, 
               height = 9 , 
               units = "in" , 
               limitsize = FALSE)
