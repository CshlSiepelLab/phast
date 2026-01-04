# src/lib/CMakeLists.txt

# List your library sources here (mirror what your Makefile uses)
add_library(phast_core STATIC
base/phast_complex_matrix.c
base/phast_complex_vector.c
base/phast_eigen.c
base/phast_hashtable.c
base/phast_list_of_lists.c
base/phast_lists.c
base/phast_markov_matrix.c
base/phast_matrix.c
base/phast_memory_handler.c
base/phast_misc.c
base/phast_numerical_opt.c
base/phast_prob_matrix.c
base/phast_prob_vector.c
base/phast_stringsplus.c
base/phast_vector.c
feature/phast_bed.c
feature/phast_category_map.c
feature/phast_genepred.c
feature/phast_gff.c
feature/phast_wig.c
hmm/phast_em.c
hmm/phast_hmm.c
motif/phast_motif.c
motif/phast_tfbs.c
msa/phast_local_alignment.c
msa/phast_maf_block.c
msa/phast_maf.c
msa/phast_msa.c
msa/phast_multi_msa.c
msa/phast_sufficient_stats.c
phylo/phast_bgc_hmm.c
phylo/phast_cons.c
phylo/phast_dgamma.c
phylo/phast_fit_column.c
phylo/phast_fit_em.c
phylo/phast_fit_feature.c
phylo/phast_indel_history.c
phylo/phast_indel_mod.c
phylo/phast_phylo_fit.c
phylo/phast_phylo_p_print.c
phylo/phast_phylo_p.c
phylo/phast_subst_distrib.c
phylo/phast_subst_mods.c
phylo/phast_tree_likelihoods.c
phylo/phast_tree_model.c
phylo/phast_trees.c
phylo_hmm/phast_gap_patterns.c
phylo_hmm/phast_phylo_hmm.c
)

set_target_properties(phast_core PROPERTIES
    OUTPUT_NAME phast
)

target_include_directories(phast_core
    PUBLIC
        "${PHAST_INCLUDE_DIR}"    
        "${PCRE_INCLUDE_DIR}"   
)

# PHAST linear algebra backend (chosen at top-level)
target_link_libraries(phast_core PUBLIC phast::linalg)

# also link to PCRE
target_link_libraries(phast_core PUBLIC phast::linalg phast_pcre)

# Install the library
install(TARGETS phast_core
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
)
