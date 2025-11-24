# one executable per program, linked against vine + PHAST

# 1. Object libraries for shared helpers
add_library(pbs_helpers OBJECT
  prequel/phast_pbs_code.c
  prequel/phast_simplex_grid.c
)

add_library(bd_phylo_hmm_helpers OBJECT
  dless/phast_bd_phylo_hmm.c
)

target_include_directories(pbs_helpers
    PUBLIC
        "${PHAST_INCLUDE_DIR}"   # for <phast/misc.h>, etc.
        "${PCRE_INCLUDE_DIR}"    # if phast_pbs_code.c needs pcre.h too
)

target_include_directories(bd_phylo_hmm_helpers
    PUBLIC
        "${PHAST_INCLUDE_DIR}"   # for <phast/misc.h>, etc.
        "${PCRE_INCLUDE_DIR}"    # if phast_pbs_code.c needs pcre.h too
)

if(APPLE)
#    target_link_libraries(pbs_helpers PUBLIC "-framework Accelerate")
#    target_link_libraries(bd_phylo_hmm_helpers PUBLIC "-framework Accelerate")
    target_compile_definitions(pbs_helpers PRIVATE VECLIB ACCELERATE_NEW_LAPACK)
    target_compile_definitions(bd_phylo_hmm_helpers PRIVATE VECLIB ACCELERATE_NEW_LAPACK)
endif()

set(PHAST_PROG_SOURCES
dless/dless.c
dless/dlessP.c
exoniphy/exoniphy.c
phastCons/phastCons.c
phastMotif/phastMotif.c
phastOdds/phastOdds.c
phyloBoot/phyloBoot.c
phyloFit/phyloFit.c
phyloP/phyloP.c
prequel/pbsDecode.c
prequel/pbsEncode.c
prequel/pbsScoreMatrix.c
prequel/pbsTrain.c
prequel/prequel.c
util/all_dists.c
util/base_evolve.c
util/chooseLines.c
util/clean_genes.c
util/consEntropy.c
util/convert_coords.c
util/display_rate_matrix.c
util/draw_tree.c
util/eval_predictions.c
util/hmm_train.c
util/hmm_tweak.c
util/hmm_view.c
util/indelFit.c
util/indelHistory.c
util/maf_parse.c
util/makeHKY.c
util/modFreqs.c
util/msa_diff.c
util/msa_split.c
util/msa_view.c
util/phast.c
util/phastBias.c
util/refeature.c
util/stringiphy.c
util/tree_doctor.c
util/treeGen.c
)

set(PHAST_EXECUTABLES)

foreach(src ${PHAST_PROG_SOURCES})
  # Derive the executable name from the filename (e.g., phyloFit.c → phyloFit)
  get_filename_component(tgt "${src}" NAME_WE)

  add_executable(${tgt} "${src}")

  # Common link: everything links against the core PHAST library
  target_link_libraries(${tgt} PRIVATE phast_core phast_pcre)

  # Attach shared helper objects where needed
  if(tgt STREQUAL "pbsDecode"
     OR tgt STREQUAL "pbsEncode"
     OR tgt STREQUAL "pbsScoreMatrix"
     OR tgt STREQUAL "pbsTrain"
     OR tgt STREQUAL "prequel")
    target_sources(${tgt} PRIVATE $<TARGET_OBJECTS:pbs_helpers>)
  endif()

  if(tgt STREQUAL "dless" OR tgt STREQUAL "dlessP")
    target_sources(${tgt} PRIVATE $<TARGET_OBJECTS:bd_phylo_hmm_helpers>)
  endif()

  list(APPEND PHAST_EXECUTABLES ${tgt})
endforeach()

# Install all progs into bin/
install(TARGETS
  ${PHAST_EXECUTABLES}
  RUNTIME DESTINATION bin
)
