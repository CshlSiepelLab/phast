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
    target_compile_definitions(pbs_helpers PRIVATE VECLIB ACCELERATE_NEW_LAPACK)
    target_compile_definitions(bd_phylo_hmm_helpers PRIVATE VECLIB ACCELERATE_NEW_LAPACK)
endif()

#
# ---- Help file generation (.help from .help_src) -----------------------
#
# Convert all src/**/.help_src to .help using src/munge-help.sh.
#
file(GLOB_RECURSE HELP_SRC_FILES
  "${CMAKE_SOURCE_DIR}/src/*.help_src"
)

set(GENERATED_HELP_FILES)
foreach(hs ${HELP_SRC_FILES})
  get_filename_component(hs_dir "${hs}" DIRECTORY)
  get_filename_component(hs_we "${hs}" NAME_WE)
  set(out "${hs_dir}/${hs_we}.help")
  list(APPEND GENERATED_HELP_FILES "${out}")

  add_custom_command(
    OUTPUT "${out}"
    COMMAND bash -c
            "bash '${CMAKE_SOURCE_DIR}/src/munge-help.sh' \
             '${hs}' > '${out}'"
    DEPENDS "${hs}" "${CMAKE_SOURCE_DIR}/src/munge-help.sh"
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    VERBATIM
  )
endforeach()

if(GENERATED_HELP_FILES)
  add_custom_target(generate_help ALL
    DEPENDS ${GENERATED_HELP_FILES}
  )
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

  # Ensure help files are generated before compiling the program
  if(TARGET generate_help)
    add_dependencies(${tgt} generate_help)
  endif()

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
