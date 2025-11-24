add_library(phast_pcre STATIC
pcre_chartables.c
pcre_compile.c
pcre_config.c
pcre_dfa_exec.c
pcre_exec.c
pcre_fullinfo.c
pcre_get.c
pcre_globals.c
pcre_info.c
pcre_maketables.c
pcre_newline.c
pcre_ord2utf8.c
pcre_refcount.c
pcre_study.c
pcre_tables.c
pcre_try_flipped.c
pcre_ucd.c
pcre_valid_utf8.c
pcre_version.c
pcre_xclass.c
)

# Mimics "-I." from the Makefile
target_include_directories(phast_pcre
    PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}"
)

# Mimics "-DHAVE_CONFIG_H"
target_compile_definitions(phast_pcre
    PRIVATE
        HAVE_CONFIG_H
)

# "CFLAGS += -fPIC" just for pcre
set_target_properties(phast_pcre PROPERTIES
    POSITION_INDEPENDENT_CODE ON
)
