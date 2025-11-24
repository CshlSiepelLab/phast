cmake_minimum_required(VERSION 3.16)

project(phast
  VERSION 1.9
  LANGUAGES C)

# Use C99 
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

add_compile_definitions(PHAST_VERSION=\"${PROJECT_VERSION}\")
add_compile_definitions(PHAST_HOME=\"${CMAKE_SOURCE_DIR}\")

# these are used everywhere
add_compile_options(
    -g -fno-inline -Wall -fno-strict-aliasing
)

if(APPLE)
    add_compile_definitions(VECLIB ACCELERATE_NEW_LAPACK)
#   - VECLIB tells external_libs.h to use Accelerate
#   - ACCELERATE_NEW_LAPACK selects the new (non-CLAPACK) interface
else()
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    add_compile_definitions(PHAST_LAPACK_GENERIC)
endif()


# ---- Project layout ----------------------------------------------------

# public headers 
set(PHAST_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")

add_subdirectory(src)
add_subdirectory(src/lib)
add_subdirectory(src/lib/pcre)

# ---- PCRE handling  ----------------------------------------------------
# use system / Homebrew PCRE
find_path(PCRE_INCLUDE_DIR pcre.h)
find_library(PCRE_LIBRARY pcre)

if (NOT PCRE_INCLUDE_DIR OR NOT PCRE_LIBRARY)
    message(FATAL_ERROR "Requested system PCRE but could not find headers or library")



# ---- Installation (for Homebrew, etc.) ---------------------------------
# Install public headers
install(DIRECTORY "${PHAST_INCLUDE_DIR}/"
        DESTINATION include
        FILES_MATCHING PATTERN "*.h")
