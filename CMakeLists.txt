cmake_minimum_required(VERSION 3.16)

project(phast
  VERSION 1.7
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


# ---- Project layout ----------------------------------------------------

# public headers 
set(PHAST_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")
set(PCRE_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/src/lib/pcre")

add_subdirectory(src)
add_subdirectory(src/lib)
add_subdirectory(src/lib/pcre)

# ---- Installation (for Homebrew, etc.) ---------------------------------
# Install public headers
install(DIRECTORY "${PHAST_INCLUDE_DIR}/"
        DESTINATION include
        FILES_MATCHING PATTERN "*.h")

