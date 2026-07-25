cmake_minimum_required(VERSION 3.18)

project(phast
  VERSION 1.9.5
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

# ---- PCRE handling: require system  -----------------------------------

# Try to find a system/Homebrew PCRE first
find_path(PCRE_INCLUDE_DIR pcre.h)
find_library(PCRE_LIBRARY pcre)
if(NOT PCRE_INCLUDE_DIR OR NOT PCRE_LIBRARY)
  message(FATAL_ERROR "PCRE not found. Install libpcre development package.")
endif()

message(STATUS "Using system PCRE: ${PCRE_LIBRARY}")

# Export-friendly: do NOT depend on an external CMake target like pcre::pcre.
add_library(phast_pcre INTERFACE)
target_include_directories(phast_pcre INTERFACE "${PCRE_INCLUDE_DIR}")
target_link_libraries(phast_pcre INTERFACE "${PCRE_LIBRARY}")

add_library(phast::pcre ALIAS phast_pcre)

# ---- Linear algebra backend (single source of truth) -------------------

add_library(phast_linalg INTERFACE)
add_library(phast::linalg ALIAS phast_linalg)

if(APPLE)
    # Use Accelerate framework on macOS
    find_library(ACCELERATE_FRAMEWORK Accelerate REQUIRED)
    target_link_libraries(phast_linalg INTERFACE ${ACCELERATE_FRAMEWORK})

    # Keep your existing compile-time switches, but attach them to the linalg "interface"
    target_compile_definitions(phast_linalg INTERFACE VECLIB ACCELERATE_NEW_LAPACK)

    # Record the backend so the installed package config knows not to re-resolve BLAS.
    set(PHAST_BLAS_BACKEND "Accelerate")
else()
    # Use system BLAS/LAPACK elsewhere.
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)

    # Link the imported targets (CMake >= 3.18) rather than the ${*_LIBRARIES} absolute
    # paths. Absolute paths get baked into the exported target and break downstream
    # consumers whose environment resolves BLAS differently (e.g. conda-forge ships
    # only libopenblas.so.0, not the unversioned libopenblas.so). With imported targets,
    # the installed phastConfig re-resolves BLAS in the consumer's environment via
    # find_dependency() below.
    target_link_libraries(phast_linalg INTERFACE BLAS::BLAS LAPACK::LAPACK)
    target_compile_definitions(phast_linalg INTERFACE PHAST_USE_SYSTEM_LAPACK)
    add_compile_definitions(PHAST_USE_SYSTEM_LAPACK)

    set(PHAST_BLAS_BACKEND "System")
endif()


# ---- Project layout ----------------------------------------------------

# public headers 
set(PHAST_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")

add_subdirectory(src)
add_subdirectory(src/lib)

# ---- Installation (for Homebrew, etc.) ---------------------------------
# Install public headers
install(DIRECTORY "${PHAST_INCLUDE_DIR}/"
        DESTINATION include
        FILES_MATCHING PATTERN "*.h")

# Export these interface targets too, because phast_core depends on them
install(TARGETS phast_linalg phast_pcre
  EXPORT phastTargets
)

include(CMakePackageConfigHelpers)

# Where to install the config files
set(PHAST_CMAKE_DIR "lib/cmake/phast")

# Install the exported targets file
install(EXPORT phastTargets
  NAMESPACE phast::
  DESTINATION ${PHAST_CMAKE_DIR}
)

# Create and install phastConfig.cmake and phastConfigVersion.cmake
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/phastConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/phastConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/phastConfig.cmake"
  INSTALL_DESTINATION ${PHAST_CMAKE_DIR}
)

install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/phastConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/phastConfigVersion.cmake"
  DESTINATION ${PHAST_CMAKE_DIR}
)
