cmake_minimum_required(VERSION 3.16)

project(phast
  VERSION 1.9.1
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

# ---- PCRE handling: prefer system, fallback to bundled ------------------

# Try to find a system/Homebrew PCRE first
find_path(PCRE_INCLUDE_DIR pcre.h)
find_library(PCRE_LIBRARY pcre)
if(NOT PCRE_INCLUDE_DIR OR NOT PCRE_LIBRARY)
  message(FATAL_ERROR "PCRE not found. Install libpcre development package.")
endif()

message(STATUS "Using system PCRE: ${PCRE_LIBRARY}")

add_library(pcre::pcre UNKNOWN IMPORTED)
set_target_properties(pcre::pcre PROPERTIES
  IMPORTED_LOCATION "${PCRE_LIBRARY}"
  INTERFACE_INCLUDE_DIRECTORIES "${PCRE_INCLUDE_DIR}"
)

add_library(phast_pcre INTERFACE)
target_link_libraries(phast_pcre INTERFACE pcre::pcre)
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
else()
    # Use system BLAS/LAPACK elsewhere
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)

    # CMake 3.16 often provides *_LIBRARIES variables reliably (targets are less consistent)
    target_link_libraries(phast_linalg INTERFACE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    target_compile_definitions(phast_linalg INTERFACE PHAST_USE_SYSTEM_LAPACK)
    add_compile_definitions(PHAST_USE_SYSTEM_LAPACK)
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
