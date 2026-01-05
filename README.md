# PHAST

**PHAST (Phylogenetic Analysis with Space/Time models)** is a C-based
toolkit for likelihood-based phylogenetic analysis with support for the
identification of conserved elements and other large-scale comparative
genomics analyses.

## Building from source

PHAST now uses a **standard CMake-based build system**.

### Dependencies
- BLAS and LAPACK
- PCRE

(Any system-provided or package-manager-provided implementations are fine.)

### Build and install

```sh
cmake -S . -B build
cmake --build build
cmake --install build
```

By default, binaries are installed to the system prefix (e.g. /usr/local/bin
or a package-manager–controlled prefix).

### Documentation and support

For usage of individual tools, run them with --help.
For questions or bug reports, please use the GitHub issue tracker.

### License

PHAST is released under the BSD 3-Clause license. See LICENSE.txt for
details.
