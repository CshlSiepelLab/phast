make
R CMD SHLIB ../*/*.o
mv ../base/complex_matrix.so ~/phast/lib/librphast.so