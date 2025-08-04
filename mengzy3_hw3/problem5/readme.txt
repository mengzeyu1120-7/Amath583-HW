To compile the template BLAS library and test code:

1. Compile Source Files into Object Files：
# Compile non-templated functions (double precision)
g++ -std=c++17 -c -fPIC ref_daxpy.cpp -o ref_daxpy.o
g++ -std=c++17 -c -fPIC ref_dgemv.cpp -o ref_dgemv.o
g++ -std=c++17 -c -fPIC ref_dgemm.cpp -o ref_dgemm.o

# Compile templated functions (generic types)
g++ -std=c++17 -c -fPIC ref_axpyt.cpp -o ref_axpyt.o
g++ -std=c++17 -c -fPIC ref_gemvt.cpp -o ref_gemvt.o
g++ -std=c++17 -c -fPIC ref_gemmt.cpp -o ref_gemmt.o

2. Create the Shared Library (librefBLAS.so)：
g++ -shared -o librefBLAS.so ref_axpyt.o ref_gemvt.o ref_gemmt.o

3. Compile the Test Program：
g++ -std=c++17 test_templates.cpp -L. -lrefBLAS -o template_test

4. run the test programme：
LD_LIBRARY_PATH=. ./template_test