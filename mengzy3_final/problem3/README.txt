NOTE: for AMATH*83 final problem 3
I installed the png16 library for portable network graphics image analysis in C++ on Hyak. 

** You do not require an salloc or sbatch script for this problem - just do it from your /gscratch/amath folder! **   

In the problem directory you will find the files:
128 -rw-r--r-- 1 k8r amath   2972 Jun  5 21:39 README.txt
128 -rw-r--r-- 1 k8r amath  21352 Jun  5 21:43 input_data.png
256 -rw-r--r-- 1 k8r amath 207157 Jun  5 21:43 jki-mm.png
128 -rw-r--r-- 1 k8r amath    490 Jun  5 21:43 cmpl.png-gs
128 -rw-r--r-- 1 k8r amath   3262 Jun  5 21:43 compare_png.cpp
128 -rw-r--r-- 1 k8r amath   6142 Jun  5 21:43 cpp-png-gs.cpp
  1 -rw-r--r-- 1 k8r amath    261 Jun  5 21:43 grayscaleThreaded.cpp

Take a look at the simple compile script:

[k8r@klone-login03 spring2025-final]$ cat cmpl.png-gs 
echo "threaded transform"
g++ -c -std=c++17 -I/gscratch/amath/roche-libs/include/libpng16/ cpp-png-gs.cpp;  
g++ -c -std=c++17 -I/gscratch/amath/roche-libs/include/libpng16/ grayscaleThreaded.cpp
g++ -o xcpp-png-gs cpp-png-gs.o grayscaleThreaded.o -L/gscratch/amath/roche-libs/lib/ -lpng -lm -pthread
echo "correctness"
g++ -c -std=c++17 -I/gscratch/amath/roche-libs/include/libpng16/ compare_png.cpp;  g++ -o xcompare_png compare_png.o -L/gscratch/amath/roche-libs/lib/ -lpng -lm -pthread

Execute the compile script:

[k8r@klone-login03 spring2025-final]$ source cmpl.png-gs 
threaded transform
correctness

Execute the transformation code:

[k8r@klone-login03 spring2025-final]$ ./xcpp-png-gs 
Usage: ./xcpp-png-gs input.png output.png [numThreads]

[k8r@klone-login03 spring2025-final]$ ./xcpp-png-gs input_data.png output_data.png 1
Sequential grayscale time: 1.06775 ms
Threaded grayscale time (1 threads): 4e-05 ms
Wrote output: output_data_seq.png
Wrote output: output_data_thread.png

See that the output file was generated:

[k8r@klone-login03 spring2025-final]$ ls -lstr 
total 3713
 128 -rw-r--r-- 1 k8r amath    3262 Jun  5 20:47 compare_png.cpp
 128 -rw-r--r-- 1 k8r amath   21352 Jun  5 20:47 input_data.png
 256 -rw-r--r-- 1 k8r amath  207157 Jun  5 20:47 jki-mm.png
 128 -rw-r--r-- 1 k8r amath     490 Jun  5 21:31 cmpl.png-gs
 128 -rw-r--r-- 1 k8r amath    6142 Jun  5 21:32 cpp-png-gs.cpp
   1 -rw-r--r-- 1 k8r amath     261 Jun  5 21:34 grayscaleThreaded.cpp
 128 -rw-r--r-- 1 k8r amath   74624 Jun  5 21:34 cpp-png-gs.o
 128 -rw-r--r-- 1 k8r amath    2376 Jun  5 21:34 grayscaleThreaded.o
 128 -rwxr-xr-x 1 k8r amath   61632 Jun  5 21:34 xcpp-png-gs
 128 -rw-r--r-- 1 k8r amath   33440 Jun  5 21:34 compare_png.o
 128 -rwxr-xr-x 1 k8r amath   37768 Jun  5 21:34 xcompare_png
 128 -rw-r--r-- 1 k8r amath   18728 Jun  5 21:35 output_data_seq.png
 128 -rw-r--r-- 1 k8r amath   19061 Jun  5 21:35 output_data_thread.png

Compare the sequential result image to the threaded result image - they should be the same. At the start, the images will clearly differ since you have not implemented the threaded function yet. 

[k8r@klone-login03 spring2025-final]$ ./xcompare_png output_data_seq.png output_data_thread.png 
Difference at (y=37, x=301), channel=0: 1 vs 2
Difference at (y=37, x=301), channel=1: 1 vs 2
Difference at (y=37, x=301), channel=2: 1 vs 2
Difference at (y=37, x=302), channel=0: 15 vs 16
Difference at (y=37, x=302), channel=1: 15 vs 16
Difference at (y=37, x=302), channel=2: 15 vs 16
...
Difference at (y=465, x=292), channel=2: 252 vs 253
Difference at (y=465, x=293), channel=0: 61 vs 62
Difference at (y=465, x=293), channel=1: 61 vs 62
Difference at (y=465, x=293), channel=2: 61 vs 62
Difference at (y=466, x=289), channel=0: 3 vs 4
Difference at (y=466, x=289), channel=1: 3 vs 4
Difference at (y=466, x=289), channel=2: 3 vs 4
❌ Images differ.

BUT, once you get it right, the result will read: 

[k8r@klone-login03 spring2025-final]$ ./xcpp-png-gs input_data.png output_data.png 6
Sequential grayscale time: 1.07587 ms
Threaded grayscale time (6 threads): 0.677422 ms
Wrote output: output_data_seq.png
Wrote output: output_data_thread.png
[k8r@klone-login03 spring2025-final]$ ./xcompare_png output_data_seq.png output_data_thread.png
✅ Images are identical.

