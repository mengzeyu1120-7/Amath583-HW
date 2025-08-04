NOTE: for AMATH*83 final problem 4 

Here, I am gifting you a messy C code you can learn to use. 

It is your job to compile it and execute a strong scaling study with the program on either N=4096 or N=8192 for np = 1,2,4,8,16 and corresponding to rectangular grid sizes of pxq = (1x1),(1x2),(2x2),(2x4),(4x4). For instance your run commands might look like:

 mpirun -np 1  ./xkrazeb 1 1 48 48 4096
 mpirun -np 2  ./xkrazeb 1 2 48 48 4096
 mpirun -np 4  ./xkrazeb 2 2 48 48 4096
 mpirun -np 8  ./xkrazeb 2 4 48 48 4096
 mpirun -np 16 ./xkrazeb 4 4 48 48 4096

You must request the correct resources for this study. For the 4096 study, I used: 

	salloc -A niac -p compute -N 4 --ntasks=16 --cpus-per-task=1 --mem-per-cpu=2G --time=00:45:00

But, you must use one of our amath partitions (don't forget that hyakalloc is your friend).

Prior to running the codes, you will need to load the correct modules. Here is what I did:

module load gcc intel/oneAPI/2023.2.1

Finally, prior to running the code, you will need to compile it. I provided a cheater script for that. Here is how I compiled the code: 

source cmpl.krazeb 

It will create a program called xkrazeb -use it for the problem.

