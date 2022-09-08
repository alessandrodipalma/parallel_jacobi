# parallel_jacobi
Parallel implementation of the dp algorithm in C++ standard library, FastFlow and OpenMP.

# Building instructions

To build the project run these commands on a terminal:

```shell
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
```

# Run
The main provides access to all the variants implemented.
You will need to specify:
- **smaller number of workers** on which the problem will be executed
- **largest number of workers** on which the problem will be executed. If these first two parameters are equal, the experiment will run only one time.
- **size**: how big the matrix will be.
- **max iterations**: how many iterations the solver will perform.

Then you will need to set a series of boolean flags to select on which variant you want to run the experiment, in the following order:
- JaFarm in STD C++
- JaFarm in FastFlow
- JaFarm in OpenMP
- JaRed in OpenMP
- JaRed in FastFlow
- Jacobi with only reduce

An example of the command is:
```shell
./build/parallel_jacobi 1 12 512 10000 0 1 0 1 0 0
```
that runs an experiment with 1 to 12 workers, with a system of size 512x512, runs on maximum 10000 iterates, and repeats for Jared FastFlow and JaFarm FastFlow.

