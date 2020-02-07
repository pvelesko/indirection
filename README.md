# Indirection Benchmark for Complex Datatypes
This benchmark is for evaluating the performance of different data structures and ways of performing operations on complex datatypes with one level of indirection. 
Such patterns are common in unstructured grids. Complex data types add another layer of complexity due to being composed of 2 scalar parts - real and imaginary. 
````
for i in num_uniques:
	cmplx0 = values0[indices0[i]];
	cmplx1 = values1[indices1[i]];
	sum += complx0 * cmplx1;
````
## Data Access Pattern
Since there is one level of indirection, the actual data pattern is not known at compile time. 
Consider two extremes: 
1. Index matrix is serial - 0, 1, 2, 3, ...
2. Index matrix is completely random - 64, 23, 1, 33 ...

In real world applications the true layout will fall somewhere in between these two extremes. For instance, triangular domain decomposition and resulting mesh traversal will usually travel between neighboring nodes such that often times 1 or 2 elements are adjacent. 
## Data Types 
### AoS
std::complex is a data structure with 2 scalars - real and imaginary.
Array of std::complex is an Array of Structure data type with the following layout in memory:

![enter image description here](https://lh3.googleusercontent.com/D3aV96d4XomRyY33u_N8KpQEzBOOubq6AbYqQ5BhUofjYn8kxd0r7lwZw7WeKjP1GVer36Hn0L5wPg)

Consider two extremes:
1. Scalar - get the first pair of <real, imag> values, and reuse the cacheline of subsequent operations
2. Random - you get the first pair <real, imag> and discard the rest of the cacheline. 

Pros: guaranteed to get at least 2 elements out of a single load. 
Cons: Have to do shuffles; lower vectorization efficiency

### SoA
Structure of Arrays layout keeps real and imaginary parts separate

![enter image description here](https://lh3.googleusercontent.com/F5OXhuB1r99Ajqknmyw7F-baABCijrCMk8UyEJn-8DWgTHa36nTRWib_Fj3S9GepMkYtRmOhxbsyOQ)

Consider two extremes:
1. Scalar - You generate two steaming loads, potentially hitting more cache banks, achieving higher bandwidth.
2. Random - You get two elements from two cache lines and discard the rest. 

Pros: potentially higher BW, more efficient vectorization
Cons: penalty for completely random access

### AoS vs SoA
While the are advantages to both layouts, such indexing operations are extremely memory bound so vectorization is unlikely to be helpful. 
For large arrays that don't fit in cache, two arrays can result in more TLB penalties. 
All in all, AoS seems like a better layout expect in extremely optimistic cases with minimal random access. 

# The Benchmark
As a proxy to randomness this benchmark takes a % of cache line to fill with sequential elements as input. 
The index matrix is populated in two steps: first fill n-sequential elements, store the remaining elements into a separate array R.
Go over the remaining indices that were left blank and fill with values from R at random. 

For example, for a 64-byte cache line containing 8-byte floats can store 8 values. 75% cache fill ratio would produce the following index matrix
Step 1:
````
indx =  {0, 1, 2,  3,  4,  5,  X1, X2,
         8, 9, 10, 11, 12, 13, X3, X4}
````
 Step 2:
````
indx =  {0, 1, 2,  3,  4,  5,  14, 7,
         8, 9, 10, 11, 12, 13, 6, 15}
````


## Notes
### Vectorization
Benchmark also includes OpenMP threading and forced vectorization using simd. There are 4 calc functions that were using to explore if there's any difference in how complex arithmetic is expressed.
`#pragma omp simd reduction(+: psi_r, psi_i)` successfully forces vectorization with independent real and complex parts.

Forcing vectorization on `std::complex` is also possible but requires the user to provide a reduction operator 
`#pragma omp declare reduction` 
More info: [http://pages.tacc.utexas.edu/~eijkhout/pcse/html/omp-reduction.html](http://pages.tacc.utexas.edu/~eijkhout/pcse/html/omp-reduction.html)
### OpenMP
OpenMP was introduced to primarily exploit hyperthreading. 
For this, Intel OpenMP environment variables were set to limit the application to 1 physical core and 1, 2, or 4 hyperthreads (this was run on KNL): 
`export KMP_HW_SUBSET=1c@2,4t` 
`export OMP_NUM_THREADS=4`
Use 1 core (core #2 specifically, to avoid clashing with OS)  and 4 hyperthreads. 

Scheduling was also explored with no gains above the default settings. 
>  TODO: Explore nested openmp with interlaced operations schedule(static, 1)



