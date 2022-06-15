## A Custom linear algebra library

---

### Motivation
<p>The reason for building this project was to implement the common linear algebra algorithms that were introduced during the introductory linear algebra course during my freshman year</p>

<br>

### Demo
<br>

$$2*x+5*y+z=20$$
$$4*x+5*y+z=10$$
$$6x+5y=0$$
$$A = \begin{pmatrix} 2 & 5 & 1\\ 4 & 5 & 1\\ 6 & 5 & 0\end{pmatrix}$$
$$B = \begin{pmatrix} 20 \\ 10 \\ 0\end{pmatrix}$$

Solution => $ x=-5,y=6,z=0$


<br/>
<br/>

![./images/gif.gif](./images/gif.gif)

### What can it do?

---

**Data structures**:
  - Matrix(```mat```)
  - ```c
    {
        unsigned int num_rows;
        unsigned int num_cols;
        double **data;
        int isSquare;
    }
    ```
  -LUP Matrix(```mat_lup```)
  - ```c
    {
      mat *L;
      mat *U;
      mat *P;
      unsigned int num_permutations;
    }
    ```

**API**:

- Matrix addition 
  - ```mat* mat_add(mat* a,mat* b)```
- Matrix subtraction
  - ```mat* mat_subtract(mat* a,mat* b)```
- Matrix multiplication
  - ```mat* mat_mul_naive(mat* a,mat* b)``` 
- Finding Row Echelon Form of a Matrix
  - ```mat* mat_ref(mat* a) ```
- Finding Reduced Row Echelon Form of a Matrix
  - ```mat* mat_rref(mat* a)```
- LU decomposition
  - ```mat_lup* mat_LU(mat* a) ```
- Solving linear equations
  - ```mat* solve_linear_LU(mat* a,mat* b)```
- Finding the inverse
  - Using LU method
    - ```mat* mat_inverse(mat* a)```
- Finding the determinant
  - ```double mat_determinant(mat* a)```


### Background and References

---

- I was completely new to C when I started with this project.
- I had done a freshmen introduction to linear algebra so I was kind of familiar with the mathematical concepts.
- References:
  - [3b1b Linear Algebra playlist](https://www.youtube.com/watch?v=kjBOesZCoqc&list=PL0-GT3co4r2y2YErbmuJw2L5tW4Ew2O5B)
  - [Linear Algebra playlist from Khan Academy](https://www.khanacademy.org/math/linear-algebra)
  - [Linear Algebra playlist from MathTheBeautiful](https://www.youtube.com/watch?v=odV3oJOpE8s&list=PLlXfTHzgMRUIqYrutsFXCOmiqKUgOgGJ5)
  - The C programming language(Dennis Ritchie and Brian Kernighan)

### Can this be used in production?

---

<p>Definitely not. This project was build for educational purposes only and this is definitely not a library that one should use for their projects.</p>

### Future of this project

---

- This project is far from complete.
- Some new features to be implemented:
  - A CLI interface for users
  - Replacing current makefile with CMake
  - Computing eigenvalues and eignevectors
  - Adding faster algorithms and better data structures
    - Addition of strassen's algorithm for matrix multiplicaton
    - Dense and sparse matrices
    - More...
