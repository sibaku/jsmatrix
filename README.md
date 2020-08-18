# JSMatrix

This is a library for matrix and vector computations for JavaScript. Besides the general features such as multiplication, it also provides methods for solving equations and computing decompositions. Currently LU decomposition with partial pivoting, eigenvalue decomposition, SVD and bidiagonalization are implemented.
Importantly, various kinds of views, such as for rows, columns, blocks, transpose and diagonals can be used to manipulate data more naturally. 

All matrices support either typed storage or generic arrays. Special matrices can be implemented by yourself by providing some object with the AbstractMat interface specified in the documentation. Documentation can be found in the docs folder and to view without cloning under https://raw.githack.com/sibaku/jsmatrix/master/doc/index.html

## Usage

The following shows code examples for the basic usage. The documentation provides more in-depth information about the different classes and functions.
You can find the code under examples/basicUsage with additional Latex output in your browser

```javascript
// Basic operations

// Create a new matrix filled with zeros
// This can be done either with calling zeros(...) with a type or using the appropriate type constructor
// Here we create a float array
var a0 = MatF32.zeros(4, 5);

// We can create a string version of this via m.toString
console.log(m.toString(a0));
// 0 0 0 0 0
// 0 0 0 0 0
// 0 0 0 0 0
// 0 0 0 0 0

// We can extract views of the rows and columns
// Extract the third row
var ra0 = m.row(a0, 2);
// Extract the first column
var ca0 = m.col(a0, 0);

// Fill the row with a fixed value
m.fill(ra0, 1.5);
console.log(m.toString(a0));
// 0   0   0   0   0  
// 0   0   0   0   0  
// 1.5 1.5 1.5 1.5 1.5
// 0   0   0   0   0  

// Create a vector and insert it into the column
var v0 = VecF32.from([1, 3, 7, 11]);
// insert the vector
m.insert(ca0, v0);
console.log(m.toString(a0));
// 1  0   0   0   0  
// 3  0   0   0   0  
// 7  1.5 1.5 1.5 1.5
// 11 0   0   0   0  

// A vector is just a matrix with 1 column
// If we want to insert a vector into a row, we need to transpose it
m.insert(m.row(a0, 1), m.transpose(VecF32.from([-1, -2, -3, -4, -5])));
console.log(m.toString(a0));
// 1  0   0   0   0  
// -1 -2  -3  -4  -5 
// 7  1.5 1.5 1.5 1.5
// 11 0   0   0   0  

// Set the values on the diagonal to 1
// The diagonal of the 4x5 matrix is of length 4
// Diagonals are vectors, that is matrices with one column
m.insert(m.diag(a0), VecF32.ones(4));
console.log(m.toString(a0));
// 1  0   0  0   0  
// -1 1   -3 -4  -5 
// 7  1.5 1  1.5 1.5
// 11 0   0  1   0      

// Many methods allow for a last parameter to store the results in
// This can be used to do operations in place in some cases and preallocate buffers to avoid temporaries
// Add values to ca0
m.add(ca0, VecF32.from([-3, -9, 17, 2]), ca0);
console.log(m.toString(a0));
// -2  0   0  0   0  
// -10 1   -3 -4  -5 
// 24  1.5 1  1.5 1.5
// 13  0   0  1   0  

// Accessing values directly is done via the at member method
console.log(a0.at(2, 3));
// 1.5
// The second parameter always defaults to 0, so vectors can be accessed with one index
console.log(m.col(a0, 0).at(1));
// -10

// Setting values follows the same pattern
a0.set(19, 1, 1);
console.log(m.toString(a0));
// -2  0   0  0   0  
// -10 19  -3 -4  -5 
// 24  1.5 1  1.5 1.5
// 13  0   0  1   0  

m.diag(a0).set(20, 2);
console.log(m.toString(a0));
// -2  0   0  0   0  
// -10 19  -3 -4  -5 
// 24  1.5 20 1.5 1.5
// 13  0   0  1   0 

// Create a square matrix and set it to the identity
var a1 = m.setId(MatF32.uninitialized(4, 4));
// Compute its determinant
var deta1 = m.det(a1);
console.log(deta1);
// 1

// Create a new matrix from values
// Values have to be column major
var a2 = MatF32.from(
    [
        1, -2, -3, 0, // first column
        0, 1, 0, 0, // second column
        1, 2, 1, 0, // third column
        2, 0, 2, 1 // fourth column
    ], 4, 4);

// Compute its inverse
var a2i = m.inv(a2);
console.log(m.toString(a2i));
// 0.25 0 -0.25 0 
// -1   1 -1    4 
// 0.75 0 0.25  -2
// 0    0 0     1 

// Multiplying the inverse with the original should be equal to the identity, so we check by summing up the absolute differences
var a3 = m.mult(a2, a2i);
// When we do not need to modify entries, we can use a more economic representation, such as a diagonal matrix
// This only stores the diagonal entries
var diff = m.absSum(m.sub(a3, DiagonalF32.id(4, 4)));
console.log(diff);
// 0

// We can solve an equation Ax = b for b
// Solve will not check, if the matrix is invertible, so for square matrices you need to be sure, they are
var b0 = VecF32.from([0, 2, 4, -2]);
var x0 = m.solve(a2, b0);
console.log(m.toString(x0));
// -1 
// -10
// 5  
// -2 

// This should be the same as a2i*b0, but it is preferrable to use solve instead of computing an inverse
var x1 = m.mult(a2i, b0);
console.log(m.toString(x1));
// -1 
// -10
// 5  
// -2 

// Solve will also work for rectangular matrices, but will only give the closest approximate solution, such that
// |Ax - b| will be minimal
var a4 = MatF32.ones(2, 5);
console.log(m.toString(a4));
// 1 1 1 1 1
// 1 1 1 1 1

var b1 = VecF32.from([0, 1]);
var x2 = m.solve(a4, b1);
console.log(m.toString(x2));
// 0.09999998658895493
// 0.09999999403953552
// 0.09999999403953552
// 0.09999999403953552
// 0.09999999403953552

console.log(m.toString(m.mult(a4, x2)));
// 0.4999999701976776
// 0.4999999701976776

var a5 = MatF32.ones(5, 2);
console.log(m.toString(a5));
// 1 1
// 1 1
// 1 1
// 1 1
// 1 1

var b2 = VecF32.from([0, 1, 2, 3, 4]);
var x3 = m.solve(a5, b2);
console.log(m.toString(x3));
// 0.9999999403953552
// 0.9999999403953552

console.log(m.toString(m.mult(a5, x3)));
// 1.9999998807907104
// 1.9999998807907104
// 1.9999998807907104
// 1.9999998807907104
// 1.9999998807907104

// If you need to compute solutions multiple times, it is advisable to decompose the matrix beforehand
// For square invertible matrices, you can use a pivoted LU decomposition
var plua2 = m.computePLUD(a2);
// The PLUD is used internally for solve as well
// You can either directly use the solve member method or call solvePLU
// This will give the same result as before
var x4 = plua2.solve(b0);
console.log(m.toString(x4));
// -1 
// -10
// 5  
// -2 

// If you need the solution that handles singular or rectangular matrices, you can use the singular value decomposition (SVD)
var svda4 = m.computeSVD(a4);
// This is used internally for solving rectangular matrices and will give the same result as before
var x5 = svda4.solve(b1);
console.log(m.toString(x5));
// 0.09999998658895493
// 0.09999999403953552
// 0.09999999403953552
// 0.09999999403953552
// 0.09999999403953552

// Additional helper functions
// Copying a matrix to one of the same type
var a6 = m.copy(a4);
console.log(m.toString(a6));
// 1 1 1 1 1
// 1 1 1 1 1

// Creating an unitialized matrix with the same type as some other matrix
// Optionally, a new size can be chosen
var a7 = m.similar(a4, 4, 3);
console.log(m.toString(a7));
// 0 0 0
// 0 0 0
// 0 0 0
// 0 0 0

// Setting many values
// Set identity
m.setId(a7);
console.log(m.toString(a7));
// 1 0 0
// 0 1 0
// 0 0 1
// 0 0 0

// Similar zeros, ones, and some values
m.setZero(a7);
console.log(m.toString(a7));
// 0 0 0
// 0 0 0
// 0 0 0
// 0 0 0

m.setOne(a7);
console.log(m.toString(a7));
// 1 1 1
// 1 1 1
// 1 1 1
// 1 1 1

m.fill(a7, 42);
console.log(m.toString(a7));
// 42 42 42
// 42 42 42
// 42 42 42
// 42 42 42


// Map and reduce

// Map each entry to its value minus its row and column index
// Not providing the output as the last argument would create a new matrix
// The output is returned, in this case a7, so a8 = a7
var a8 = m.map(a7, (value, row, col) => value - row - col, a7);
console.log(m.toString(a8));
// 42 41 40
// 41 40 39
// 40 39 38
// 39 38 37

// We can also map to other types, for example, we can convert each entry to a string
// For that we need to provide an appropriate output
var a9 = MatAny.uninitialized(a8.rows(), a8.cols());
m.map(a8, value => "v: " + value, a9);
console.log(m.toString(a9));
// v: 42 v: 41 v: 40
// v: 41 v: 40 v: 39
// v: 40 v: 39 v: 38
// v: 39 v: 38 v: 37

// Reductions provide ways to reduce a matrix expression
// Some reductions, such as sum, min and max are already provided
// Here we find the maximum absolute value of a matrix
var a10 = MatF32.from([
    -10, 2,
    4, 9
], 2, 2);

var maxa10 = m.reduce(a10, (accum, v) => Math.max(accum, Math.abs(v)), 0);
console.log(maxa10);
// 10

// The reduce operation can also take an initial value, which is useful if you don't want to operate on the same data as in the array
// Here we concatenate all entries to a string
var a10s = m.reduce(a10, (accum, v) => accum + v + ";", "");
console.log(a10s);
// -10;2;4;9;

// Reductions can also be comptued per column or row, as seen in the more elaborate SVD example
```

## Example - Rendering and calculating a normal from points

You can find a more complicated example in the examples folder under exampleSVD. Here, geometric functions such as for rotation or creating a camera with a view and perspective matrix are used together with some more advanced reductions and the SVD.

A random point cloud located roughly in a rotated plane is generated. The SVD is used to compute an approximation of the normal of that plane from the random points.

At the end, the points and normal are visualized in a HTML canvas similar to graphic APIs such as OpenGL.

### Preview
![3D Points rendered in 2D and visualized normal computed from SVD](https://github.com/sibaku/jsmatrix/blob/master/examples/exampleSVD/example_svd_normal.jpg?raw=true)
