/* global renderMathInElement */
import {
    MatF32,
    VecF32,
    MatAny,
    DiagonalF32,
} from '../../jsmatrix.js';

import * as m from '../../jsmatrix.js';

// Helper function to map a matrix to latex
function toLtx(a) {
    const mstr = m.map(a, x => x === undefined ? "\\_" : x instanceof Number ? x.toPrecision(3) : x, MatAny.uninitialized(a.rows(), a.cols()));

    const rows = m.rowreduce(mstr, x => m.toArray(x).join(" & "), MatAny.uninitialized(a.rows(), 1));

    return '\\begin{pmatrix}' + m.toArray(rows).join("\\\\") + '\\end{pmatrix}';

}

// Append a div to the document and update its math content
function write(s) {
    const el = document.createElement('div');
    el.innerText = s;
    document.body.appendChild(
        el
    );
    renderMathInElement(el);
}

window.onload = function () {

    // Basic operations

    // Create a new matrix filled with zeros
    // This can be done either with calling zeros(...) with a type or using the appropriate type constructor
    // Here we create a float array
    var a0 = MatF32.zeros(4, 5);

    write(`Create a zero matrix of size \\((4,5)\\):\n  \\(\\mathbf{A} = ${toLtx(a0)}\\)`);
    // We can create a string version of this via m.toString
    console.log(m.toString(a0));
    // 0 0 0 0 0
    // 0 0 0 0 0
    // 0 0 0 0 0
    // 0 0 0 0 0

    // We can extract views of the rows and columns
    // Extract the third row
    var ra0 = m.row(a0, 2);
    write(`Extracting the 3rd row as an expression:\n \\(${toLtx(ra0)}\\)`);
    // Extract the first column
    var ca0 = m.col(a0, 0);
    write(`Extracting the 1st column as an expression:\n \\(${toLtx(ca0)}\\)`);

    // Fill the row with a fixed value
    m.fill(ra0, 1.5);
    write(`Filling just the extracted row with the value \\(1.5\\):\n \\(\\mathbf{A} = ${toLtx(a0)}\\)`);

    console.log(m.toString(a0));
    // 0   0   0   0   0  
    // 0   0   0   0   0  
    // 1.5 1.5 1.5 1.5 1.5
    // 0   0   0   0   0  

    // Create a vector and insert it into the column
    var v0 = VecF32.from([1, 3, 7, 11]);
    write(`Create a vector:\n \\(${toLtx(v0)}\\)`);

    // insert the vector
    m.insert(ca0, v0);

    write(`Insert that vector into the extracted column:\n \\(\\mathbf{A} = ${toLtx(a0)}\\)`);

    console.log(m.toString(a0));
    // 1  0   0   0   0  
    // 3  0   0   0   0  
    // 7  1.5 1.5 1.5 1.5
    // 11 0   0   0   0  

    // A vector is just a matrix with 1 column
    // If we want to insert a vector into a row, we need to transpose it
    m.insert(m.row(a0, 1), m.transpose(VecF32.from([-1, -2, -3, -4, -5])));
    write(`Insert a (column) vector into the second row by transposing it:\n \\(\\mathbf{A} = ${toLtx(a0)}\\)`);

    console.log(m.toString(a0));
    // 1  0   0   0   0  
    // -1 -2  -3  -4  -5 
    // 7  1.5 1.5 1.5 1.5
    // 11 0   0   0   0  

    // Set the values on the diagonal to 1
    // The diagonal of the 4x5 matrix is of length 4
    // Diagonals are vectors, that is matrices with one column
    m.insert(m.diag(a0), VecF32.ones(4));
    write(`Insert a vector of ones into the diagonal:\n \\(\\mathbf{A} = ${toLtx(a0)}\\)`);

    console.log(m.toString(a0));
    // 1  0   0  0   0  
    // -1 1   -3 -4  -5 
    // 7  1.5 1  1.5 1.5
    // 11 0   0  1   0      

    // Many methods allow for a last parameter to store the results in
    // This can be used to do operations in place in some cases and preallocate buffers to avoid temporaries
    // Add values to ca0
    m.add(ca0, VecF32.from([-3, -9, 17, 2]), ca0);
    write(`Add the vector \\(${toLtx(m.transpose(VecF32.from([-3, -9, 17, 2])))}^T \\) to another and store the result in the first one.In this case, the first vector is the view of the first row, thus changing the matrix in place. This can also be used to pre allocate and reuse memory: \n \\(\\mathbf{A} = ${toLtx(a0)} \\)`);

    console.log(m.toString(a0));
    // -2  0   0  0   0  
    // -10 1   -3 -4  -5 
    // 24  1.5 1  1.5 1.5
    // 13  0   0  1   0  

    // Accessing values directly is done via the at member method
    console.log(a0.at(2, 3));
    write(`Elements can be accessed with the at method. \\(\\mathbf{A}_{2,3} = ${a0.at(2, 3)}\\)`);

    // 1.5
    // The second parameter always defaults to 0, so vectors can be accessed with one index
    console.log(m.col(a0, 0).at(1));
    // -10

    // Setting values follows the same pattern
    a0.set(19, 1, 1);
    write(`Inserting is done with the set method. Insert 19 at \\(\\mathbf{A}_{1,1}\\):\n \\(\\mathbf{A} = ${toLtx(a0)}\\)`);

    console.log(m.toString(a0));
    // -2  0   0  0   0  
    // -10 19  -3 -4  -5 
    // 24  1.5 1  1.5 1.5
    // 13  0   0  1   0  

    m.diag(a0).set(20, 2);
    write(`Vectors are just matrices with one row. Set and at will default to 0 as their last parameter, so you can index a vector with just one index. Here we set the third entry on the diagonal to 20:\n \\(\\mathbf{A} = ${toLtx(a0)}\\)`);

    console.log(m.toString(a0));
    // -2  0   0  0   0  
    // -10 19  -3 -4  -5 
    // 24  1.5 20 1.5 1.5
    // 13  0   0  1   0 

    // Create a square matrix and set it to the identity
    // Alternatively, use the MatF32.id method
    var a1 = m.setId(MatF32.uninitialized(4, 4));
    write(`Create a new identity matrix:\n \\(\\mathbf{B} = ${toLtx(a1)}\\)`);

    // Compute its determinant
    var deta1 = m.det(a1);
    write(`Compute its determinant with the det function:\n \\(\\det \\mathbf{B} = ${deta1}\\)`);

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

    write(`Create a new matrix from data. Data is stored in column major order:\n \\(\\mathbf{C} = ${toLtx(a2)}\\)`);

    // Compute its inverse
    var a2i = m.inv(a2);
    write(`Compute its inverse with the inv function:\n \\(\\mathbf{C}^{-1} = ${toLtx(a2i)}\\)`);

    console.log(m.toString(a2i));
    // 0.25 0 -0.25 0 
    // -1   1 -1    4 
    // 0.75 0 0.25  -2
    // 0    0 0     1 

    // Multiplying the inverse with the original should be equal to the identity, so we check by summing up the absolute differences
    var a3 = m.mult(a2, a2i);

    write(`Sanity check for the inverse:\n \\(\\mathbf{C}\\mathbf{C}^{-1} = ${toLtx(a3)}\\)`);

    // When we do not need to modify entries, we can use a more economic representation, such as a diagonal matrix
    // This only stores the diagonal entries
    var diff = m.absSum(m.sub(a3, DiagonalF32.id(4, 4)));
    write(`Sometimes, a read only matrix is enough. Here we compute the absolute sum of differences between the previous computation and the identity matrix, which is here just stored as a diagonal matrix expression, that only stores the diagonal entries. The sum of absolute differences is: \\(${diff}\\)`);

    console.log(diff);
    // 0

    // We can solve an equation Ax = b for b
    // Solve will not check, if the matrix is invertible, so for square matrices you need to be sure, they are
    var b0 = VecF32.from([0, 2, 4, -2]);
    var x0 = m.solve(a2, b0);
    write(`Solving with the solve method.\n\\(\\mathbf{C}\\mathbf{x} = \\mathbf{b}\\)\n\\(${toLtx(a2)}\\mathbf{x} = ${toLtx(b0)}\\):\n \\(\\mathbf{x}= ${toLtx(x0)}\\)`);

    console.log(m.toString(x0));
    // -1 
    // -10
    // 5  
    // -2 

    // This should be the same as a2i*b0, but it is preferrable to use solve instead of computing an inverse
    var x1 = m.mult(a2i, b0);
    write(`Comparing that result with solving with the inverse:\n \\(\\mathbf{C}^{-1}\\mathbf{b} = ${toLtx(x1)}\\)`);

    console.log(m.toString(x1));
    // -1 
    // -10
    // 5  
    // -2 

    // Solve will also work for rectangular matrices, but will only give the closest approximate solution, such that
    // |Ax - b| will be minimal
    var a4 = MatF32.ones(2, 5);
    write(`Create a rectangular matrix:\n \\(\\mathbf{D} = ${toLtx(a4)}\\)`);

    console.log(m.toString(a4));
    // 1 1 1 1 1
    // 1 1 1 1 1

    var b1 = VecF32.from([0, 1]);
    var x2 = m.solve(a4, b1);
    write(`Solving also works for general rectangular matrices. In that case, the SVD is used to find the closest approximate solution with respect to the 2-norm`
        + `Solve for \\(\\mathbf{b}=${toLtx(m.transpose(b1))}^T\\):\n \\(\\mathbf{x} = ${toLtx(x2)}\\)`);

    console.log(m.toString(x2));
    // 0.09999998658895493
    // 0.09999999403953552
    // 0.09999999403953552
    // 0.09999999403953552
    // 0.09999999403953552

    write(`Check solution:\n \\(\\mathbf{D}\\mathbf{x}=\\)\n\\(${toLtx(a4)}${toLtx(x2)} = ${toLtx(m.mult(a4, x2))}\\)`);
    console.log(m.toString(m.mult(a4, x2)));
    // 0.4999999701976776
    // 0.4999999701976776

    var a5 = MatF32.ones(5, 2);
    write(`Create a differently shaped rectangular matrix:\n \\(\\mathbf{E} = ${toLtx(a5)}\\)`);

    console.log(m.toString(a5));
    // 1 1
    // 1 1
    // 1 1
    // 1 1
    // 1 1

    var b2 = VecF32.from([0, 1, 2, 3, 4]);
    var x3 = m.solve(a5, b2);
    write(`Solve for \\(\\mathbf{b}=${toLtx(m.transpose(b2))}^T\\):\n \\(\\mathbf{x} = ${toLtx(x3)}\\)`);

    console.log(m.toString(x3));
    // 0.9999999403953552
    // 0.9999999403953552

    write(`Check solution:\n \\(\\mathbf{D}\\mathbf{x}=\\)\n\\(${toLtx(a5)}${toLtx(x3)} = ${toLtx(m.mult(a5, x3))}\\)`);
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
    write(`If you need to solve multiple times with the same matrix, consider computing a decompositon before. For square invertible matrices, a good idea is to use the pivoted LU decomposition by calling computePLUD or PLUD.compute:\n`
        + `We solve \\(\\mathbf{C}\\mathbf{x} = \\mathbf{b}\\)\n\\(${toLtx(a2)}\\mathbf{x} = ${toLtx(b0)}\\)\n`
        + ` \\(\\mathbf{x}= ${toLtx(x4)}\\)`);

    console.log(m.toString(x4));
    // -1 
    // -10
    // 5  
    // -2 

    // If you need the solution that handles singular or rectangular matrices, you can use the singular value decomposition (SVD)
    var svda4 = m.computeSVD(a4);

    // This is used internally for solving rectangular matrices and will give the same result as before
    var x5 = svda4.solve(b1);
    write(`For general and possibly singular matrices, you can use a singular value decomposition. This is done with the computeSVD function or SVD.compute:\n`
        + `We solve \\(\\mathbf{E}\\mathbf{x} = \\mathbf{b}\\)\n\\(${toLtx(a2)}\\mathbf{x} = ${toLtx(b1)}\\)\n`
        + ` \\(\\mathbf{x}= ${toLtx(x5)}\\)`);
    console.log(m.toString(x5));
    // 0.09999998658895493
    // 0.09999999403953552
    // 0.09999999403953552
    // 0.09999999403953552
    // 0.09999999403953552

    // Additional helper functions
    // Copying a matrix to one of the same type
    var a6 = m.copy(a4);
    write(`Use the copy function to create a copy of a matrix with the same type:\n \\(\\mathbf{F} = ${toLtx(a6)}\\)`);
    console.log(m.toString(a6));
    // 1 1 1 1 1
    // 1 1 1 1 1

    // Creating an unitialized matrix with the same type as some other matrix
    // Optionally, a new size can be chosen
    var a7 = m.similar(a4, 4, 3);
    write(`Use the similar function to create an uninitialized matrix with the same type as the input but possibly different dimensions:\n \\(\\mathbf{G} = ${toLtx(a7)}\\)`);
    console.log(m.toString(a7));
    // 0 0 0
    // 0 0 0
    // 0 0 0
    // 0 0 0

    // Setting many values
    // Set identity
    m.setId(a7);
    write(`Use the setId function to set a matrix to identity:\n \\(\\mathbf{G} = ${toLtx(a7)}\\)`);
    console.log(m.toString(a7));
    // 1 0 0
    // 0 1 0
    // 0 0 1
    // 0 0 0

    // Similar zeros, ones, and some values
    m.setZero(a7);
    write(`Use the setZero function to set a matrix to zero:\n \\(\\mathbf{G} = ${toLtx(a7)}\\)`);
    console.log(m.toString(a7));
    // 0 0 0
    // 0 0 0
    // 0 0 0
    // 0 0 0

    m.setOne(a7);
    write(`Use the setOne function to set a matrix to ones:\n \\(\\mathbf{G} = ${toLtx(a7)}\\)`);
    console.log(m.toString(a7));
    // 1 1 1
    // 1 1 1
    // 1 1 1
    // 1 1 1

    m.fill(a7, 42);
    write(`Use the fill function to set a matrix to a specified value:\n \\(\\mathbf{G} = ${toLtx(a7)}\\)`);
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
    write(`Use the map function to map each element of a matrix to a new one, if needed a new matrix is created`
        + ` Here we subtract the row and colum index from each value:\n \\(\\mathbf{G} = ${toLtx(a7)}\\)`);
    console.log(m.toString(a8));
    // 42 41 40
    // 41 40 39
    // 40 39 38
    // 39 38 37

    // We can also map to other types, for example, we can convert each entry to a string
    // For that we need to provide an appropriate output
    var a9 = MatAny.uninitialized(a8.rows(), a8.cols());
    write(`We can also use non numeric types. These are captured by using a non typed Array or by using the MatAny factory`
        + `Here we create a generic unitialized matrix:\n \\(\\mathbf{H} = ${toLtx(a9)}\\)`);

    m.map(a8, value => "v: " + value, a9);
    write(`We can then use the generic matrix as an output for mappings. Here we prefix every value of \\(\\mathbf{G}\\) with the string "v: " :\n \\(\\mathbf{H} = ${toLtx(a9)}\\)`);
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
    write(`Reductions provide ways to reduce a matrix expression.` +
        `Some reductions, such as sum, min and max are already provided.` +
        `Here we find the maximum absolute value of a matrix:\n \\(\\mathbf{I} = ${toLtx(a10)}\\)`
        + `\nThe maximum is: \\(${maxa10}\\)`);
    console.log(maxa10);
    // 10

    // The reduce operation can also take an initial value, which is useful if you don't want to operate on the same data as in the array
    // Here we concatenate all entries to a string
    var a10s = m.reduce(a10, (accum, v) => accum + v + ";", "");
    write(`The reduce operation can also take an initial value, which is useful if you don't want to operate on the same data as in the array.`
        + `Here we concatenate all entries to a string: ${a10s}`);
    console.log(a10s);
    // -10;2;4;9;

    // Reductions can also be comptued per column or row, as seen in the more elaborate SVD example

};