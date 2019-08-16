// Linter comments
/* global Float32Array, Float64Array, Int8Array, Uint8Array, Uint8ClampedArray,Int16Array,Uint16Array,Int32Array,Uint32Array,BigInt64Array,BigUint64Array */
'use strict';

/**
 * Storage types used for matrices
 */
var f32 = Float32Array;
var f64 = Float64Array;
var i8 = Int8Array;
var ui8 = Uint8Array;
var ui8c = Uint8ClampedArray;
var i16 = Int16Array;
var ui16 = Uint16Array;
var i32 = Int32Array;
var ui32 = Uint32Array;
var i64 = BigInt64Array;
var ui64 = BigUint64Array;
var any = Array;

//*****************************************************
function optional(v, defaultValue) {
    return v !== undefined ? v : defaultValue;
}
//*****************************************************
/** void => int
 * 
 * Generic function
    @name intVoidFunction
    @function
    @returns {number}
*/
//*****************************************************
/** (i,j) => * 
 * 
 * Function to access elements
    @name atFunction
    @function
    @param {number} i - The row
    @param {number} [j=0] - The column
    @returns {*} The matrix element at location (i,j)
*/
//*****************************************************
/** (v,i,j) => AbstractMat 
 * 
 * Function to set element values. 
    @name setFunction
    @function
    @param {number} v - The new value
    @param {number} i - The row
    @param {number} [j=0] - The column
    @returns {AbstractMat} this
*/
//*****************************************************
/** void => function
 * 
 * Generic function
    @name funcFunction
    @function
    @returns {function}
*/
//*****************************************************
/**
 * An abstract matrix
 * @typedef {Object} AbstractMat
 * @property {intVoidFunction} rows - The number of rows
 * @property {intVoidFunction} cols - The number of cols
 * @property {atFunction} at - The function to retrieve elements from. Access is only valid in the ranges i=[0,rows-1] and j=[0,cols-1]
 * @property {setFunction} set - The function to set elements. Access is only valid in the ranges i=[0,rows-1] and j=[0,cols-1]. Some matrices may not allow setting certain entries
 * @property {funcFunction} type - The function to get the underlying type of the matrix. This should either be Array or any TypedArray
 */
//*****************************************************
/**
 * An array-based matrix class. It can work with strided data.
 * 
 * @extends AbstractMat
 */
class Mat {
    /**
     * 
     * @param {Array | TypedArray} data Input Array. Can be typed or generic.
     * @param {number} rows Number of rows
     * @param {number} cols Number of colums
     * @param {number} [innerStride=1] Units from one element to the next. By default 1
     * @param {number} [outerStride=rows] Units from one column to the next. By default the number of rows
     */
    constructor(data, rows, cols, innerStride, outerStride) {
        /** @private */
        this._data = data;
        /** @private */
        this._rows = rows;
        /** @private */
        this._cols = cols;
        /** @private */
        this._innerStride = optional(innerStride, 1);
        /** @private */
        this._outerStride = optional(outerStride, this._rows);
    }
    /**
    * 
    * @param {Array | TypedArray} data Input Array. Can be typed or generic.
    * @param {number} rows Number of rows
    * @param {number} cols Number of colums
    * @param {number} [innerStride=1] Units from one element to the next. By default 1
    * @param {number} [outerStride=rows] Units from one column to the next. By default the number of rows
    *
    * @returns {Mat} A new 
    * @see constructor
    */
    static new(data, rows, cols, innerStride, outerStride) {
        return new Mat(data, rows, cols, innerStride, outerStride);
    }

    rows() {
        return this._rows;
    }
    cols() {
        return this._cols;
    }


    index(i, j) {
        j = j !== undefined ? j : 0;

        return this._innerStride * (i + j * this._outerStride);
    }
    at(i, j) {

        var idx = this.index(i, j);

        if (idx > this._data.length) {
            throw "Trying to access data outside of range";
        }

        return this._data[idx];
    }

    set(v, i, j) {
        var idx = this.index(i, j);
        if (idx > this._data.length) {
            throw "Trying to access data outside of range";
        }
        this._data[idx] = v;
        return this;
    }

    type() {
        return this._data.constructor;
    }


}
//*****************************************************
/**
 * Helper to create special matrices with a given underlying array type
 */
class TypedMatFactory {
    /**
     * 
     * @param {function} type - The base type to be used. This should be either Array or a TypedArray
     */
    constructor(type) {
        this._type = type;
    }
    /**
     * Creates a new TypedMatFactory object
     * 
     * @param {function} type - The base type to be used. This should be either Array or a TypedArray
     * @returns {TypedMatFactory}
     * 
     * @see TypedMatFactory.constructor
     */
    static new(type) {
        return new TypedMatFactory(type);
    }

    static newFromMat(m) {
        return new TypedMatFactory(m.type !== undefined && m.type instanceof Function ? m.type() : Array);
    }
    uninitialized(rows, cols) {
        var n = rows * cols;
        var data = new this._type(n);
        return Mat.new(data, rows, cols);
    }

    from(data, rows, cols) {
        var n = rows * cols;
        var datat = new this._type(n);

        for (let i = 0; i < n; i++) {
            datat[i] = data[i];
        }
        return Mat.new(data, rows, cols);
    }
    all(v, rows, cols) {
        var n = rows * cols;
        var data = new this._type(n);
        for (let i = 0; i < n; i++) {
            data[i] = v;
        }
        return Mat.new(data, rows, cols);
    }

    zeros(rows, cols) {
        return this.all(0, rows, cols);
    }

    ones(rows, cols) {
        return this.all(1, rows, cols);
    }

    id(rows, cols) {
        cols = optional(cols, rows);
        var m = this.zeros(rows, cols);

        var dm = diag(m);
        fill(dm, 1);

        return m;
    }

    copy(m) {
        var c = m.cols();
        var r = m.rows();
        var n = c * r;

        var data = new this._type(n);

        var idx = 0;
        for (let j = 0; j < c; j++) {
            for (let i = 0; i < r; i++) {
                data[idx] = m.at(i, j);
                idx++;
            }
        }
        return Mat.new(data, r, c);
    }

    rand(rows, cols) {
        var r = rows;
        var c = cols;
        var n = c * r;

        var data = new this._type(n);

        var idx = 0;
        for (let j = 0; j < c; j++) {
            for (let i = 0; i < r; i++) {
                data[idx] = Math.random();
                idx++;
            }
        }
        return Mat.new(data, r, c);
    }
}
//*****************************************************
/**
 * A TypedMatFactory to create Float32Array based matrices
 */
var MatF32 = new TypedMatFactory(Float32Array);
/**
 * A TypedMatFactory to create Float64Array based matrices
 */
var MatF64 = new TypedMatFactory(Float64Array);
/**
 * A TypedMatFactory to create generic Array based matrices
 */
var MatAny = new TypedMatFactory(Array);
/**
 * A TypedMatFactory to create Int8Array based matrices
 */
var MatI8 = new TypedMatFactory(Int8Array);
/**
 * A TypedMatFactory to create Uint8Array based matrices
 */
var MatUI8 = new TypedMatFactory(Uint8Array);
/**
 * A TypedMatFactory to create Uint8ClampedArray based matrices
 */
var MatUI8Clamped = new TypedMatFactory(Uint8ClampedArray);

/**
 * A TypedMatFactory to create Int16Array based matrices
 */
var MatI16 = new TypedMatFactory(Int16Array);
/**
 * A TypedMatFactory to create Uint16Array based matrices
 */
var MatUI16 = new TypedMatFactory(Uint16Array);

/**
 * A TypedMatFactory to create Int32Array based matrices
 */
var MatI32 = new TypedMatFactory(Int32Array);
/**
 * A TypedMatFactory to create Uint32Array based matrices
 */
var MatUI32 = new TypedMatFactory(Uint32Array);

/**
 * A TypedMatFactory to create BigInt64Array based matrices
 */
var MatI64 = new TypedMatFactory(BigInt64Array);
/**
 * A TypedMatFactory to create BigUint64Array based matrices
 */
var MatUI64 = new TypedMatFactory(BigUint64Array);
//*****************************************************
/**
 * Helper to create special vectors with a given underlying array type
 */
class TypedVecFactory {
    /**
     * 
     * @param {function} type - The base type to be used. This should be either Array or a TypedArray
     */
    constructor(type) {
        /** @private */
        this._type = type;
        /** @private */
        this._factory = TypedMatFactory.new(type);
    }
    /**
     * Creates a new TypedVecFactory object
     * 
     * @param {function} type - The base type to be used. This should be either Array or a TypedArray
     * @returns {TypedVecFactory}
     * 
     * @see TypedVecFactory.constructor
     */
    static new(type) {
        return new TypedVecFactory(type);
    }

    static newFromMat(m) {
        return new TypedVecFactory(m.type !== undefined && m.type instanceof Function ? m.type() : Array);
    }

    from(data, n) {
        n = optional(n, data.length);
        return this._factory.from(data, n, 1);
    }
    uninitialized(n) {
        return this._factory.uninitialized(n, 1);
    }

    all(v, n) {
        return this._factory.all(v, n, 1);
    }

    zeros(n) {
        return this.all(0, n);
    }

    ones(n) {
        return this.all(1, n);
    }

    id(n) {
        return this._factory.id(n);
    }

    copy(m) {

        if (!isVec(m)) {
            throw "Vector factory can't copy from non-vector";
        }

        return this._factory.copy(m);
    }

    rand(n) {

        return this._factory.rand(n, 1);
    }
}
//*****************************************************
/**
 * A TypedVecFactory to create Float32Array based vectors
 */
var VecF32 = new TypedVecFactory(Float32Array);
/**
 * A TypedVecFactory to create Float64Array based vectors
 */
var VecF64 = new TypedVecFactory(Float64Array);
/**
 * A TypedVecFactory to create generic Array based vectors
 */
var VecAny = new TypedVecFactory(Array);
/**
 * A TypedVecFactory to create Int8Array based vectors
 */
var VecI8 = new TypedVecFactory(Int8Array);
/**
 * A TypedVecFactory to create Uint8Array based vectors
 */
var VecUI8 = new TypedVecFactory(Uint8Array);
/**
 * A TypedVecFactory to create Uint8ClampedArray based vectors
 */
var VecUI8Clamped = new TypedVecFactory(Uint8ClampedArray);

/**
 * A TypedVecFactory to create Int16Array based vectors
 */
var VecI16 = new TypedVecFactory(Int16Array);
/**
 * A TypedVecFactory to create Uint16Array based vectors
 */
var VecUI16 = new TypedVecFactory(Uint16Array);

/**
 * A TypedVecFactory to create Int32Array based vectors
 */
var VecI32 = new TypedVecFactory(Int32Array);
/**
 * A TypedVecFactory to create Uint32Array based vectors
 */
var VecUI32 = new TypedVecFactory(Uint32Array);

/**
 * A TypedVecFactory to create BigInt64Array based vectors
 */
var VecI64 = new TypedVecFactory(BigInt64Array);
/**
 * A TypedVecFactory to create BigUint64Array based vectors
 */
var VecUI64 = new TypedVecFactory(BigUint64Array);
//*****************************************************
/**
 * A diagonal matrix that only stores the diagonal as data.
 * 
 * @extends AbstractMat
 */
class Diagonal {
    constructor(data, rows, cols, stride) {
        this._data = data;
        this._rows = rows !== undefined ? rows : this._data.length;
        this._cols = cols !== undefined ? cols : this._rows;
        this._stride = stride !== undefined ? stride : 1;

        if (this._data.length < this._stride * Math.min(this._rows, this._cols)) {
            throw "Not enough data for diagonal construction";
        }
    }

    static new(data, rows, cols, stride) {
        return new Diagonal(data, rows, cols, stride);
    }

    static fromMat(m, rows, cols, stride) {
        if (!isVec(m)) {
            throw "Diagonal only takes vectors as input";
        }

        var data = toArray(m);
        return this.new(data, rows, cols, stride);
    }

    rows() {
        return this._rows;
    }

    cols() {
        return this._cols;
    }
    at(i, j) {
        if (i !== j) {
            return 0;
        }

        return this._data[this._stride * i];
    }

    set(v, i, j) {
        if (i !== j) {
            throw "Trying to set offdiagonal element of diagonal matrix";
        }
        this._data[this._stride * i] = v;
        return this;
    }

    type() {
        return this._data.constructor;
    }
}
//*****************************************************
/**
 * Helper to create special diagonal matrices with a given underlying array type
 */
class TypedDiagonalFactory {
    constructor(type) {
        this._type = type;
    }

    static new(type) {
        return new TypedDiagonalFactory(type);
    }

    static newFromMat(m) {
        return new TypedDiagonalFactory(m.type !== undefined && m.type instanceof Function ? m.type() : Array);
    }


    from(data, rows, cols) {
        cols = optional(cols, rows);
        var n = Math.min(rows, cols);
        var datat = new this._type(n);
        for (let i = 0; i < n; i++) {
            datat[i] = data[i];
        }

        return Diagonal.new(data, rows, cols);
    }
    uninitialized(rows, cols) {
        cols = cols !== undefined ? cols : rows;
        var n = Math.min(rows, cols);
        var data = new this._type(n);
        return Diagonal.new(data, rows, cols);
    }
    all(v, rows, cols) {
        cols = optional(cols, rows);
        var n = Math.min(rows, cols);
        var data = new this._type(n);
        for (let i = 0; i < n; i++) {
            data[i] = v;
        }
        return Diagonal.new(data, rows, cols);
    }
    id(rows, cols) {
        return this.all(1, rows, cols);
    }

    rand(rows, cols) {
        cols = cols !== undefined ? cols : rows;
        var n = Math.min(rows, cols);
        var data = new this._type(n);
        for (let i = 0; i < data.length; i++) {
            data[i] = Math.random();
        }
        return Diagonal.new(data, rows, cols);
    }
}
//*****************************************************
/**
 * A TypedDiagonalFactory to create Float32Array based diagonal matrices
 */
var DiagonalF32 = new TypedDiagonalFactory(Float32Array);

/**
 * A TypedDiagonalFactory to create Float64Array based diagonal matrices
 */
var DiagonalF64 = new TypedDiagonalFactory(Float64Array);

/**
 * A TypedDiagonalFactory to create generic Array based diagonal matrices
 */
var DiagonalAny = new TypedDiagonalFactory(Array);

/**
 * A TypedDiagonalFactory to create Int8Array based matrices
 */
var DiagonalI8 = new TypedDiagonalFactory(Int8Array);
/**
 * A TypedDiagonalFactory to create Uint8Array based matrices
 */
var DiagonalUI8 = new TypedDiagonalFactory(Uint8Array);
/**
 * A TypedDiagonalFactory to create Uint8ClampedArray based matrices
 */
var DiagonalUI8Clamped = new TypedDiagonalFactory(Uint8ClampedArray);

/**
 * A TypedDiagonalFactory to create Int16Array based matrices
 */
var DiagonalI16 = new TypedDiagonalFactory(Int16Array);
/**
 * A TypedDiagonalFactory to create Uint16Array based matrices
 */
var DiagonalUI16 = new TypedDiagonalFactory(Uint16Array);

/**
 * A TypedDiagonalFactory to create Int32Array based matrices
 */
var DiagonalI32 = new TypedDiagonalFactory(Int32Array);
/**
 * A TypedDiagonalFactory to create Uint32Array based matrices
 */
var DiagonalUI32 = new TypedDiagonalFactory(Uint32Array);

/**
 * A TypedDiagonalFactory to create BigInt64Array based matrices
 */
var DiagonalI64 = new TypedDiagonalFactory(BigInt64Array);
/**
 * A TypedDiagonalFactory to create BigUint64Array based matrices
 */
var DiagonalUI64 = new TypedDiagonalFactory(BigUint64Array);
//*****************************************************
/**
 * A view of the diagonal of a matrix
 * 
 * @extends AbstractMat
 */
class DiagonalView {
    /**
     * 
     * @param {*} m 
     */
    constructor(m) {
        this._m = m;
        this._size = Math.min(this._m.rows(), this._m.cols());
    }

    static new(m) {
        return new DiagonalView(m);
    }

    rows() {
        return this._size;
    }

    cols() {
        // row vector
        return 1;
    }

    at(i) {
        return this._m.at(i, i);
    }

    set(v, i) {
        this._m.set(v, i, i);
        return this;
    }

    type() {
        return this._m.type();
    }
}
//*****************************************************
/**
 * A view of a block of a matrix
 * 
 * @extends AbstractMat
 */
class BlockView {
    constructor(m, i, j, rows, cols) {
        this._m = m;
        this._i = i;
        this._j = j;
        this._rows = rows;
        this._cols = cols;
    }
    static new(m, i, j, rows, cols) {
        return new BlockView(m, i, j, rows, cols);
    }

    rows() {
        return this._rows;
    }

    cols() {
        return this._cols;
    }

    at(i, j) {
        j = j !== undefined ? j : 0;
        return this._m.at(i + this._i, j + this._j);
    }

    set(v, i, j) {
        j = j !== undefined ? j : 0;
        return this._m.set(v, i + this._i, j + this._j);
    }

    type() {
        return this._m.type();
    }
}
//*****************************************************
/**
 * A view of a column of a matrix. This will be a column vector
 * 
 * @extends BlockView
 */
class ColumnView extends BlockView {
    constructor(m, j) {
        var r = m.rows();
        super(m, 0, j, r, 1);
    }

    static new(m, j) {
        return new ColumnView(m, j);
    }
}
//*****************************************************
/**
 * A view of a row of a matrix. This will be a row vector
 * 
 * @extends BlockView
 */
class RowView extends BlockView {
    constructor(m, i) {
        var c = m.cols();
        super(m, i, 0, 1, c);
    }

    static new(m, i) {
        return new RowView(m, i);
    }
}
//*****************************************************
/**
 * Transpose view of a matrix. Functions as a transpose without copy
 * 
 * @extends AbstractMat
 */
class TransposeView {
    constructor(m) {
        this._m = m;
    }
    static new(m) {
        return new TransposeView(m);
    }
    rows() {
        return this._m.cols();
    }

    cols() {
        return this._m.rows();
    }

    at(i, j) {
        j = j !== undefined ? j : 0;
        return this._m.at(j, i);
    }

    set(v, i, j) {
        j = j !== undefined ? j : 0;
        return this._m.set(v, j, i);
    }

    type() {
        return this._m.type();
    }
}
//*****************************************************
/**
 * Specifies the type of triangular matrix.
 * 
 * All elements not included in the triangular part are considered zero
 */
var TriangularMode = {
    /**
     * Upper diagonal matrix including the diagonal
     */
    UPPER: 0,
    /**
     * Upper diagonal matrix including the diagonal
     */
    LOWER: 1,
    /**
     * Strictly upper diagonal matrix excluding the diagonal
     */
    STRICTLY_UPPER: 2,
    /**
     * Strictly lower diagonal matrix excluding the diagonal
     */
    STRICTLY_LOWER: 3,
    /**
     * Upper diagonal matrix. The diagonal entries entries are 1
     */
    UNIT_UPPER: 4,
    /**
    * Lower diagonal matrix. The diagonal entries entries are 1
     */
    UNIT_LOWER: 5
};
//*****************************************************
/**
 * A view of a triangular portion of a given matrix.
 * 
 * The type of portion may be specified by one of the values in TriangularMode. 
 * Entries that are part of the view may be changed in the underlying matrix. 
 * When specifying TriangularMode.UNIT_UPPER or TriangularMode.UNIT_LOWER, the diagonal accessor will return 1, but setting the diagonal is not possible.
 * @see TriangularMode
 * 
 * @extends AbstractMat
 */
class TriangularView {
    constructor(m, mode) {
        this._m = m;
        this._mode = mode;
    }

    static new(m, mode) {
        return new TriangularView(m, mode);
    }

    rows() {
        return this._m.rows();
    }

    cols() {
        return this._m.cols();
    }

    checkIndex(i, j) {
        return (this._mode === TriangularMode.UPPER && i <= j) ||
            (this._mode === TriangularMode.LOWER && i >= j) ||
            ((this._mode === TriangularMode.STRICTLY_UPPER || this._mode === TriangularMode.UNIT_UPPER) && i < j) ||
            ((this._mode === TriangularMode.STRICTLY_LOWER || this._mode === TriangularMode.UNIT_LOWER) && i > j);
    }

    at(i, j) {
        switch (this._mode) {
            case TriangularMode.UPPER:
                if (i > j) {
                    return 0;
                }
                else {
                    return this._m.at(i, j);
                }
            case TriangularMode.LOWER:
                if (i < j) {
                    return 0;
                }
                else {
                    return this._m.at(i, j);
                }
            case TriangularMode.STRICTLY_UPPER:
                if (i >= j) {
                    return 0;
                }
                else {
                    return this._m.at(i, j);
                }
            case TriangularMode.STRICTLY_LOWER:
                if (i <= j) {
                    return 0;
                }
                else {
                    return this._m.at(i, j);
                }

            case TriangularMode.UNIT_LOWER:
                if (i < j) {
                    return 0;
                }
                else if (i === j) {
                    return 1;
                }
                else {
                    return this._m.at(i, j);
                }
            case TriangularMode.UNIT_UPPER:
                if (i >= j) {
                    return 0;
                }
                else if (i === j) {
                    return 1;
                }
                else {
                    return this._m.at(i, j);
                }
            default:
                break;
        }

    }
    set(v, i, j) {
        if (!this.checkIndex(i, j)) {
            throw "Trying to set value outside triangular part";
        }

        this._m.set(v, i, j);
        return this;
    }

    type() {
        return this._m.type();
    }

    solve(b, out) {
        if (this._mode === TriangularMode.UPPER ||
            this._mode === TriangularMode.STRICTLY_UPPER ||
            this._mode === TriangularMode.UNIT_UPPER) {
            return solveUpperTriagonal(this, b, out);
        }
        else {
            return solveLowerTriagonal(this, b, out);
        }
    }
}
//*****************************************************
/**
 * Represents a row permutation matrix without explicit storage of the matrix.
 * 
 * Permutations are represented by a permutation table given as an array.
 * The underlying type can be set explicitely.
 * 
 * @extends AbstractMat
 */
class RowPermutation {
    constructor(data, rows, cols, type) {
        /** @private */
        this._data = data;
        /** @private */
        this._rows = optional(rows, this._data.length);
        /** @private */
        this._cols = optional(cols, this._rows);
        this._type = optional(type, currentDefaultType);
    }

    static new(data, rows, cols, type) {
        return new RowPermutation(data, rows, cols, type);
    }

    at(i, j) {
        // permute rows then map to identity
        i = this._data[i];
        j = optional(j, 0);

        return i === j ? 1 : 0;
    }

    set(/**v,i,j */) {
        throw "Permuation matrix is immutable";
    }

    rows() {
        return this._cols;
    }

    cols() {
        return this._cols;
    }

    type() {
        return this._type;
    }
}
//*****************************************************
/**
 * A view of the minor submatrix of a given matrix.
 * 
 * The minor is the matrix with one deleted row and column
 * 
 * @extends AbstractMat
 */
class MinorView {
    constructor(m, i, j) {
        if (i >= m.rows() || j >= m.cols()) {
            throw "Trying to remove row/column outside of source";
        }

        /** @private */
        this._m = m;
        /** @private */
        this._i = i;
        /** @private */
        this._j = j;
    }
    static new(m, i, j) {
        return new MinorView(m, i, j);
    }

    rows() {
        return this._m.rows() - 1;
    }
    cols() {
        return this._m.cols() - 1;
    }

    at(i, j) {
        // skip row, if it was deleted
        if (i >= this._i) {
            i += 1;
        }
        // skip col, if it was deleted
        if (j >= this._j) {
            j += 1;
        }

        return this._m.at(i, j);
    }

    set(v, i, j) {
        // skip row, if it was deleted
        if (i >= this._i) {
            i += 1;
        }
        // skip col, if it was deleted
        if (j >= this._j) {
            j += 1;
        }

        return this._m.set(v, i, j);
    }

    type() {
        return this._m.type();
    }

}
//*****************************************************
/**
 * A padded view to embed another matrix in.
 * 
 * Its size needs to be bigger then the size of the underlying matrix. Offdiagonal and diagonal constant values may be specified
 * 
 * @extends AbstractMat
 */
class PaddedView {
    constructor(m, rows, cols, offDiagonalValue, diagonalValue) {
        /** @private */
        this._m = m;
        /** @private */
        this._rows = rows;
        /** @private */
        this._cols = cols;

        if (this._rows < m.rows() || this._cols < m.cols()) {
            throw "Padded View needs to be at least as big as source";
        }
        /** @private */
        this._offDiagonalValue = optional(offDiagonalValue, 0);
        /** @private */
        this._diagonalValue = optional(diagonalValue, 1);
    }

    static new(m, rows, cols, offDiagonalValue, diagonalValue) {
        return new PaddedView(m, rows, cols, offDiagonalValue, diagonalValue);
    }

    rows() {
        return this._rows;
    }

    cols() {
        return this._cols;
    }

    at(i, j) {
        if (i < this._m.rows() && j < this._colsm.cols()) {
            return this._m.at(i, j);
        }

        if (i === j) {
            return this._diagonalValue;
        } else {
            return this._offDiagonalValue;
        }
    }
    set(v, i, j) {
        if (i < this._m.rows() && j < this._colsm.cols()) {
            return this._m.set(v, i, j);
        }

        throw "Trying to set values in padding";
    }

    type() {
        return this._m.type();
    }
}
//*****************************************************
/**
 * Creates an unitialized matrix similar to the given one
 * 
 * The copy will have the same size and underlying type as m
 * @param {AbstractMat} m - The input matrix
 * @returns {AbstractMat} An unitialized matrix of the same size and type as m
 */
function similar(m, rows, cols) {
    rows = rows !== undefined ? rows : m.rows();
    cols = cols !== undefined ? cols : m.cols();


    return TypedMatFactory.new(m.type()).uninitialized(rows, cols);
}
//*****************************************************
/**
 * Creates a copy of the given matrix.
 * 
 * The copy will have the same underlying type as m
 * @param {AbstractMat} m - The input matrix
 * @returns {Mat} A copy of m
 */
function copy(m) {
    return TypedMatFactory.new(m.type()).copy(m);
}
//*****************************************************
/**
 * Checks whether a matrix is a vector.
 * 
 * Only column vectors are considered vectors, that is matrices with one column
 * 
 * @param {AbstractMat} m - The matrix to check
 * @returns {boolean} True, if the matrix is a column vector, false otherwise
 */
function isVec(m) {
    return m.cols() === 1;
}/**
 * Checks whether a matrix is a row vector, that is matrices with one row
 * 
 * @param {AbstractMat} m - The matrix to check
 * @returns {boolean} True, if the matrix is a row vector, false otherwise
 */
function isRowVec(m) {
    return m.rows() === 1;
}
//*****************************************************
//*****************************************************
/**
 * Adds a value to each entry in the matrix
 * @param {AbstractMat} a - First matrix
 * @param {number | string} v - The value to add
 * @returns {AbstractMat} a
 */
function addScalar(a, v) {

    var r = a.rows();
    var c = a.cols();

    for (let j = 0; j < c; j++) {
        for (let i = 0; i < r; i++) {
            a.set(a.at(i, j) + v, i, j);
        }
    }

    return a;
}
//*****************************************************
/**
 * Adds matrix b to matrix a.
 * @param {AbstractMat} a - First matrix
 * @param {AbstractMat} b - Second matrix
 * @param {AbstractMat} [out] - The output matrix. If not specified, a new matrix will be created.
 * @returns {AbstractMat} a+b
 */
function add(a, b, out) {


    if (a.cols() !== b.cols() || a.rows() != b.rows()) {
        throw "Trying to add input of different sizes";
    }
    var r = a.rows();
    var c = a.cols();

    out = out !== undefined ? out : similar(a);

    for (let j = 0; j < c; j++) {
        for (let i = 0; i < r; i++) {
            let v = a.at(i, j) + b.at(i, j);
            out.set(v, i, j);
        }
    }

    return out;
}
//*****************************************************
/**
 * Subtracts matrix b from matrix a.
 * @param {AbstractMat} a - First matrix
 * @param {AbstractMat} b - Second matrix
 * @param {AbstractMat} [out] - The output matrix. If not specified, a new matrix will be created
 * @returns {AbstractMat} a-b
 */
function sub(a, b, out) {

    if (a.cols() !== b.cols() || a.rows() != b.rows()) {
        throw "Trying to add input of different sizes";
    }

    var r = a.rows();
    var c = a.cols();

    out = out !== undefined ? out : similar(a);

    for (let j = 0; j < c; j++) {
        for (let i = 0; i < r; i++) {
            let v = a.at(i, j) - b.at(i, j);
            out.set(v, i, j);
        }
    }

    return out;
}
//*****************************************************
/**
 * Sets all entries in a matrix equal to a given value
 * 
 * @param {AbstractMat} m - The matrix to change
 * @param {*} v - The value
 * @returns {AbstractMat} m
 */
function fill(a, v) {
    var r = a.rows();
    var c = a.cols();
    for (let j = 0; j < c; j++) {
        for (let i = 0; i < r; i++) {
            a.set(v, i, j);
        }
    }

    return a;
}
//*****************************************************
/**
 * Computes a componenwise multiplication of two matrices of the same size
 * 
 * Computes (i,j) => a.at(i,j)*b.at(i,j)
 * 
 * @param {AbstractMat} a - The first input matrix
 * @param {AbstractMat} b - The second input matrix
 * @param {AbstractMat} out - The output matrix. If not specified, a new matrix will be created
 */
function cwiseMult(a, b, out) {

    return map(a, (v, row, col) => v * b.at(row, col), out);
}
//*****************************************************
/**
 * Computes a componenwise division of two matrices of the same size
 * 
 * Computes (i,j) => a.at(i,j)/b.at(i,j)
 * 
 * @param {AbstractMat} a - The first input matrix
 * @param {AbstractMat} b - The second input matrix
 * @param {AbstractMat} out - The output matrix. If not specified, a new matrix will be created
 */
function cwiseDiv(a, b, out) {

    return map(a, (v, row, col) => v / b.at(row, col), out);
}
//*****************************************************
/**
 * Creates a new DiagonalView of a given matrix
 * 
 * @param {AbstractMat} m - The input matrix
 * @returns {DiagonalView} A DiagonalView of m
 * 
 * @see DiagonalView
 */
function diag(m) {
    return DiagonalView.new(m);
}
//*****************************************************
/**
 * Creates a new BlockView of a given matrix
 * 
 * @param {AbstractMat} m - The input matrix 
 * @param {number} i - The start row
 * @param {number} j - The start column
 * @param {number} rows - The number of rows
 * @param {number} cols - The number of columns
 * @returns {BlockView} A BlockView of m
 * 
 * @see BlockView
 */
function block(m, i, j, rows, cols) {
    rows = optional(rows, m.rows() - i);
    cols = optional(cols, m.cols() - j);
    return BlockView.new(m, i, j, rows, cols);
}
//*****************************************************
/**
 * Creates a new RowView for the i-th row of a given matrix
 * 
 * @param {AbstractMat} m - The input matrix 
 * @param {number} i - The row
 * @returns {ColumnView} A RowView of the i-th row of m
 * 
 * @see RowView
 */
function row(m, i) {
    return RowView.new(m, i);
}
//*****************************************************
/**
 * Creates a new ColumnView for the j-th column of a given matrix
 * 
 * @param {AbstractMat} m - The input matrix 
 * @param {number} j - The column
 * @returns {ColumnView} A ColumnView of the j-th column of m
 * 
 * @see ColumnView
 */
function col(m, j) {
    return ColumnView.new(m, j);
}
//*****************************************************
/**
 * Multiplies a MxN matrix A with a NxP matrix B resulting in a MxP matrix C, so C = A*B
 * 
 * @param {AbstractMat} a - Matrix to multiply on the left
 * @param {AbstractMat} b - Matrix to multiply on the right
 * @param {AbstractMat} [out] - The output matrix. If not specified, a new matrix will be created. Should not be the same as a or b
 * @returns {AbstractMat} The output matrix a*b
 */
function mult(a, b, out) {
    var n = a.rows();
    var m = a.cols();

    if (m !== b.rows()) {
        throw "Incompatible matrix dimensions";
    }
    var p = b.cols();

    out = out !== undefined ? out : similar(a, n, p);

    if (out.rows() !== n || out.cols() !== p) {
        throw "Output dimension does not match input";
    }

    for (let j = 0; j < p; j++) {
        for (let i = 0; i < n; i++) {
            var s = 0;
            for (let k = 0; k < m; k++) {
                s += a.at(i, k) * b.at(k, j);
            }
            out.set(s, i, j);
        }
    }

    return out;
}
//*****************************************************
/**
 * Negates all entries in a given matrix
 * 
 * @param {AbstractMat} m - The source matrix
 * @param {AbstractMat} [out] - The output matrix. If not specified, a new matrix will be created. May be the same as m
 * @returns {AbstractMat} The output matrix 
 */
function neg(m, out) {
    var r = m.rows();
    var c = m.cols();

    out = out !== undefined ? out : similar(m);

    for (let j = 0; j < c; j++) {
        for (let i = 0; i < r; i++) {
            out.set(-m.at(i, j), i, j);
        }
    }

    return out;
}
//*****************************************************
/**
 * Maps each element to its absolute value
 * 
 * @param {AbstractMat} a - The input matrix
 * @param {AbstractMat} [out] - The output matrix. If not specified, a new matrix will be created
 * @returns {AbstractMat} The output matrix
 */
function abs(a, out) {
    return map(a, x => Math.abs(x), out);
}
//*****************************************************
/**
 * Sets all entries in a matrix to zero
 * 
 * @param {AbstractMat} m - The matrix to change
 * @returns {AbstractMat} m
 */
function setZero(m) {
    return fill(m, 0);
}
//*****************************************************
/**
 * Sets all entries in a matrix to one
 * 
 * @param {AbstractMat} m - The matrix to change
 * @returns {AbstractMat} m
 */
function setOne(m) {
    return fill(m, 1);
}
//*****************************************************
/**
 * Sets the given MxN matrix to the identity. 
 * 
 * Each entriy (i,j) will be replace by 1, if i=j and 0 else
 * 
 * @param {AbstractMat} m - The matrix to change 
 * @returns {AbstractMat} m
 */
function setId(m) {
    var r = m.rows();
    var c = m.cols();

    for (let j = 0; j < c; j++) {
        for (let i = 0; i < r; i++) {
            let v = i === j ? 1 : 0;
            m.set(v, i, j);
        }
    }

    return m;

}
//*****************************************************
/**
 * Create a transpose view of the given matrix
 * 
 * @param {AbstractMat} m - The input matrix
 * @returns {TransposeView} A TransposeView of the input
 */
function transpose(m) {
    return TransposeView.new(m);
}
//*****************************************************
/**
 * Maps each element of the given matrix to a new value computed by the given function
 * 
 * @param {AbstractMat} m - The input matrix
 * @param {function} f - (value,row,col) => *. Function to map a value at position (row,col) to a new value
 * @param {AbstractMat} [out] - The output matrix. If none is specified, a new matrix will be created
 */
function map(m, f, out) {
    var r = m.rows();
    var c = m.cols();

    out = out !== undefined ? out : similar(m);

    for (let j = 0; j < c; j++) {
        for (let i = 0; i < r; i++) {
            out.set(f(m.at(i, j), i, j), i, j);
        }
    }

    return out;
}
//*****************************************************
/**
 * Creates a copy of the given matrix with a new underlying storage type
 * 
 * @param {AbstractMat} m - The input matrix
 * @param {*} type - The type type of data storage to be used for the converted matrix
 * @returns {Mat} A new matrix with the given type and contents converted from the input matrix
 */
function convert(m, type) {
    var out = TypedMatFactory.new(type).uninitialized(m.rows(), m.cols());

    return map(m, x => x, out);

}
//*****************************************************
/**
 * Reduces a matrix to a single value.
 * 
 * The reducer function is called for each element. The return value is passed to the next invocation as the accumulator.
 * 
 * If initValue is given, this will be the initial value of the accumulator.
 * If initValue is not specified the initial value of the accumulator will be the matrix's first element. 
 * In that case, the first element the reducer is called with is the second matrix element.
 * 
 * @param {AbstractMat} m - The input matrix
 * @param {*} f - (accum,v,row,col,m) => *. The reduce function
 * @param {*} [initValue] Initial value used for accumulation
 * @returns {*} The accumulated result
 */
function reduce(m, f, initValue) {
    var no_init = initValue === undefined;

    var r = m.rows();
    var c = m.cols();
    if (no_init && r * c === 0) {
        throw "Calling reduce on empty input without initial value";
    }

    var accum = no_init ? m.at(0, 0) : initValue;


    for (let j = 0; j < c; j++) {
        for (let i = 0; i < r; i++) {
            if (no_init) {
                // skip first, if no initvalue is given
                no_init = false;
                continue;
            }
            accum = f(accum, m.at(i, j), i, j, m);
        }
    }

    return accum;
}
//*****************************************************
/**
 * Computes the sum of all matrix elements
 * @param {AbstractMat} m - The input matrix
 * @returns {number} The sum of all elements
 */
function sum(m) {
    return reduce(m, (acc, v) => acc + v);
}
//*****************************************************
/**
 * Computes the absolute sum of all matrix elements
 * @param {AbstractMat} m - The input matrix
 * @returns {number} The absolute sum of all elements
 */
function absSum(m) {
    return reduce(m, (acc, v) => acc + Math.abs(v), 0);
}
//*****************************************************
/**
 * Computes the sum of all squared matrix elements
 * @param {AbstractMat} m - The input matrix
 * @returns {number} The sum of all squared elements
 */
function sqrSum(m) {
    return reduce(m, (acc, v) => acc + v * v, 0);
}
//*****************************************************
/**
 * Computes the trace of the given matrix.
 * 
 * The trace is the sum of the diagonal elements
 * 
 * @param {AbstractMat} m - The input matrix
 * @returns {number} The trace of the matrix
 */
function trace(m) {
    return sum(diag(m));
}
//*****************************************************
/**
 * Computes the rank of a matrix
 * 
 * The rank is computed using the SVD
 * 
 * @param {AbstractMat} m - The input matrix
 * @returns {number} The rank of the matrix
 */
function rank(m) {
    var svd = computeSVD(m);
    return rankSVD(svd);
}
//*****************************************************
/**
 * Computes the rank of a matrix given its SVD
 * 
 * @param {SVD} m - The SVD of a matrix
 * @returns {number} The rank of the matrix
 */
function rankSVD(svd) {
    // count non-zero singular values
    var s = svd.S;

    for (let i = 0; i < s.rows(); i++) {
        if (s.at(i) < 1E-7) {
            return i;
        }
    }

    return s.rows();
}
//*****************************************************
/**
 * Computes the condition number of a matrix
 * 
 * It is computed as the ration sigma_{max}/sigma_{min}, where the sigmas are the singular values of the matrix
 * 
 * @param {AbstractMat} m - The input matrix
 * @returns {number} The condition number
 */
function cond(m) {
    var svd = computeSVD(m);
    return condSVD(svd);
}
//*****************************************************
/**
 * Computes the condition number of a matrix given its SVD
 * 
 * It is computed as the ration sigma_{max}/sigma_{min}, where the sigmas are the singular values of the matrix
 * 
 * @param {SVD} m - The SVD of the input matrix
 * @returns {number} The condition number
 */
function condSVD(svd) {

    var s = svd.S;
    var smax = s.at(0);
    var smin = s.at(s.rows() - 1);

    return smax / smin;

}
//*****************************************************
/**
 * Computes the maximum of all matrix elements
 * @param {AbstractMat} m - The input matrix
 * @returns {number} The maximum of all elements
 */
function max(m) {
    return reduce(m, (acc, v) => Math.max(acc, v));
}
//*****************************************************
/**
 * Computes the minimum of all matrix elements
 * @param {AbstractMat} m - The input matrix
 * @returns {number} The minimum of all elements
 */
function min(m) {
    return reduce(m, (acc, v) => Math.min(acc, v));
}
//*****************************************************
/**
 * @typedef {Object} ArgResult
 * @property {*} v - The value
 * @property {number} row - The row in which v was found
 * @property {number} col - The column in which v was found
 */
//*****************************************************
/**
 * Computes the location of the maximum of all matrix elements
 * @param {AbstractMat} m - The input matrix
 * @returns {ArgResult} The maximum of all elements and where it is
 */
function argmax(m) {
    return reduce(m, (acc, v, i, j) => {
        if (v > acc.v) {
            acc.v = v;
            acc.row = i;
            acc.col = j;
        }
        return acc;
    }, { v: -Infinity, row: -1, col: -1 });
}
//*****************************************************
/**
 * Computes the location of the minimum of all matrix elements
 * @param {AbstractMat} m - The input matrix
 * @returns {ArgResult} The minimum of all elements and where it is
 */
function argmin(m) {
    return reduce(m, (acc, v, i, j) => {
        if (v < acc.v) {
            acc.v = v;
            acc.row = i;
            acc.col = j;
        }
        return acc;
    }, { v: Infinity, row: -1, col: -1 });
}
//*****************************************************
/**
 * Produces an Array from all values in the matrix in column-major order
 * @param {AbstractMat} m - The input matrix
 * @returns {Array} An array with all matrix entries
 */
function toArray(m) {
    var c = m.cols();
    var r = m.rows();
    var result = [];
    for (let j = 0; j < c; j++) {
        for (let i = 0; i < r; i++) {
            result.push(m.at(i, j));
        }
    }

    return result;
}
//*****************************************************
/**
 * Creates an identiy matrix with the currently set default type
 * 
 * @param {number} rows - The number of rows
 * @param {number} [cols=rows] - The number of cols  
 * @returns {Diagonal} The identity matrix
 */
function id(rows, cols, type) {
    type = optional(type, currentDefaultType);
    return TypedDiagonalFactory.new(type).id(rows, cols);
}
//*****************************************************
/**
 * Creates a zero matrix with the currently set default type
 * 
 * @param {number} rows - The number of rows
 * @param {number} cols - The number of cols 
 * @returns {Mat} The zero matrix
 */
function zeros(rows, cols, type) {
    type = optional(type, currentDefaultType);
    return TypedMatFactory.new(type).zeros(rows, cols);
}
//*****************************************************
/**
 * Creates a random matrix
 * 
 * @param {number} rows - The number of rows
 * @param {number} cols - The number of cols 
 * @returns {Mat} The zero matrix
 */
function rand(rows, cols, type) {
    type = optional(type, currentDefaultType);
    return TypedMatFactory.new(type).rand(rows, cols);
}
//*****************************************************
/**
 * Creates a ones matrix with the currently set default type
 * 
 * @param {number} rows - The number of rows
 * @param {number} cols - The number of cols 
 * @returns {Mat} The ones matrix
 */
function ones(rows, cols, type) {
    type = optional(type, currentDefaultType);
    return TypedMatFactory.new(type).ones(rows, cols);
}
//*****************************************************
/**
 * Creates a new uninitialized matrix with the currently set default type
 * 
 * @param {number} rows - The number of rows
 * @param {number} cols - The number of cols 
 * @returns {Mat} An unitialized matrix
 */
function mat(rows, cols, type) {
    type = optional(type, currentDefaultType);
    return TypedMatFactory.new(type).uninitialized(rows, cols);
}
//*****************************************************
/**
 * Creates a Mat from data with a given underlying type
 * @param {Array | TypedArray} data - The input data
 * @param {number} rows - The rows 
 * @param {number} cols - The columns
 * @param {Object} [type] - The underlying type. Defaults to the current default type
 * @returns {Mat} A mat with a typed copy of the given data
 */
function from(data, rows, cols, type) {
    type = optional(type, currentDefaultType);
    return TypedMatFactory.new(type).copy(Mat.new(data, rows, cols));
}
//*****************************************************
/**
 * Creates a vector from data with a given underlying type
 * @param {Array | TypedArray} data - The input data
 * @param {Object} [type] - The underlying type. Defaults to the current default type
 * @returns {Mat} A vector with a typed copy of the given data
 */
function vecFrom(data, type) {
    var n = data.length;
    type = optional(type, currentDefaultType);
    return TypedMatFactory.new(type).copy(Mat.new(data, n, 1));
}
//*****************************************************
/**
 * Creates a new uninitialized vector with the currently set default type
 * 
 * @param {number} n - The number of elements
 * @returns {Mat} An unitialized vector
 */
function vec(n, type) {
    type = optional(type, currentDefaultType);
    return TypedMatFactory.new(type).uninitialized(n, 1);
}
//*****************************************************
/**
 * Reduces a matrix along its columns.
 * 
 * A MxN Matrix will be reduced to a 1xM row vector.
 * The callback function is used to compute the reduced result of a column that is written in the corresponding output slot
 * 
 * @param {Abstractmat} m - The input matrix
 * @param {function} f - Function (col,j) => *. Will be called for each column. The return value will be written in the output at position j
 * @param {AbstractMat} [out] - The output matrix. If not specified, a new matrix will be created 
 * @returns {AbstractMat} The output
 */
function colreduce(m, f, out) {
    out = out !== undefined ? out : similar(m, 1, m.cols());

    for (let j = 0; j < m.cols(); j++) {
        out.set(f(col(m, j), j), 0, j);
    }

    return out;
}
//*****************************************************
/**
 * Calls a function for each column.
 * 
 * This allows for changing the columns.
 * 
 * @param {Abstractmat} m - The input matrix
 * @param {function} f - Function (col,j) => void. Will be called for each column.
 * @returns {AbstractMat} m
 */
function colwise(m, f) {

    for (let j = 0; j < m.cols(); j++) {
        f(col(m, j), j);
    }

    return m;
}
//*****************************************************
/**
 * Computes the norm of the given matrix.
 * 
 * This computes the 2 p-norm. For vectors, this corresponds to the euclidean norm and for general matrices the frobenius norm
 * @see norm2
 * 
 * @param {AbstractMat} m - The input matrix
 * @returns {number} The norm
 */
function norm(m) {
    return norm2(m);
}
//*****************************************************
/**
 * Computes the 2 p-norm of the given matrix.
 * 
 * For vectors, this corresponds to the euclidean norm and for general matrices the frobenius norm
 * 
 * @param {AbstractMat} m - The input matrix
 * @returns {number} The 2 p-norm
 */
function norm2(m) {
    return Math.sqrt(norm2Squared(m));
}
//*****************************************************
/**
 * Computes the squared 2 p-norm of the given matrix.
 * 
 * This equals the square of the 2 p-norm
 * 
 * @param {AbstractMat} m - The input matrix
 * @returns {number} The squared 2 p-norm
 */
function norm2Squared(m) {
    return reduce(m, (acc, v) => acc + v * v, 0.0);
}
//*****************************************************
/**
 * Alias for norm2
 * @see norm2
 */
var normFrobenius = norm2;
//*****************************************************
/**
 * Computes the 1 p-norm of the given matrix.
 * 
 * This is the taxicab norm
 * 
 * @param {AbstractMat} m - The input matrix
 * @returns {number} The 1 p-norm
 */
function norm1(m) {
    return reduce(m, (acc, v) => acc + Math.abs(v), 0.0);
}
//*****************************************************
/**
 * Computes the Inf p-norm of the given matrix.
 * 
 * This is the maximum norm
 * 
 * @param {AbstractMat} m - The input matrix
 * @returns {number} The 1 p-norm
 */
function normInf(m) {
    return reduce(m, (acc, v) => Math.max(acc, Math.abs(v)), 0.0);
}
//*****************************************************
/**
 * Normalizes the given matrix with respect to the 2 p-norm
 * 
 * This computes for each element m_{i,j} = m_{i,j}/norm(m)
 * @param {AbstractMat} m - The input matrix
 * @param {AbstractMat} [out] - The output matrix. If not specified, a new matrix will be created
 * @returns {AbstractMat} The output
 */
function normalize(m, out) {
    var n = norm(m);

    out = out !== undefined ? out : copy(m);

    return scale(m, 1.0 / n, out);
}
//*****************************************************
/**
 * Reduces a matrix along its rows.
 * 
 * A MxN Matrix will be reduced to a Mx1 vector.
 * The callback function is used to compute the reduced result of a row that is written in the corresponding output slot
 * 
 * @param {Abstractmat} m - The input matrix
 * @param {function} f - Function (row,i) => *. Will be called for each row. The return value will be written in the output at position i
 * @param {AbstractMat} [out] - The output matrix. If not specified, a new matrix will be created 
 * @returns {AbstractMat} The output
 */
function rowreduce(m, f, out) {
    out = out !== undefined ? out : similar(m, m.rows(), 1);

    for (let i = 0; i < m.rows(); i++) {
        out.set(f(row(m, i), i), i, 0);
    }

    return out;
}
//*****************************************************
/**
 * Calls a function for each row.
 * 
 * This allows for changing the rows.
 * 
 * @param {Abstractmat} m - The input matrix
 * @param {function} f - Function (row,i) => void. Will be called for each row.
 * @returns {AbstractMat} m
 */
function rowwise(m, f) {

    for (let i = 0; i < m.rows(); i++) {
        f(row(m, i), i);
    }

    return m;
}
//*****************************************************
/**
 * Copies the contents for b into a.
 * 
 * Both matrices need to be of the same size
 * 
 * @param {AbstractMat} a - The matrix to write into
 * @param {AbstractMat} b - The matrix to be read from
 * @returns {AbstractMat} a
 */
function insert(a, b) {
    if (a.rows() !== b.rows() || a.cols() !== b.cols()) {
        throw "Insertion failed: Source and target have different sizes";
    }

    var c = a.cols();
    var r = b.rows();

    for (let j = 0; j < c; j++) {
        for (let i = 0; i < r; i++) {
            a.set(b.at(i, j), i, j);
        }
    }

    return a;

}
//*****************************************************
/**
 * Swaps two rows of a matrix
 * 
 * @param {AbstractMat} a - The input matrix
 * @param {number} row0 - The first row
 * @param {number} row1 - The second row
 * @returns {AbstractMat} a
 */
function swapRow(a, row0, row1) {
    var r0 = copy(row(a, row0));
    insert(row(a, row0), row(a, row1));
    insert(row(a, row1), r0);

    return a;
}
//*****************************************************
/**
 * Swaps two columns of a matrix
 * 
 * @param {AbstractMat} a - The input matrix
 * @param {number} col0 - The first column
 * @param {number} col1 - The second column
 * @returns {AbstractMat} a
 */
function swapCol(a, col0, col1) {
    var c0 = copy(col(a, col0));
    insert(col(a, col0), col(a, col1));
    insert(col(a, col1), c0);

    return a;
}
//*****************************************************
/**
 * Computes the determinant of a square matrix
 * 
 * @param {AbstractMat} m - The input matrix
 * @returns {number} The determinant
 */
function det(m) {
    if (m.rows() !== m.cols()) {
        throw "Determinant only defined for square sources";
    }

    // check if m has a dedicated det function
    if (m.det && m.det instanceof Function) {
        return m.det();
    }

    var n = m.rows();

    if (n === 1) {
        return m.at(0, 0);
    }
    else if (n === 2) {
        return m.at(0, 0) * m.at(1, 1) - m.at(0, 1) * m.at(1, 0);
    } else if (n === 3) {
        var a = m.at(0, 0);
        var b = m.at(0, 1);
        var c = m.at(0, 2);

        var d = m.at(1, 0);
        var e = m.at(1, 1);
        var f = m.at(1, 2);

        var g = m.at(2, 0);
        var h = m.at(2, 1);
        var i = m.at(2, 2);


        return a * e * i + b * f * g + c * d * h - c * e * g - b * d * i - a * f * h;

    }
    else {

        // TODO use more efficient formulation

        var lu = computePLUD(m);

        // singular matrix
        if (!lu) {
            return 0.0;
        }
        // lower diagonal is 1s -> determinant is 1
        // total determinant is therefore just product of u diagonal
        var s = (lu.numSwaps % 2 === 0 ? 1 : -1) * reduce(diag(lu.U), (acc, v) => acc * v);
        // every invertable matrix can be PLU decomposed
        // if there are any NaNs, the matrix was not invertible
        if (isNaN(s)) {
            s = 0.0;
        }
        return s;
    }
}
//*****************************************************
/**
 * Solves mx = b for x, where m is a lower triagonal matrix.
 * 
 * @param {AbstractMat} m - The lower triagonal matrix
 * @param {AbstractMat} b - The matrix to solve for
 * @param {AbstractMat} [out] - The output matrix. If not specified, a new one will be created
 * @returns {AbstractMat} The result
 */
function solveLowerTriagonal(m, b, out) {
    var rows = m.rows();
    var cols = m.cols();

    var brows = b.rows();
    var bcols = b.cols();


    if (rows !== cols) {
        throw "Lower triagonal matrix not square";
    }

    if (rows !== brows) {
        throw "b does not have correct size";
    }

    out = out !== undefined ? out : similar(b);

    for (let bc = 0; bc < bcols; bc++) {
        for (let i = 0; i < rows; i++) {
            var sumOffDiag = 0;
            for (let j = 0; j < i; j++) {
                sumOffDiag += m.at(i, j) * out.at(j, bc);
            }

            out.set((b.at(i, bc) - sumOffDiag) / m.at(i, i), i, bc);
        }
    }

    return out;
}
//*****************************************************
/**
 * Solves mx = b for x, where m is a upper triagonal matrix.
 * 
 * @param {AbstractMat} m - The upper triagonal matrix
 * @param {AbstractMat} b - The matrix to solve for
 * @param {AbstractMat} [out] - The output matrix. If not specified, a new one will be created
 * @returns {AbstractMat} The result
 */
function solveUpperTriagonal(m, b, out) {
    var rows = m.rows();
    var cols = m.cols();

    var brows = b.rows();
    var bcols = b.cols();


    if (rows !== cols) {
        throw "Lower triagonal matrix not square";
    }

    if (rows !== brows) {
        throw "b does not have correct size";
    }

    out = out !== undefined ? out : similar(b);

    for (let bc = 0; bc < bcols; bc++) {
        for (let i = rows - 1; i >= 0; i--) {
            var sumOffDiag = 0;
            for (let j = rows - 1; j > i; j--) {
                sumOffDiag += m.at(i, j) * out.at(j, bc);
            }

            out.set((b.at(i, bc) - sumOffDiag) / m.at(i, i), i, bc);
        }
    }

    return out;
}
//*****************************************************
/**
 * Checks, if the given matrix is square
 * 
 * @param {AbstractMat} a - The input matrix
 * @returns {boolean} True, if the matrix is square, false otherwise
 */
function isSquare(a) {
    return a.rows() === a.cols();
}
//*****************************************************
/**
 * Solves the linear system ax = b for x
 * 
 * This will not check, if the matrix is invertable and assume it is in the square case.
 * 
 * Solving for square matrices is accomplished via a PLU decomposition. If that fails due to singularities,
 * the result will be computed with an SVD.
 * Rectangular matrices are solved using the SVD, thus not solving exactly, but minimizing |Ax - b|_2
 * 
 * @param {AbstractMat} a - The linear system
 * @param {AbstractMat} b - The matrix to solve for
 * @param {AbstractMat} [out] - The output matrix. If not specified, a new matrix will be created
 * @returns {AbstractMat} The output
 */
function solve(a, b, out) {

    // check if a has a dedicated solve function
    if (a.solve && a.solve instanceof Function) {
        return a.solve(b, out);
    }
    if (isSquare(a)) {
        // try plu
        var plu = computePLUD(a);

        // singular matrix solve with svd instead
        if (!plu) {
            var svd = computeSVD(a);
            return svd.solve(b, out);
        }
        return plu.solve(b, out);
    }

    // for non square matrices -> use svd

    var svd = computeSVD(a);
    return svd.solve(b, out);
}
//*****************************************************
/**
 * Solves the linear system ax = b for x using an already computed PLU decomposition.
 * 
 * This will not check, if the matrix is invertable and assume it is
 * 
 * @param {PLUD} plu - The PLU decomposition
 * @param {AbstractMat} b - The matrix to solve for
 * @param {AbstractMat} [out] - The output matrix. If not specified, a new matrix will be created
 * @returns {AbstractMat} The output
 */
function solvePLU(plu, b, out) {
    var pb = mult(plu.P, b);

    var y = solveLowerTriagonal(plu.L, pb);
    var x = solveUpperTriagonal(plu.U, y, out);

    return x;

}
//*****************************************************
/**
 * Solves the linear system Ax = b for x using an already computed SVD decomposition.
 * 
 * This will find x, such that |Ax - b|_2 is minimized
 * 
 * @param {SVD} svd - The SVD
 * @param {AbstractMat} b - The matrix to solve for
 * @param {AbstractMat} [out] - The output matrix. If not specified, a new matrix will be created
 * @returns {AbstractMat} The output
 */
function solveSVD(svd, b, out) {
    var U = svd.U;
    var V = svd.V;
    var S = svd.S;

    var M = U.rows();
    var N = V.rows();

    var P = b.cols();

    let r = 0;

    out = out !== undefined ? out : similar(b, N, P);
    setZero(out);

    for (let j = 0; j < P; j++) {
        var bj = col(b, j);
        var zj = col(out, j);

        for (r = 0; r < S.rows(); r++) {
            var sigma = S.at(r);
            // TODO relative epsilon?
            if (sigma <= 1E-7) {
                break;
            }

            var zr = dot(col(U, r), bj) / sigma;
            add(zj, scale(col(V, r), zr), zj);
        }

    }


    return out;
}
//*****************************************************
/**
 * Computes the inverse of a square matrix
 * 
 * This will not check, whether a matrix is invertable
 * 
 * @param {AbstractMat} m - The input matrix
 * @returns {AbstractMat} The inverse
 */
function inv(m) {
    if (m.rows() !== m.cols()) {
        throw "Inverse only defined for square sources";
    }

    // check if m has a dedicated inv function
    if (m.inv && m.inv instanceof Function) {
        return m.inv();
    }

    var n = m.rows();

    if (n === 1) {
        return similar(m).set(1 / m.at(0, 0), 0, 0);
    } else if (n === 2) {
        let a = m.at(0, 0);
        let b = m.at(0, 1);
        let c = m.at(1, 0);
        let d = m.at(1, 1);

        let f = 1 / (a * d - b * c);

        let result = similar(m);
        result.set(f * d, 0, 0);
        result.set(f * a, 1, 1);
        result.set(-f * c, 1, 0);
        result.set(-f * b, 0, 1);
        return result;
    }
    else if (n === 3) {
        let a = m.at(0, 0);
        let b = m.at(0, 1);
        let c = m.at(0, 2);

        let d = m.at(1, 0);
        let e = m.at(1, 1);
        let f = m.at(1, 2);

        let g = m.at(2, 0);
        let h = m.at(2, 1);
        let i = m.at(2, 2);

        let det = a * e * i + b * f * g + c * d * h - c * e * g - b * d * i - a * f * h;
        let factor = 1 / det;
        let A = (e * i - f * h);
        let B = -(d * i - f * g);
        let C = (d * h - e * g);
        let D = -(b * i - c * h);
        let E = (a * i - c * g);
        let F = -(a * h - b * g);
        let G = (b * f - c * e);
        let H = -(a * f - c * d);
        let I = (a * e - b * d);

        let result = similar(m);

        result.set(factor * A, 0, 0);
        result.set(factor * B, 1, 0);
        result.set(factor * C, 2, 0);

        result.set(factor * D, 0, 1);
        result.set(factor * E, 1, 1);
        result.set(factor * F, 2, 1);

        result.set(factor * G, 0, 2);
        result.set(factor * H, 1, 2);
        result.set(factor * I, 2, 2);


        return result;
    }
    else {

        // Blockwise inversion
        // Old
        // let na = Math.max(n - 1, Math.floor(n / 2));

        // let A = block(m, 0, 0, na, na);
        // let B = block(m, 0, na, na, 1);
        // let C = block(m, na, 0, 1, na);
        // let D = block(m, na, na, 1, 1);



        // let Ai = inv(A);
        // let CAi = mult(C, Ai);
        // let DCA1B = sub(D, mult(CAi, B));
        // let DCA1Bi = inv(DCA1B);

        // var tl = mult(Ai, add(DiagonalF.id(na), mult(B, mult(DCA1Bi, CAi))));
        // var tr = mult(Ai, mult(B, DCA1Bi));
        // neg(tr, tr);

        // var bl = mult(DCA1Bi, CAi);
        // neg(bl, bl);

        // var br = DCA1Bi;

        // var result = similar(m);

        // insert(block(result, 0, 0, na, na), tl);
        // insert(block(result, 0, na, na, 1), tr);
        // insert(block(result, na, 0, 1, na), bl);
        // insert(block(result, na, na, 1, 1), br);

        // Invert using PLU decomposition
        let id = TypedDiagonalFactory.newFromMat(m).id(n);
        let result = solve(m, id);

        return result;

    }

}
//*****************************************************
/**
 * LU decomposition with partial pivoting
 */
class PLUD {
    /**
     * 
     * @param {AbstractMat} P - Permutation matrix
     * @param {AbstractMat} L - Lower unit triangular matrix
     * @param {AbstractMat} U - Upper triangular matrix
     * @param {number} numSwaps - The number of swaps in the permutation matrix
     */
    constructor(P, L, U, numSwaps) {
        this.P = P;
        this.L = L;
        this.U = U;
        this.numSwaps = numSwaps;
    }
    /**
    * Creates a new PLUD object
    * 
    * @param {AbstractMat} P - Permutation matrix
    * @param {AbstractMat} L - Lower unit triangular matrix
    * @param {AbstractMat} U - Upper triangular matrix
    * @param {number} numSwaps - The number of swaps in the permutation matrix
    * 
    * @returns {PLUD} - A new PLUD object
    * 
    * @see constructor
    */
    static new(P, L, U, numSwaps) {
        return new PLUD(P, L, U, numSwaps);
    }

    /**
     * Solves Ax=b for x, where A is the matrix represented as this decomposition
     * 
     * @param {AbstractMat} b - The matrix to solve for
     * @param {AbstractMat} [out] - Optional output
     */
    solve(b, out) {
        return solvePLU(this, b, out);
    }

    toMat(out) {
        var M = this.P.rows();
        var N = this.U.cols();

        out = out !== undefined ? out : similar(this.U, M, N);

        if (out.rows() !== M || out.cols() !== N) {
            throw "Output has wrong size";
        }

        return mult(transpose(this.P), mult(this.L, this.U), out);
    }


}
//*****************************************************
/**
 * Computes an LU decomposition with partial pivoting.
 * 
 * The decomposition of a Matrix A is defined by: P * A = L * U
 * 
 * P is a permutation matrix
 * L is an upper triangular matrix
 * U is a lower triangular matrix
 * 
 * A can be reconstructed by A = P^T * L * U
 * 
 * @param {AbstractMat} a Input matrix
 * @param {AbstractMat} [out] Output whre LU will be compactly stored. If not specified a new matrix will be created. 
 * @returns {PLUD} The LU decomposition with partial pivoting
 */
function computePLUD(a, out) {

    var n = a.rows();
    var r = a.cols();
    var rnmin = Math.min(a.rows(), a.cols());

    var permutes = TypedMatFactory.new(Int32Array).uninitialized(n, 1);
    map(permutes, (v, i) => i, permutes);

    if (out === undefined) {
        out = similar(a);

    }
    insert(out, a);
    var numSwaps = 0;

    for (let k = 0; k < rnmin; k++) {

        // find largest value in colum for partial pivot

        var maxIndex = k;
        var maxEl = Math.abs(out.at(maxIndex, k));
        for (let l = k + 1; l < n; l++) {
            if (Math.abs(out.at(l, k)) > maxEl) {
                maxIndex = l;
                maxEl = Math.abs(out.at(l, k));
            }
        }

        // swap row k with maxIndex
        if (maxIndex != k) {
            numSwaps++;
            swapRow(permutes, maxIndex, k);
            swapRow(out, maxIndex, k);
        }

        // Algorithm from "Matrix computations"

        var outkk = out.at(k, k);

        // singularity detected
        if (Math.abs(outkk) < 1E-7) {
            return null;
        }
        var subcol = subvec(col(out, k), k + 1);
        scale(subcol, 1.0 / outkk, subcol);

        // update lower block
        // case n > r
        if (k < r) {
            var rowRho = subrowvec(row(out, k), k + 1);

            for (let rho = k + 1; rho < n; rho++) {
                var subrow = subrowvec(row(out, rho), k + 1);
                var arhok = out.at(rho, k);
                sub(subrow, scale(rowRho, arhok), subrow);

            }
        }
    }

    var blockL = block(out, 0, 0, a.rows(), rnmin);
    var blockU = block(out, 0, 0, rnmin, a.cols());
    var L = TriangularView.new(blockL, TriangularMode.UNIT_LOWER);
    var U = TriangularView.new(blockU, TriangularMode.UPPER);
    var P = RowPermutation.new(toArray(permutes), permutes.rows(), permutes.rows(), a.type());
    return PLUD.new(P, L, U, numSwaps);
}
//*****************************************************
/**
 * Scales an input matrix by some value
 * 
 * @param {AbstractMat} a - The input matrix
 * @param {number} v - The value to scale by
 * @param {AbstractMat} [out] - Output of the operation. If not specified, a new matrix will be created
 * @returns {AbstractMat} The result of a*v
 */
function scale(a, v, out) {
    out = out !== undefined ? out : similar(a);

    for (let j = 0; j < a.cols(); j++) {
        for (let i = 0; i < a.rows(); i++) {
            out.set(a.at(i, j) * v, i, j);
        }
    }

    return out;
}
//*****************************************************
// Over-/underflow safe (a^2 + b^2)^(1/2)
function hypot(a, b) {
    var absa = Math.abs(a);
    var absb = Math.abs(b);

    var sqr = function (x) {
        return x * x;
    };
    if (absa > absb) {

        return absa * Math.sqrt(1.0 + sqr(absb / absa));
    }

    return absb === 0.0 ? 0.0 : absb * Math.sqrt(1.0 + sqr(absa / absb));
}
//*****************************************************
function householderVector(x, out) {
    if (!isVec(x)) {
        throw "Householder transform needs to operate on vector";
    }

    let v = out !== undefined ? out : copy(x);
    var m = x.rows();

    // compute squared length of subvector starting at i = 1

    var sigma = 0.0;
    for (let i = 1; i < m; i++) {
        let vi = x.at(i);
        sigma += vi * vi;
    }

    var x0 = x.at(0);
    v.set(1.0, 0);


    if (sigma === 0.0) {
        return { beta: 0.0, v: v };
    }

    var my = Math.sqrt(x0 * x0 + sigma);

    if (x0 <= 0.0) {
        v.set(x0 - my, 0);
    }
    else {
        v.set(-sigma / (x0 + my), 0);
    }


    var v0 = v.at(0);
    var v02 = v0 * v0;
    var beta = 2.0 * v02 / (sigma + v02);

    scale(v, 1.0 / v0, v);
    return { beta: beta, v: v };

}
//*****************************************************
/**
 * Computes sum_{j}^{cols}sum_{i}^{rows} a[i,j]*b[i,j]
 * 
 * For vectors this is the standard dot product
 * 
 * @param {AbstractMat} a - The first matrix
 * @param {AbstractMat} b - The second matrix
 * @returns The sum of all multiplied corresponding elements in a and b
 */
function dot(a, b) {
    var r = a.rows();
    var c = a.cols();

    if (r !== b.rows() || c !== b.cols()) {
        throw "Inputs must match in size for dot product";
    }

    var s = 0.0;

    for (let j = 0; j < c; j++) {
        for (let i = 0; i < r; i++) {
            s += a.at(i, j) * b.at(i, j);
        }
    }

    return s;
}
//*****************************************************
/**
 * Constructs a subvector view for a given vector
 * 
 * @param {AbstractMat} v - The base vector
 * @param {number} start - The start index
 * @param {number} [rows] - The number of rows. When not specified, the rows will be the remaining one in the base vector
 * @returns {BlockView} A subvector view of the base vector
 */
function subvec(v, start, rows) {
    if (!isVec(v)) {
        throw "Input for subvec needs to be a vector";
    }
    rows = rows !== undefined ? rows : v.rows() - start;
    return block(v, start, 0, rows, 1);
}
//*****************************************************
/**
 * Constructs a subvector view for a given row vector
 * 
 * @param {AbstractMat} v - The base row vector
 * @param {number} start - The start index
 * @param {number} [cols] - The number of columns. When not specified, the columns will be the remaining one in the base vector
 * @returns {BlockView} A subvector view of the base vector
 */
function subrowvec(v, start, cols) {
    if (!isRowVec(v)) {
        throw "Input for subvec needs to be a row vector";
    }
    cols = cols !== undefined ? cols : v.cols() - start;
    return block(v, 0, start, 1, cols);
}
//*****************************************************
// Adapted from Matrix Multiplication 5.1
function applyHouseholderLeft(beta, v, a, out) {
    out = out !== undefined ? out : copy(a);
    if (beta === 0.0) {
        return out;
    }

    var r = a.rows();
    var c = a.cols();

    // apply per col
    for (let j = 0; j < c; j++) {
        // compute v^T * A[:,j]

        var aj = col(out, j);
        let wj = aj.at(0);

        wj += dot(subvec(v, 1), subvec(aj, 1));

        var bwj = beta * wj;

        // apply to column: (aj)_i' = (aj)_i - (w_j*beta)*v_i
        {
            let a0j = aj.at(0);
            aj.set(a0j - beta * wj, 0);
        }
        for (let i = 1; i < r; i++) {
            aj.set(aj.at(i) - bwj * v.at(i), i);
        }

    }

    return out;

}
//*****************************************************
// Adapted from Matrix Multiplication 5.1
function applyHouseholderRight(beta, v, a, out) {
    out = out !== undefined ? out : copy(a);
    if (beta === 0.0) {
        return out;
    }

    var r = a.rows();
    var c = a.cols();

    // apply per row
    for (let i = 0; i < r; i++) {
        var ai = transpose(row(out, i));

        let wi = ai.at(0);

        wi += dot(subvec(v, 1), subvec(ai, 1));

        var bwi = beta * wi;

        // apply to entries in row
        {
            let a0i = ai.at(0);
            ai.set(a0i - beta * wi, 0);
        }
        for (let j = 1; j < c; j++) {
            ai.set(ai.at(j) - bwi * v.at(j), j);
        }
    }

    return out;
}
//*****************************************************
/**
 * Computes an upper bidiagonal decomposition.
 * 
 * This decomposition algorithm is only defined for rows >= cols.
 * 
 * The decomposition is defined by: B = U^T * A * V. A can then be computed from B with A = U * B * V^T
 * B is an upper diagonal matrix. U and V are products of Householder matrices, with U = U_0 *... * U_n and V = V_0 * ... * V_{n-3}
 * 
 * The result is stored in a matrix the same size as the input. 
 * The essential parts of the U_i are stored in the columns below the diagonal, where only zeroes would be. 
 * The essential parts of V_i  are analogously stored in the rows above the bidiagonal. The matrices themselves can be retrieved via the unpackUBV function 
 * 
 * @param {AbstractMat} a - The input matrix with a.rows() >= a.cols()
 * @param {AbstractMat} [out] - Output matrix. If none is specified, a new matrix is created. a itself can be used
 * @returns {AbstractMat} The packed result of the bidiagonalization
 */
function computeUBVD(a, out) {
    var M = a.rows();
    var N = a.cols();

    if (M < N) {
        throw "Biadiagonal only implemented for M >= N";
    }
    if (out !== undefined) {
        insert(out, a);
        a = out;
    }
    else {
        a = copy(a);
    }

    for (let j = 0; j < N; j++) {
        let cj = subvec(col(a, j), j);
        let house = householderVector(cj);

        var blockj = block(a, j, j, M - j, N - j);
        applyHouseholderLeft(house.beta, house.v, blockj, blockj);

        insert(subvec(cj, 1), subvec(house.v, 1, M - j - 1));

        if (j < N - 2) {
            let rj = subvec(transpose(row(a, j)), j + 1);
            let blockj = block(a, j, j + 1, M - j, N - (j + 1));
            let house = householderVector(rj);
            applyHouseholderRight(house.beta, house.v, blockj, blockj);

            insert(subvec(rj, 1), subvec(house.v, 1, N - j - 2));
        }
    }


    return a;
}
//*****************************************************
/**
 * Represents a UBV decomposition of some matrix A.
 * 
 * B is an upper bidiagonal matrix with B = U^T * A * V
 */
class UBVD {

    constructor(U, Vt, B) {
        this.U = U;
        this.Vt = Vt;
        this.B = B;
    }

    static new(U, Vt, B) {
        return new UBVD(U, Vt, B);
    }
}
//*****************************************************
/**
 * Reconstructs U,B and V^T matrices from a packed bidiagonalization decomposition.
 * 
 * @see decomposeBidiag
 * 
 * @param {AbstractMat} ubv The packed bidiagonal decomposition
 * @returns {UBVD} Object containg the U,B and V^T matrices  
 */
function unpackUBV(ubv) {

    let M = ubv.rows();
    let N = ubv.cols();
    var fac = TypedMatFactory.newFromMat(ubv);
    var B = fac.zeros(ubv.rows(), ubv.cols());

    insert(diag(B), diag(ubv));
    // superdiag
    for (let i = 0; i < N - 1; i++) {
        B.set(ubv.at(i, i + 1), i, i + 1);
    }

    var U = fac.zeros(M, M);
    fill(subvec(diag(U), 0, N), 1.0);

    let house = fac.zeros(M, 1);


    for (let j = N - 1; j >= 0; j--) {
        // extract essential
        let hv = subvec(house, j);
        hv.set(1.0, 0);

        let cj = subvec(col(ubv, j), j);
        let cjs = subvec(cj, 1);
        insert(subvec(hv, 1), cjs);

        let sigma = dot(cjs, cjs);
        let betaj = sigma === 0.0 ? 0.0 : 2.0 / (1.0 + sigma);
        let blockj = block(U, j, j, M - j, N - j);

        applyHouseholderLeft(betaj, hv, blockj, blockj);
    }

    var Vt = fac.id(N);
    house = fac.zeros(N - 1, 1);
    for (let i = N - 3; i >= 0; i--) {
        let hv = subvec(house, i);
        hv.set(1.0, 0);

        let ri = subvec(transpose(row(ubv, i)), i + 1);
        let ris = subvec(ri, 1);
        insert(subvec(hv, 1), ris);

        let sigma = dot(ris, ris);
        let betai = sigma === 0.0 ? 0.0 : 2.0 / (1.0 + sigma);
        let blocki = block(Vt, i + 1, i + 1, N - (i + 1), N - (i + 1));
        applyHouseholderRight(betai, hv, blocki, blocki);
    }

    return UBVD.new(U, Vt, B);
}
//*****************************************************
// Adapted from GSL library
function createGivens(a, b) {
    var c, s;
    if (b === 0.0) {
        c = 1.0;
        s = 0.0;
    }
    else if (Math.abs(b) > Math.abs(a)) {
        let t = -a / b;
        let s1 = 1.0 / Math.sqrt(1 + t * t);
        s = s1;
        c = s1 * t;
    }
    else {
        let t = -b / a;
        let c1 = 1.0 / Math.sqrt(1 + t * t);
        c = c1;
        s = c1 * t;
    }

    return { c: c, s: s };
}
//*****************************************************
// Adapted from GSL library
function trailingEigenvalue(d, f) {
    var n = d.rows();

    var da = d.at(n - 2);
    var db = d.at(n - 1);
    var fa = (n > 2) ? f.at(n - 3) : 0.0;
    var fb = f.at(n - 2);

    var ta = da * da + fa * fa;
    var tb = db * db + fb * fb;
    var tab = da * fb;

    var dt = (ta - tb) / 2.0;

    var mu;

    if (dt >= 0) {
        mu = tb - (tab * tab) / (dt + hypot(dt, tab));
    }
    else {
        mu = tb + (tab * tab) / ((-dt) + hypot(dt, tab));
    }

    return mu;
}
//*****************************************************
// Adapted from GSL library
function createSchur(d0, f0, d1) {
    var apq = 2.0 * d0 * f0;
    var c, s;
    if (d0 === 0.0 || f0 === 0.0) {
        c = 1.0;
        s = 0.0;
        return { c: c, s: s };
    }

    // TODO implement rescaling

    if (apq !== 0.0) {
        let t;
        let tau = (f0 * f0 + (d1 + d0) * (d1 - d0)) / apq;

        if (tau >= 0.0) {
            t = 1.0 / (tau + hypot(1.0, tau));
        }
        else {
            t = -1.0 / (-tau + hypot(1.0, tau));
        }

        c = 1.0 / hypot(1.0, t);
        s = t * c;
    }
    else {
        c = 1.0;
        s = 0.0;
    }

    return { c: c, s: s };
}
//*****************************************************
// Adapted from GSL library
function svd2(d, f, U, V) {
    var i;
    var c, s, a11, a12, a21, a22;

    var M = U.rows();
    var N = V.rows();

    var d0 = d.at(0);
    var f0 = f.at(0);

    var d1 = d.at(1);

    if (d0 === 0.0) {
        /* Eliminate off-diagonal element in [0,f0;0,d1] to make [d,0;0,0] */

        let giv = createGivens(f0, d1);
        c = giv.c;
        s = giv.s;

        /* compute B <= G^T B X,  where X = [0,1;1,0] */

        d.set(c * f0 - s * d1, 0);
        f.set(s * f0 + c * d1, 0);
        d.set(0.0, 1);

        /* Compute U <= U G */

        for (i = 0; i < M; i++) {
            let Uip = U.at(i, 0);
            let Uiq = U.at(i, 1);

            U.set(c * Uip - s * Uiq, i, 0);
            U.set(s * Uip + c * Uiq, i, 1);
        }

        /* Compute V <= V X */

        swapCol(V, 0, 1);

        return;
    }
    else if (d1 === 0.0) {
        /* Eliminate off-diagonal element in [d0,f0;0,0] */

        let giv = createGivens(d0, f0);
        c = giv.c;
        s = giv.s;
        /* compute B <= B G */

        d.set(d0 * c - f0 * s, 0);
        f.set(0.0, 0);

        /* Compute V <= V G */

        for (i = 0; i < N; i++) {
            let Vip = V.at(i, 0);
            let Viq = V.at(i, 1);

            V.set(c * Vip - s * Viq, i, 0);
            V.set(s * Vip + c * Viq, i, 1);
        }

        return;
    }
    else {
        /* Make columns orthogonal, A = [d0, f0; 0, d1] * G */

        var sh = createSchur(d0, f0, d1);
        c = sh.c;
        s = sh.s;
        /* compute B <= B G */

        a11 = c * d0 - s * f0;
        a21 = - s * d1;

        a12 = s * d0 + c * f0;
        a22 = c * d1;

        /* Compute V <= V G */

        for (i = 0; i < N; i++) {
            let Vip = V.at(i, 0);
            let Viq = V.at(i, 1);

            V.set(c * Vip - s * Viq, i, 0);
            V.set(s * Vip + c * Viq, i, 1);
        }

        /* Eliminate off-diagonal elements, bring column with largest
           norm to first column */

        if (hypot(a11, a21) < hypot(a12, a22)) {
            let t1, t2;

            /* B <= B X */

            t1 = a11; a11 = a12; a12 = t1;
            t2 = a21; a21 = a22; a22 = t2;

            /* V <= V X */
            swapCol(V, 0, 1);
        }

        let giv = createGivens(a11, a21);
        c = giv.c;
        s = giv.s;

        /* compute B <= G^T B */

        d.set(c * a11 - s * a21, 0);
        f.set(c * a12 - s * a22, 0);
        d.set(s * a12 + c * a22, 1);

        /* Compute U <= U G */

        for (i = 0; i < M; i++) {
            let Uip = U.at(i, 0);
            let Uiq = U.at(i, 1);

            U.set(c * Uip - s * Uiq, i, 0);
            U.set(s * Uip + c * Uiq, i, 1);
        }

        return;
    }
}
//*****************************************************
// Adapted from GSL library
function chaseOutIntermediateZero(d, f, U, k0) {

    var M = U.rows();

    var n = d.rows();
    var c, s;
    var x, y;
    var k;

    x = f.at(k0);
    y = d.at(k0 + 1);

    for (k = k0; k < n - 1; k++) {
        var giv = createGivens(y, -x);
        c = giv.c;
        s = giv.s;
        /* Compute U <= U G */


        {
            let i;

            for (i = 0; i < M; i++) {
                let Uip = U.at(i, k0);
                let Uiq = U.at(i, k + 1);

                U.set(c * Uip - s * Uiq, i, k0);
                U.set(s * Uip + c * Uiq, i, k + 1);
            }
        }

        /* compute B <= G^T B */

        d.set(s * x + c * y, k + 1);

        if (k === k0) {
            f.set(c * x - s * y, k);
        }

        if (k < n - 2) {
            let z = f.at(k + 1);
            f.set(c * z, k + 1);

            x = -s * z;
            y = d.at(k + 2);
        }
    }
}
//*****************************************************
// Adapted from GSL library
function chaseOutTrailingZero(d, f, V) {

    var N = V.rows();
    var n = d.rows();
    var c, s;
    var x, y;
    var k;

    x = d.at(n - 2);
    y = f.at(n - 2);

    for (k = n - 1; k-- > 0;) {
        var giv = createGivens(x, y);
        c = giv.c;
        s = giv.s;
        /* Compute V <= V G where G = [c, s ; -s, c] */

        {
            let i;

            for (i = 0; i < N; i++) {
                let Vip = V.at(i, k);
                let Viq = V.at(i, n - 1);

                V.set(c * Vip - s * Viq, i, k);
                V.set(s * Vip + c * Viq, i, n - 1);
            }
        }

        /* compute B <= B G */

        d.set(c * x - s * y, k);

        if (k === n - 2) {

            f.set(s * x + c * y, k);

        }
        if (k > 0) {
            let z = f.at(k - 1);
            f.set(c * z, k - 1);

            x = d.at(k - 1);
            y = s * z;
        }
    }
}
//*****************************************************
// Adapted from GSL library
function qrstep(d, f, U, V) {

    var M = U.rows();
    var N = V.rows();
    var n = d.rows();
    var y, z;
    var ak, bk, zk, ap, bp, aq;
    var i, k;

    if (n === 1)
        return;  /* shouldn't happen */

    /* Compute 2x2 svd directly */

    if (n === 2) {
        svd2(d, f, U, V);
        return;
    }

    /* Chase out any zeroes on the diagonal */

    for (i = 0; i < n - 1; i++) {
        let d_i = d.at(i);

        if (d_i === 0.0) {
            chaseOutIntermediateZero(d, f, U, i);
            return;
        }
    }

    /* Chase out any zero at the end of the diagonal */

    {
        let d_nm1 = d.at(n - 1);

        if (d_nm1 === 0.0) {
            chaseOutTrailingZero(d, f, V);
            return;
        }
    }


    /* Apply QR reduction steps to the diagonal and offdiagonal */

    {
        let d0 = d.at(0);
        let f0 = f.at(0);

        let d1 = d.at(1);

        {
            let mu = trailingEigenvalue(d, f);

            y = d0 * d0 - mu;
            z = d0 * f0;
        }

        /* Set up the recurrence for Givens rotations on a bidiagonal matrix */

        ak = 0;
        bk = 0;

        ap = d0;
        bp = f0;

        aq = d1;
    }

    for (k = 0; k < n - 1; k++) {
        let c, s;
        var giv = createGivens(y, z);
        c = giv.c;
        s = giv.s;

        /* Compute V <= V G */


        for (i = 0; i < N; i++) {
            let Vip = V.at(i, k);
            let Viq = V.at(i, k + 1);

            V.set(c * Vip - s * Viq, i, k);
            V.set(s * Vip + c * Viq, i, k + 1);
        }


        /* compute B <= B G */

        {
            let bk1 = c * bk - s * z;

            let ap1 = c * ap - s * bp;
            let bp1 = s * ap + c * bp;
            let zp1 = -s * aq;

            let aq1 = c * aq;

            if (k > 0) {
                f.set(bk1, k - 1);
            }

            ak = ap1;
            bk = bp1;
            zk = zp1;

            ap = aq1;

            if (k < n - 2) {
                bp = f.at(k + 1);
            }
            else {
                bp = 0.0;
            }

            y = ak;
            z = zk;
        }

        let giv2 = createGivens(y, z);
        c = giv2.c;
        s = giv2.s;
        /* Compute U <= U G */


        for (i = 0; i < M; i++) {
            let Uip = U.at(i, k);
            let Uiq = U.at(i, k + 1);

            U.set(c * Uip - s * Uiq, i, k);
            U.set(s * Uip + c * Uiq, i, k + 1);
        }

        /* compute B <= G^T B */

        {
            let ak1 = c * ak - s * zk;
            let bk1 = c * bk - s * ap;
            let zk1 = -s * bp;

            let ap1 = s * bk + c * ap;
            let bp1 = c * bp;

            d.set(ak1, k);

            ak = ak1;
            bk = bk1;
            zk = zk1;

            ap = ap1;
            bp = bp1;

            if (k < n - 2) {
                aq = d.at(k + 2);
            }
            else {
                aq = 0.0;
            }

            y = bk;
            z = zk;
        }
    }

    f.set(bk, n - 2);
    d.set(ap, n - 1);
}
//*****************************************************
// Adapted from GSL library
// Implements the unit round-off from "Matrix Computations" Algorithm 8.6.2
function chopSmallElements(d, f) {
    var N = d.rows();
    var d_i = d.at(0);

    for (let i = 0; i < N - 1; i++) {
        let f_i = f.at(i);
        let d_ip1 = d.at(i + 1);

        // TODO Use better epsilon
        if (Math.abs(f_i) < 1E-7 * (Math.abs(d_i) + Math.abs(d_ip1))) {
            f.set(0.0, i);
        }
        d_i = d_ip1;
    }

}
//*****************************************************
/**
 * A singular value decomposition U * S * V^T
 */
class SVD {

    /**
     * 
     * @param {AbstractMat} U - The U matrix
     * @param {AbstractMat} S - The singular values in a vector
     * @param {AbstractMat} V - The V matrix
     */
    constructor(U, S, V) {
        /** U */
        this.U = U;
        /** S matrix represented as a vector containing the singular values */
        this.S = S;
        /** V  */
        this.V = V;
        /** V^T */
        this.Vt = transpose(this.V);

    }
    /**
     * Creates a new SVD object
     * @param {AbstractMat} U - The U matrix
     * @param {AbstractMat} S - The singular values in a vector
     * @param {AbstractMat} V - The V matrix
     * @returns {SVD} A new SVD object
     * 
     * @see constructor
     */
    static new(U, S, V) {
        return new SVD(U, S, V);
    }

    /**
     * Finds the vector x that minimizes ||Ax - b||_2, where A is the matrix decomposed into this SVD.
     * @see solveSVD
     * @param {AbstractMat} b - The vector to solve for
     * @param {AbstractMat} [out] - Optional output vector
     */
    solve(b, out) {
        return solveSVD(this, b, out);
    }

    /**
     * Computes U * S * V^T and produces the matrix this decomposition represents
     * 
     * @param {AbstractMat} [out] - The output matrix. If not specified, a new matrix will be created
     * @returns {AbstractMat} The result U * S * V^T
     */
    toMat(out) {
        var U = this.U;
        var V = this.V;
        var S = this.S;

        var M = U.rows();
        var N = V.cols();
        out = out !== undefined ? out : similar(U, M, N);

        if (out.rows() !== M || out.cols() !== N) {
            throw "Output matrix has wrong dimensions";
        }

        mult(U, mult(Diagonal.new(toArray(S), M, N), transpose(V)), out);

        return out;

    }

    /**
     * Creates the pseudo inverse of this SVD
     * If the matrix has full rank, this will equal the inverse
     * 
     * @returns {SVD} The SVD of the inverse
     */
    inv() {
        var U = copy(this.U);
        var V = copy(this.V);
        var S = copy(this.S);

        // replace every non zero singular value with reciprocal
        for (let i = 0; i < S.rows(); i++) {
            let si = S.at(i);
            // TODO better epsilon
            if (si < 1E-7) {
                // arrived at first zero singular value
                // since they are ordered we are finished
                break;
            }

            S.set(1.0 / si, i);

        }

        return SVD.new(V, S, U);
    }
}
//*****************************************************
/**
 * Computes the singular value decomposition of a general matrix.
 * 
 * A MxN matrix A can be decomposed into U * S * V^T
 * 
 * U is a unitary MxM matrix
 * V is a unitary NxN matrix
 * S is a MxN diagonal matrix
 * 
 * The diagonal entries of S are the singular values of A and the columns of U and V are the corresponding left and right singular vectors.
 * The singular values and their singular vectors are sorted in ascending order. S is represented as a vector containing the singular values.
 * 
 * Partially adapted from GSL library and based on the SVD algorithm from "Matrix Computations" in section 8.6
 * @param {AbstractMat} m - The input matrix
 * @returns {SVD} The SVD
 */
function computeSVD(m) {

    // in case of m < n -> transpose temporarily
    var isTransposed = false;

    if (m.rows() < m.cols()) {
        m = transpose(m);
        isTransposed = true;
    }

    let M = m.rows();
    let N = m.cols();

    let K = Math.min(M, N);

    var ubv = unpackUBV(computeUBVD(m));



    // diagonal
    var d = copy(diag(ubv.B));
    var fac = TypedMatFactory.newFromMat(ubv.B);
    // superdiagonal
    var work = fac.uninitialized(N, 1);

    var sd = subvec(work, 1);
    for (let i = 0; i < N - 1; i++) {
        sd.set(ubv.B.at(i, i + 1), i);
    }

    // remove very small elements
    chopSmallElements(d, sd);

    var U = copy(ubv.U);
    var V = ubv.Vt;

    // Adapted from GSL 
    var b = N - 1;
    var iter = 0;

    var a;

    while (b > 0) {
        var sdbm1 = sd.at(b - 1);

        if (sdbm1 === 0.0 || isNaN(sdbm1)) {
            b--;
            continue;
        }
        /* Find the largest unreduced block (a,b) starting from b
         and working backwards */

        a = b - 1;

        while (a > 0) {
            var sdam1 = sd.at(a - 1);

            if (sdam1 === 0.0 || isNaN(sdam1)) {
                break;
            }

            a--;
        }

        iter++;

        if (iter > 100 * N) {
            throw "SVD did not converge";
        }

        {
            var nBlock = b - a + 1;


            var dBlock = subvec(d, a, nBlock);
            var sdBlock = subvec(sd, a, nBlock - 1);

            var uBlock = block(U, 0, a, U.rows(), nBlock);

            var vBlock = block(V, 0, a, V.rows(), nBlock);

            qrstep(dBlock, sdBlock, uBlock, vBlock);
        }
    }


    for (let j = 0; j < K; j++) {
        let Sj = d.at(j);

        if (Sj < 0.0) {
            for (let i = 0; i < N; i++) {
                let Vij = V.at(i, j);
                V.set(-Vij, i, j);
            }

            d.set(-Sj, j);
        }
    }

    // sort singular values
    for (let i = 0; i < K; i++) {
        let maxSigma = argmax(subvec(d, i));

        // argmax starts from offset
        var rowMax = maxSigma.row + i;
        if (rowMax !== i) {
            // swap singular values stored in col vec
            swapRow(d, i, rowMax);
            // swap cols
            swapCol(V, i, rowMax);
            swapCol(U, i, rowMax);
        }
    }


    // this method return U,V
    // A = USV^T -> A^T = VSU^T = B
    // B = U_bS_bV_b^T 
    // => U_b = V
    // => S_b = S
    // => V_b^T = U^T
    // Since we return V instead of V^T, we have V_b = U 
    // redo transpose
    if (isTransposed) {

        return SVD.new(V, d, U);

    }
    else {
        return SVD.new(U, d, V);
    }

}
//*****************************************************
// pad strings
function pad(s, n) {
    if (s.length >= n) {
        return s;
    }

    while (s.length < n) {
        s += " ";
    }
    return s;
}
//*****************************************************
/**
 * Creates a pretty printable version of a matrix
 * @param {AbstractMat} m - The input matrix
 * @returns {string} A printable representation of the input matrix
 */
function prettyprint(m) {
    var mstr = map(m, x => x.toString(), MatAny.uninitialized(m.rows(), m.cols()));

    // find maximum length per col
    var ls = colreduce(mstr, x => reduce(x, (acc, v) => Math.max(acc, v.length), 0),
        TypedMatFactory.new(Int16Array).uninitialized(1, m.cols()));

    map(mstr, (x, i, j) => pad(x, ls.at(0, j)), mstr);
    var rows = rowreduce(mstr, x => toArray(x).join(" "), MatAny.uninitialized(m.rows(), 1));

    return toArray(rows).join("\n");

}
//*****************************************************
/**
 * Creates a pretty printable version of a matrix
 * @see prettyprint
 * @param {AbstractMat} m - The input matrix
 * @returns {string} A printable representation of the input matrix
 */
function toString(m) {
    return prettyprint(m);
}
//*****************************************************
var currentDefaultType = f32;
//*****************************************************
/**
 * Get the default storage type used for generic matrices
 * @returns {object} The current default storage type
 */
function getDefaultType() {
    return currentDefaultType;
}
//*****************************************************
/**
 * Sets the default storage type.
 * 
 * At the beginning, this will be Float32Array
 * @param {object} type - The new storage type
 */
function setDefaultType(type) {
    currentDefaultType = type;
}
//*****************************************************

/*
 * Geometry
 */

/**
 * Computes the cross product between two 3D vectors a x b
 * 
 * @param {AbstractMat} a - The first vector
 * @param {AbstractMat} b - The second vector
 * @param {AbstractMat} [out] - The output matrix. If not specified, a new matrix will be created.
 * @returns {AbstractMat} a x b
 */
function cross(a, b, out) {
    if (!isVec(a) || !isVec(b)) {
        throw "Cross product only defined on vectors";
    }
    if (a.rows() !== 3 || b.rows() !== 3) {
        throw "Cross product only defined in 3D";
    }


    var a0 = a.at(0);
    var a1 = a.at(1);
    var a2 = a.at(2);

    var b0 = b.at(0);
    var b1 = b.at(1);
    var b2 = b.at(2);

    out = out !== undefined ? out : similar(a);

    out.set(a1 * b2 - a2 * b1, 0);
    out.set(a2 * b0 - a0 * b2, 1);
    out.set(a0 * b1 - a1 * b0, 2);

    return out;

}

/**
 * Computes a 4x4 view matrix for 3D space
 * 
 * @param {AbstractMat} eye - The camera center
 * @param {AbstractMat} center - The point to look at
 * @param {AbstractMat} up - The up vector
 * @returns {Mat} The view matrix
 */
function lookAt(eye, center, up) {
    var z = sub(eye, center);
    normalize(z, z);

    var x = cross(up, z);
    normalize(x, x);

    var y = cross(z, x);

    var V = setId(similar(eye, 4, 4));

    var R = block(V, 0, 0, 3, 3);

    insert(row(R, 0), transpose(x));
    insert(row(R, 1), transpose(y));
    insert(row(R, 2), transpose(z));

    var T = setId(similar(eye, 4, 4));
    insert(block(T, 0, 3, 3, 1), neg(eye));

    return mult(V, T);

}

/**
 * Computes a 4x4 orthographic projection for 3D space
 * 
 * Note: This does not include a flip in the z-direction
 * 
 * @param {number} left - Left plane
 * @param {number} right - Right plane
 * @param {number} bottom - Bottom plane
 * @param {number} top - Top plane
 * @param {number} near - Near plane
 * @param {number} far - Far plane
 * @returns {Mat} The orthographic projection
 */
function ortho(left, right, bottom, top, near, far) {
    var P = setId(mat(4, 4, f32));

    P.set(2.0 / (right - left), 0, 0);
    P.set(2.0 / (top - bottom), 1, 1);
    P.set(2.0 / (far - near), 2, 2);

    var c = col(P, 3);
    c.set(- (left + right) / (right - left), 0);
    c.set(-(bottom + top) / (top - bottom), 1);
    c.set(- (far + near) / (far - near), 2);

    return P;
}

/**
 * Computes a 4x4 central perspective matrix.
 * 
 * This implements the pinhole camera with non-linear z-scale
 * 
 * @param {number} near -  Near plane
 * @param {number} far - Far plane 
 * @returns {Mat} The central perspective matrix
 */
function centralPerspective(near, far) {
    var p = setId(mat(4, 4, f32));
    insert(diag(p), vecFrom([near, near, near + far, 0], f32));

    p.set(1.0, 3, 2);
    p.set(-near * far, 2, 3);

    return p;
}

/**
 * Computes a 4x4 frustum projection.
 * 
 * Note: This includes a z-coordinate flip
 * 
 * @param {number} left - Left plane
 * @param {number} right - Right plane
 * @param {number} bottom - Bottom plane
 * @param {number} top - Top plane
 * @param {number} near - Near plane
 * @param {number} far - Far plane
 * @returns {Mat} The frustum projection
 */
function frustum(left, right, bottom, top, near, far) {
    var P = setId(mat(4, 4, f32));

    P.set(2.0 * near / (right - left), 0, 0);
    P.set(2.0 * near / (top - bottom), 1, 1);
    P.set(-(far + near) / (far - near), 2, 2);
    P.set(0.0, 3, 3);


    P.set((left + right) / (right - left), 0, 2);
    P.set((bottom + top) / (top - bottom), 1, 2);

    P.set(-1.0, 3, 2);

    P.set(-2.0 * far * near / (far - near), 2, 3);


    return P;
}

/**
 * Computes a 4x4 perspective matrix given a field of view
 * 
 * Note: This includes a z-coordinate flip
 * 
 * @param {number} fov - The full field of view
 * @param {number} aspect - The aspect ratio width/height
 * @param {number} near - The near plane
 * @param {number} far - The far plane
 * @returns {Mat} The perspective matrix
 */
function perspective(fov, aspect, near, far) {

    var xsize = near * Math.tan(fov * 0.5);
    var ysize = xsize / aspect;

    var O = ortho(-xsize, xsize, -ysize, ysize, near, far);

    var rtol = setId(mat(4, 4, f32));
    diag(rtol).set(-1.0, 2);

    var persp = centralPerspective(near, far);

    return mult(O, mult(persp, rtol));
}


/**
 * Computes a 4x4 viewport transform matrix
 * 
 * @param {number} x0 - The origin x 
 * @param {number} y0 - The origin y
 * @param {number} w - The width of the viewport
 * @param {number} h - The height of the viewport
 * @param {boolean} [flipy=false] - If true, the viewport will be flipped along its y-axis. This is needed if the drawing surface's coordinate system starts on the upper left
 * @returns {Mat} The viewport matrix  
 */
function viewport(x0, y0, w, h, flipy) {
    var V = setId(mat(4, 4, f32));

    var d = diag(V);
    d.set(w / 2.0, 0);

    if (flipy) {
        d.set(-h / 2.0, 1);

    }
    else {
        d.set(h / 2.0, 1);

    }
    d.set(0.5, 2);

    var c = col(V, 3);

    c.set(w / 2 + x0, 0);
    if (flipy) {
        c.set(h / 2 - y0, 1);

    }
    else {

        c.set(h / 2 + y0, 1);
    }
    c.set(0.5, 2);

    return V;
}

/**
 * Homogenizes a matrix storing points in its columns
 * 
 * For each column, this function will divide the column by its last entry
 * 
 * @param {AbstractMat} points - The input points
 * @returns {AbstractMat} points
 */
function homogenize(points) {
    colwise(points, col => {
        scale(col, 1.0 / col.at(col.rows() - 1), col);
    });
    return points;
}

/**
 * Computes a 3x3 cross product matrix K from a vector k such that for any vector v: Kv = cross(k,v)
 * 
 * @param {AbstractMat} k - A 3D vector
 * @returns {Mat} The cross prodcut matrix
 */
function crossMatrix(k) {
    if (!isVec(k) || k.rows() !== 3) {
        throw "k has to be a 3D vector to compute the cross product matrix";
    }
    var K = setZero(similar(k, 3, 3));

    K.set(-k.at(2), 0, 1);
    K.set(k.at(2), 1, 0);

    K.set(k.at(1), 0, 2);
    K.set(-k.at(1), 2, 0);

    K.set(-k.at(0), 1, 2);
    K.set(k.at(0), 2, 1);

    return K;
}

/**
 * Computes a 3x3 rotation matrix, which represents a rotaion around an axis
 * 
 * @param {AbstractMat} axis - The axis to rotate around
 * @param {number} angle - The angle to rotate in rad
 * @returns {Mat} The rotation matrix
 */
function axisAngle(axis, angle) {
    axis = copy(axis);
    normalize(axis, axis);
    // Rodrigues' formula
    var K = crossMatrix(axis);
    var K2 = mult(K, K);

    scale(K2, 1.0 - Math.cos(angle), K2);
    scale(K, Math.sin(angle), K);

    add(K, K2, K);
    var R = add(K, id(3, 3), K);

    return R;
}

/**
 * Computes a 4x4 3D rotation matrix, which represents a rotaion around an axis
 * 
 * @param {AbstractMat} axis - The axis to rotate around
 * @param {number} angle - The angle to rotate in rad
 * @returns {Mat} The rotation matrix
 */
function axisAngle4(axis, angle) {
    var R = axisAngle(axis, angle);

    var R4 = setId(similar(R, 4, 4));
    insert(block(R4, 0, 0, 3, 3), R);

    return R4;
}



/**
 * Converts degrees to radians with degrees*PI/180
 * 
 * @param {number} deg - The angle in degrees
 * @returns {number} The angle in radians
 */
function deg2rad(deg) {
    return deg * Math.PI / 180.0;
}
/**
 * Converts radians to degrees with radians*180/PI
 * 
 * @param {number} rad - The angle in radians
 * @returns {number} The angle in degrees
 */
function rad2deg(rad) {
    return rad * 180.0 / Math.PI;
}

/**
 * Creates a 4x4 translation matrix for a given translation vector
 * 
 * @param {AbstractMat} t - 3D translation vector
 * @returns {Mat} A translation matrix
 */
function translation(t) {
    if (!isVec(t) || t.rows() !== 3) {
        throw "Translation vector needs to be 3-dimensional";
    }
    var T = setId(similar(t, 4, 4));

    var subcol = subvec(col(T, 3), 0, 3);
    insert(subcol, t);

    return T;
}

/**
 * Creates a 4x4 scaling matrix for a given scaling vector. 
 * This vector contains the scaling factors for each dimension
 * 
 * @param {AbstractMat} s - 3D scaling vector
 * @returns {Mat} The scaling matrix
 */
function scaling(s) {
    if (!isVec(s) || s.rows() !== 3) {
        throw "Scaling vector needs to be 3-dimensional";
    }
    var S = similar(s, 4, 4);

    var subdiag = subvec(diag(S), 0, 3);
    insert(subdiag, s);

    return S;
}

/**
 * 3D x vector [1,0,0]
 */
var X = VecF32.from([1, 0, 0]);
/**
 * 3D y vector [0,1,0]
 */
var Y = VecF32.from([0, 1, 0]);
/**
 * 3D z vector [0,0,1]
 */
var Z = VecF32.from([0, 0, 1]);

/**
 * Homogenous 3D x vector [1,0,0,0]
 */
var Xh = VecF32.from([1, 0, 0, 0]);
/**
 * Homogenous 3D y vector [0,1,0,0]
 */
var Yh = VecF32.from([0, 1, 0, 0]);
/**
 * Homogenous 3D z vector [0,0,1,0]
 */
var Zh = VecF32.from([0, 0, 1, 0]);


export {
    f32,
    f64,
    i8,
    ui8,
    ui8c,
    i16,
    ui16,
    i32,
    ui32,
    i64,
    ui64,
    any,
    Mat,
    TypedMatFactory,
    MatF32,
    MatF64,
    MatAny,
    MatI8,
    MatUI8,
    MatUI8Clamped,
    MatI16,
    MatUI16,
    MatI32,
    MatUI32,
    MatI64,
    MatUI64,
    mat,
    VecF32,
    VecF64,
    VecAny,
    VecI8,
    VecUI8,
    VecUI8Clamped,
    VecI16,
    VecUI16,
    VecI32,
    VecUI32,
    VecI64,
    VecUI64,
    vec,
    Diagonal,
    TypedDiagonalFactory,
    DiagonalF32,
    DiagonalF64,
    DiagonalAny,
    DiagonalI8,
    DiagonalUI8,
    DiagonalUI8Clamped,
    DiagonalI16,
    DiagonalUI16,
    DiagonalI32,
    DiagonalUI32,
    DiagonalI64,
    DiagonalUI64,
    DiagonalView,
    BlockView,
    ColumnView,
    RowView,
    TransposeView,
    PaddedView,
    MinorView,
    TriangularMode,
    TriangularView,
    RowPermutation,
    isVec,
    add,
    addScalar,
    sub,
    mult,
    fill,
    insert,
    setId,
    setOne,
    setZero,
    ones,
    zeros,
    id,
    rand,
    copy,
    similar,
    from,
    vecFrom,
    diag,
    transpose,
    block,
    row,
    col,
    subvec,
    map,
    convert,
    reduce,
    sum,
    neg,
    abs,
    cwiseMult, cwiseDiv,
    dot,
    scale,
    absSum,
    sqrSum,
    trace,
    max,
    min,
    argmax,
    argmin,
    norm,
    norm2,
    norm2Squared,
    norm1,
    normInf,
    normFrobenius,
    normalize,
    toArray,
    colreduce,
    rowreduce,
    rowwise,
    colwise,
    det,
    inv,
    rank,
    cond,
    computePLUD,
    computeSVD,
    computeUBVD,
    unpackUBV,
    solve,
    solvePLU,
    solveSVD,
    solveLowerTriagonal,
    solveUpperTriagonal,
    prettyprint,
    toString,
    getDefaultType,
    setDefaultType,
    rad2deg,
    deg2rad,
    cross,
    lookAt,
    ortho,
    frustum,
    perspective,
    centralPerspective,
    viewport,
    homogenize,
    crossMatrix,
    axisAngle,
    axisAngle4,
    translation,
    scaling,
    X, Y, Z,
    Xh, Yh, Zh
};