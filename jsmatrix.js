/**
 * Storage types used for matrices
 */
const f32 = Float32Array;
const f64 = Float64Array;
const i8 = Int8Array;
const ui8 = Uint8Array;
const ui8c = Uint8ClampedArray;
const i16 = Int16Array;
const ui16 = Uint16Array;
const i32 = Int32Array;
const ui32 = Uint32Array;
const i64 = BigInt64Array;
const ui64 = BigUint64Array;
const any = Array;

//*****************************************************
let currentDefaultType = f32;

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
    constructor(data, rows, cols, innerStride = 1, outerStride = rows) {
        /** @private */
        this._data = data;
        /** @private */
        this._rows = rows;
        /** @private */
        this._cols = cols;
        /** @private */
        this._innerStride = innerStride;
        /** @private */
        this._outerStride = outerStride;
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
            throw new Error("Trying to access data outside of range");
        }

        return this._data[idx];
    }

    set(v, i, j) {
        const idx = this.index(i, j);
        if (idx > this._data.length) {
            throw new Error("Trying to access data outside of range");
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
        const n = rows * cols;
        const data = new this._type(n);
        return Mat.new(data, rows, cols);
    }

    from(data, rows, cols) {
        const n = rows * cols;
        const datat = new this._type(n);

        for (let i = 0; i < n; i++) {
            datat[i] = data[i];
        }
        return Mat.new(datat, rows, cols);
    }

    /**
     * Creates a new matrix from an array of rows
     * @param {Array<AbstractMat>} rows - An array of row vectors. All vectors need to have the same size
     * @returns {AbstractMat} A newly created matrix
     */
    fromRows(rows) {
        const r = rows.length;
        if (r < 1) {
            throw new Error("At least one row needed");
        }

        const c = rows[0].cols();

        const result = this.uninitialized(r, c);

        for (let i = 0; i < r; i++) {
            insert(row(result, i), rows[i]);
        }

        return result;

    }

    /**
     * Creates a new matrix from an array of columns
     * @param {Array<AbstractMat>} cols - An array of column vectors. All vectors need to have the same size
     * @returns {AbstractMat} A newly created matrix
     */
    fromCols(cols) {
        const c = cols.length;
        if (c < 1) {
            throw new Error("At least one column needed");
        }

        const r = cols[0].rows();

        const result = this.uninitialized(r, c);

        for (let i = 0; i < c; i++) {
            insert(col(result, i), cols[i]);
        }

        return result;

    }

    all(v, rows, cols) {
        const n = rows * cols;
        const data = new this._type(n);
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

    id(rows, cols = rows) {
        const m = this.zeros(rows, cols);

        const dm = diag(m);
        fill(dm, 1);

        return m;
    }

    copy(m) {
        const c = m.cols();
        const r = m.rows();
        const n = c * r;

        const data = new this._type(n);

        let idx = 0;
        for (let j = 0; j < c; j++) {
            for (let i = 0; i < r; i++) {
                data[idx] = m.at(i, j);
                idx++;
            }
        }
        return Mat.new(data, r, c);
    }

    rand(rows, cols) {
        const r = rows;
        const c = cols;
        const n = c * r;

        const data = new this._type(n);

        let idx = 0;
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
const MatF32 = new TypedMatFactory(Float32Array);
/**
 * A TypedMatFactory to create Float64Array based matrices
 */
const MatF64 = new TypedMatFactory(Float64Array);
/**
 * A TypedMatFactory to create generic Array based matrices
 */
const MatAny = new TypedMatFactory(Array);
/**
 * A TypedMatFactory to create Int8Array based matrices
 */
const MatI8 = new TypedMatFactory(Int8Array);
/**
 * A TypedMatFactory to create Uint8Array based matrices
 */
const MatUI8 = new TypedMatFactory(Uint8Array);
/**
 * A TypedMatFactory to create Uint8ClampedArray based matrices
 */
const MatUI8Clamped = new TypedMatFactory(Uint8ClampedArray);

/**
 * A TypedMatFactory to create Int16Array based matrices
 */
const MatI16 = new TypedMatFactory(Int16Array);
/**
 * A TypedMatFactory to create Uint16Array based matrices
 */
const MatUI16 = new TypedMatFactory(Uint16Array);

/**
 * A TypedMatFactory to create Int32Array based matrices
 */
const MatI32 = new TypedMatFactory(Int32Array);
/**
 * A TypedMatFactory to create Uint32Array based matrices
 */
const MatUI32 = new TypedMatFactory(Uint32Array);

/**
 * A TypedMatFactory to create BigInt64Array based matrices
 */
const MatI64 = new TypedMatFactory(BigInt64Array);
/**
 * A TypedMatFactory to create BigUint64Array based matrices
 */
const MatUI64 = new TypedMatFactory(BigUint64Array);

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

    from(data, n = data.length) {
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
            throw new Error("Vector factory can't copy from non-vector");
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
const VecF32 = new TypedVecFactory(Float32Array);
/**
 * A TypedVecFactory to create Float64Array based vectors
 */
const VecF64 = new TypedVecFactory(Float64Array);
/**
 * A TypedVecFactory to create generic Array based vectors
 */
const VecAny = new TypedVecFactory(Array);
/**
 * A TypedVecFactory to create Int8Array based vectors
 */
const VecI8 = new TypedVecFactory(Int8Array);
/**
 * A TypedVecFactory to create Uint8Array based vectors
 */
const VecUI8 = new TypedVecFactory(Uint8Array);
/**
 * A TypedVecFactory to create Uint8ClampedArray based vectors
 */
const VecUI8Clamped = new TypedVecFactory(Uint8ClampedArray);

/**
 * A TypedVecFactory to create Int16Array based vectors
 */
const VecI16 = new TypedVecFactory(Int16Array);
/**
 * A TypedVecFactory to create Uint16Array based vectors
 */
const VecUI16 = new TypedVecFactory(Uint16Array);

/**
 * A TypedVecFactory to create Int32Array based vectors
 */
const VecI32 = new TypedVecFactory(Int32Array);
/**
 * A TypedVecFactory to create Uint32Array based vectors
 */
const VecUI32 = new TypedVecFactory(Uint32Array);

/**
 * A TypedVecFactory to create BigInt64Array based vectors
 */
const VecI64 = new TypedVecFactory(BigInt64Array);
/**
 * A TypedVecFactory to create BigUint64Array based vectors
 */
const VecUI64 = new TypedVecFactory(BigUint64Array);

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
            throw new Error("Not enough data for diagonal construction");
        }
    }

    static new(data, rows, cols, stride) {
        return new Diagonal(data, rows, cols, stride);
    }

    static fromMat(m, rows, cols, stride) {
        if (!isVec(m)) {
            throw new Error("Diagonal only takes vectors as input");
        }

        const data = toArray(m);
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
            throw new Error("Trying to set offdiagonal element of diagonal matrix");
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


    from(data, rows, cols = rows) {
        const n = Math.min(rows, cols);
        const datat = new this._type(n);
        for (let i = 0; i < n; i++) {
            datat[i] = data[i];
        }

        return Diagonal.new(datat, rows, cols);
    }

    uninitialized(rows, cols) {
        cols = cols !== undefined ? cols : rows;
        const n = Math.min(rows, cols);
        const data = new this._type(n);
        return Diagonal.new(data, rows, cols);
    }

    all(v, rows, cols = rows) {
        const n = Math.min(rows, cols);
        const data = new this._type(n);
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
        const n = Math.min(rows, cols);
        const data = new this._type(n);
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
const DiagonalF32 = new TypedDiagonalFactory(Float32Array);

/**
 * A TypedDiagonalFactory to create Float64Array based diagonal matrices
 */
const DiagonalF64 = new TypedDiagonalFactory(Float64Array);

/**
 * A TypedDiagonalFactory to create generic Array based diagonal matrices
 */
const DiagonalAny = new TypedDiagonalFactory(Array);

/**
 * A TypedDiagonalFactory to create Int8Array based matrices
 */
const DiagonalI8 = new TypedDiagonalFactory(Int8Array);
/**
 * A TypedDiagonalFactory to create Uint8Array based matrices
 */
const DiagonalUI8 = new TypedDiagonalFactory(Uint8Array);
/**
 * A TypedDiagonalFactory to create Uint8ClampedArray based matrices
 */
const DiagonalUI8Clamped = new TypedDiagonalFactory(Uint8ClampedArray);

/**
 * A TypedDiagonalFactory to create Int16Array based matrices
 */
const DiagonalI16 = new TypedDiagonalFactory(Int16Array);
/**
 * A TypedDiagonalFactory to create Uint16Array based matrices
 */
const DiagonalUI16 = new TypedDiagonalFactory(Uint16Array);

/**
 * A TypedDiagonalFactory to create Int32Array based matrices
 */
const DiagonalI32 = new TypedDiagonalFactory(Int32Array);
/**
 * A TypedDiagonalFactory to create Uint32Array based matrices
 */
const DiagonalUI32 = new TypedDiagonalFactory(Uint32Array);

/**
 * A TypedDiagonalFactory to create BigInt64Array based matrices
 */
const DiagonalI64 = new TypedDiagonalFactory(BigInt64Array);
/**
 * A TypedDiagonalFactory to create BigUint64Array based matrices
 */
const DiagonalUI64 = new TypedDiagonalFactory(BigUint64Array);

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
        const r = m.rows();
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
        const c = m.cols();
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
const TriangularMode = {
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
                } else {
                    return this._m.at(i, j);
                }
            case TriangularMode.LOWER:
                if (i < j) {
                    return 0;
                } else {
                    return this._m.at(i, j);
                }
            case TriangularMode.STRICTLY_UPPER:
                if (i >= j) {
                    return 0;
                } else {
                    return this._m.at(i, j);
                }
            case TriangularMode.STRICTLY_LOWER:
                if (i <= j) {
                    return 0;
                } else {
                    return this._m.at(i, j);
                }

            case TriangularMode.UNIT_LOWER:
                if (i < j) {
                    return 0;
                } else if (i === j) {
                    return 1;
                } else {
                    return this._m.at(i, j);
                }
            case TriangularMode.UNIT_UPPER:
                if (i >= j) {
                    return 0;
                } else if (i === j) {
                    return 1;
                } else {
                    return this._m.at(i, j);
                }
            default:
                break;
        }

    }

    set(v, i, j) {
        if (!this.checkIndex(i, j)) {
            throw new Error("Trying to set value outside triangular part");
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
        } else {
            return solveLowerTriagonal(this, b, out);
        }
    }
}

//*****************************************************

/**
 * Represents a row permutation matrix without explicit storage of the matrix.
 *
 * Permutations are represented by a permutation table given as an array.
 * The underlying type can be set explicitly.
 *
 * @extends AbstractMat
 */
class RowPermutation {
    constructor(data, rows = data.length, cols = rows, type = currentDefaultType) {
        /** @private */
        this._data = data;
        /** @private */
        this._rows = rows;
        /** @private */
        this._cols = cols;
        this._type = type;
    }

    static new(data, rows, cols, type) {
        return new RowPermutation(data, rows, cols, type);
    }

    at(i, j = 0) {
        // permute rows then map to identity
        i = this._data[i];

        return i === j ? 1 : 0;
    }

    set(/**v,i,j */) {
        throw new Error("Permuation matrix is immutable");
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
            throw new Error("Trying to remove row/column outside of source");
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
    constructor(m, rows, cols, offDiagonalValue = 0, diagonalValue = 1) {
        /** @private */
        this._m = m;
        /** @private */
        this._rows = rows;
        /** @private */
        this._cols = cols;

        if (this._rows < m.rows() || this._cols < m.cols()) {
            throw new Error("Padded View needs to be at least as big as source");
        }
        /** @private */
        this._offDiagonalValue = offDiagonalValue;
        /** @private */
        this._diagonalValue = diagonalValue;
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
        if (i < this._m.rows() && j < this._m.cols()) {
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

        throw new Error("Trying to set values in padding");
    }

    type() {
        return this._m.type();
    }
}

//*****************************************************
/**
 * Creates an uninitialized matrix similar to the given one
 *
 * The copy will have the same size and underlying type as m
 * @param {AbstractMat} m - The input matrix
 * @param [rows] - Number of rows of the resulting matrix
 * @param [cols]- Number of cols of the resulting matrix
 * @returns {AbstractMat} An uninitialized matrix of the same size and type as m
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
}

/**
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

    const r = a.rows();
    const c = a.cols();

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


    if (a.cols() !== b.cols() || a.rows() !== b.rows()) {
        throw new Error("Trying to add input of different sizes");
    }
    const r = a.rows();
    const c = a.cols();

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

    if (a.cols() !== b.cols() || a.rows() !== b.rows()) {
        throw new Error("Trying to add input of different sizes");
    }

    const r = a.rows();
    const c = a.cols();

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
    const r = a.rows();
    const c = a.cols();
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
 * Computes a componenwise minimum of two matrices of the same size
 *
 * Computes (i,j) => Math.min(a.at(i,j),b.at(i,j))
 *
 * @param {AbstractMat} a - The first input matrix
 * @param {AbstractMat} b - The second input matrix
 * @param {AbstractMat} out - The output matrix. If not specified, a new matrix will be created
 */
function cwiseMin(a, b, out) {

    return map(a, (v, row, col) => Math.min(v, b.at(row, col)), out);
}

//*****************************************************
/**
 * Computes a componenwise maximum of two matrices of the same size
 *
 * Computes (i,j) => Math.max(a.at(i,j),b.at(i,j))
 *
 * @param {AbstractMat} a - The first input matrix
 * @param {AbstractMat} b - The second input matrix
 * @param {AbstractMat} out - The output matrix. If not specified, a new matrix will be created
 */
function cwiseMax(a, b, out) {

    return map(a, (v, row, col) => Math.max(v, b.at(row, col)), out);
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
function block(m, i, j, rows = m.rows() - i, cols = m.cols() - j) {
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
    const n = a.rows();
    const m = a.cols();

    if (m !== b.rows()) {
        throw new Error("Incompatible matrix dimensions");
    }
    const p = b.cols();

    out = out !== undefined ? out : similar(a, n, p);

    if (out.rows() !== n || out.cols() !== p) {
        throw new Error("Output dimension does not match input");
    }

    for (let j = 0; j < p; j++) {
        for (let i = 0; i < n; i++) {
            let s = 0;
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
    const r = m.rows();
    const c = m.cols();

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
    const r = m.rows();
    const c = m.cols();

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
    const r = m.rows();
    const c = m.cols();

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
    const out = TypedMatFactory.new(type).uninitialized(m.rows(), m.cols());

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
    let no_init = initValue === undefined;

    const r = m.rows();
    const c = m.cols();
    if (no_init && r * c === 0) {
        throw new Error("Calling reduce on empty input without initial value");
    }

    let accum = no_init ? m.at(0, 0) : initValue;


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
    const svd = computeSVD(m);
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
    const s = svd.S;

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
    const svd = computeSVD(m);
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

    const s = svd.S;
    const smax = s.at(0);
    const smin = s.at(s.rows() - 1);

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
    const c = m.cols();
    const r = m.rows();
    const result = [];
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
 * @param type - The type of the matrix
 * @returns {Diagonal} The identity matrix
 */
function id(rows, cols, type = currentDefaultType) {
    return TypedDiagonalFactory.new(type).id(rows, cols);
}

//*****************************************************
/**
 * Creates a zero matrix with the currently set default type
 *
 * @param {number} rows - The number of rows
 * @param {number} cols - The number of cols
 * @param type - The type of the matrix
 * @returns {Mat} The zero matrix
 */
function zeros(rows, cols, type = currentDefaultType) {
    return TypedMatFactory.new(type).zeros(rows, cols);
}

//*****************************************************
/**
 * Creates a random matrix
 *
 * @param {number} rows - The number of rows
 * @param {number} cols - The number of cols
 * @param type - The type of the matrix
 * @returns {Mat} The zero matrix
 */
function rand(rows, cols, type = currentDefaultType) {
    return TypedMatFactory.new(type).rand(rows, cols);
}

//*****************************************************
/**
 * Creates a ones matrix with the currently set default type
 *
 * @param {number} rows - The number of rows
 * @param {number} cols - The number of cols
 * @param type - The type of the matrix
 * @returns {Mat} The ones matrix
 */
function ones(rows, cols, type = currentDefaultType) {
    return TypedMatFactory.new(type).ones(rows, cols);
}

//*****************************************************
/**
 * Creates a new uninitialized matrix with the currently set default type
 *
 * @param {number} rows - The number of rows
 * @param {number} cols - The number of cols
 * @param type - The type of the matrix
 * @returns {Mat} An uninitialized matrix
 */
function mat(rows, cols, type = currentDefaultType) {
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
function from(data, rows, cols, type = currentDefaultType) {
    return TypedMatFactory.new(type).copy(Mat.new(data, rows, cols));
}

//*****************************************************
/**
 * Creates a vector from data with a given underlying type
 * @param {Array | TypedArray} data - The input data
 * @param {Object} [type] - The underlying type. Defaults to the current default type
 * @returns {Mat} A vector with a typed copy of the given data
 */
function vecFrom(data, type = currentDefaultType) {
    const n = data.length;
    return TypedMatFactory.new(type).copy(Mat.new(data, n, 1));
}

//*****************************************************
/**
 * Creates a new uninitialized vector with the currently set default type
 *
 * @param {number} n - The number of elements
 * @param type - The type of the matrix
 * @returns {Mat} An uninitialized vector
 */
function vec(n, type = currentDefaultType) {
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
const normFrobenius = norm2;

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
    const n = norm(m);

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
 * @param {Mat} b - The matrix to be read from
 * @returns {AbstractMat} a
 */
function insert(a, b) {
    if (a.rows() !== b.rows() || a.cols() !== b.cols()) {
        throw new Error("Insertion failed: Source and target have different sizes");
    }

    const c = a.cols();
    const r = b.rows();

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
    const r0 = copy(row(a, row0));
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
    const c0 = copy(col(a, col0));
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
        throw new Error("Determinant only defined for square sources");
    }

    // check if m has a dedicated det function
    if (m.det && m.det instanceof Function) {
        return m.det();
    }

    const n = m.rows();

    if (n === 1) {
        return m.at(0, 0);
    } else if (n === 2) {
        return m.at(0, 0) * m.at(1, 1) - m.at(0, 1) * m.at(1, 0);
    } else if (n === 3) {
        const a = m.at(0, 0);
        const b = m.at(0, 1);
        const c = m.at(0, 2);

        const d = m.at(1, 0);
        const e = m.at(1, 1);
        const f = m.at(1, 2);

        const g = m.at(2, 0);
        const h = m.at(2, 1);
        const i = m.at(2, 2);


        return a * e * i + b * f * g + c * d * h - c * e * g - b * d * i - a * f * h;

    } else {

        // TODO use more efficient formulation

        const lu = computePLUD(m);

        // singular matrix
        if (!lu) {
            return 0.0;
        }
        // lower diagonal is 1s -> determinant is 1
        // total determinant is therefore just product of u diagonal
        let s = (lu.numSwaps % 2 === 0 ? 1 : -1) * reduce(diag(lu.U), (acc, v) => acc * v);
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
    const rows = m.rows();
    const cols = m.cols();

    const brows = b.rows();
    const bcols = b.cols();


    if (rows !== cols) {
        throw new Error("Lower triagonal matrix not square");
    }

    if (rows !== brows) {
        throw new Error("b does not have correct size");
    }

    out = out !== undefined ? out : similar(b);

    for (let bc = 0; bc < bcols; bc++) {
        for (let i = 0; i < rows; i++) {
            let sumOffDiag = 0;
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
    const rows = m.rows();
    const cols = m.cols();

    const brows = b.rows();
    const bcols = b.cols();


    if (rows !== cols) {
        throw new Error("Lower triagonal matrix not square");
    }

    if (rows !== brows) {
        throw new Error("b does not have correct size");
    }

    out = out !== undefined ? out : similar(b);

    for (let bc = 0; bc < bcols; bc++) {
        for (let i = rows - 1; i >= 0; i--) {
            let sumOffDiag = 0;
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
        const plu = computePLUD(a);

        // singular matrix solve with svd instead
        if (!plu) {
            return computeSVD(a).solve(b, out);
        }
        return plu.solve(b, out);
    }

    // for non square matrices -> use svd

    const svd = computeSVD(a);
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
    const pb = mult(plu.P, b);

    const y = solveLowerTriagonal(plu.L, pb);
    const x = solveUpperTriagonal(plu.U, y, out);

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
    const U = svd.U;
    const V = svd.V;
    const S = svd.S;

    const N = V.rows();

    const P = b.cols();

    let r = 0;

    out = out !== undefined ? out : similar(b, N, P);
    setZero(out);

    for (let j = 0; j < P; j++) {
        const bj = col(b, j);
        const zj = col(out, j);

        for (r = 0; r < S.rows(); r++) {
            const sigma = S.at(r);
            // TODO relative epsilon?
            if (sigma <= 1E-7) {
                break;
            }

            const zr = dot(col(U, r), bj) / sigma;
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
        throw new Error("Inverse only defined for square sources");
    }

    // check if m has a dedicated inv function
    if (m.inv && m.inv instanceof Function) {
        return m.inv();
    }

    const n = m.rows();

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
    } else if (n === 3) {
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
    } else {

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
    static compute(a, out) {
        return computePLUD(a, out);
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
            throw new Error("Output has wrong size");
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
function computePLUD(a, out = similar(a)) {

    const n = a.rows();
    const r = a.cols();
    const rnmin = Math.min(a.rows(), a.cols());

    const permutes = TypedMatFactory.new(Int32Array).uninitialized(n, 1);
    map(permutes, (v, i) => i, permutes);

    insert(out, a);
    let numSwaps = 0;

    for (let k = 0; k < rnmin; k++) {

        // find largest value in colum for partial pivot

        let maxIndex = k;
        let maxEl = Math.abs(out.at(maxIndex, k));
        for (let l = k + 1; l < n; l++) {
            if (Math.abs(out.at(l, k)) > maxEl) {
                maxIndex = l;
                maxEl = Math.abs(out.at(l, k));
            }
        }

        // swap row k with maxIndex
        if (maxIndex !== k) {
            numSwaps++;
            swapRow(permutes, maxIndex, k);
            swapRow(out, maxIndex, k);
        }

        // Algorithm from "Matrix computations"

        const outkk = out.at(k, k);

        // singularity detected
        if (Math.abs(outkk) < 1E-7) {
            return null;
        }
        const subcol = subvec(col(out, k), k + 1);
        scale(subcol, 1.0 / outkk, subcol);

        // update lower block
        // case n > r
        if (k < r) {
            const rowRho = subrowvec(row(out, k), k + 1);

            for (let rho = k + 1; rho < n; rho++) {
                const subrow = subrowvec(row(out, rho), k + 1);
                const arhok = out.at(rho, k);
                sub(subrow, scale(rowRho, arhok), subrow);

            }
        }
    }

    const blockL = block(out, 0, 0, a.rows(), rnmin);
    const blockU = block(out, 0, 0, rnmin, a.cols());
    const L = TriangularView.new(blockL, TriangularMode.UNIT_LOWER);
    const U = TriangularView.new(blockU, TriangularMode.UPPER);
    const P = RowPermutation.new(toArray(permutes), permutes.rows(), permutes.rows(), a.type());
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
    const absa = Math.abs(a);
    const absb = Math.abs(b);

    const sqr = x => x * x;

    if (absa > absb) {

        return absa * Math.sqrt(1.0 + sqr(absb / absa));
    }

    return absb === 0.0 ? 0.0 : absb * Math.sqrt(1.0 + sqr(absa / absb));
}

//*****************************************************
function householderVector(x, out) {
    if (!isVec(x)) {
        throw new Error("Householder transform needs to operate on vector");
    }

    let v = out !== undefined ? out : copy(x);
    const m = x.rows();

    // compute squared length of subvector starting at i = 1

    let sigma = 0.0;
    for (let i = 1; i < m; i++) {
        let vi = x.at(i);
        sigma += vi * vi;
    }

    const x0 = x.at(0);
    v.set(1.0, 0);


    if (sigma === 0.0) {
        return { beta: 0.0, v: v };
    }

    const my = Math.sqrt(x0 * x0 + sigma);

    if (x0 <= 0.0) {
        v.set(x0 - my, 0);
    } else {
        v.set(-sigma / (x0 + my), 0);
    }


    const v0 = v.at(0);
    const v02 = v0 * v0;
    const beta = 2.0 * v02 / (sigma + v02);

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
 * @returns {Number }The sum of all multiplied corresponding elements in a and b
 */
function dot(a, b) {
    const r = a.rows();
    const c = a.cols();

    if (r !== b.rows() || c !== b.cols()) {
        throw new Error("Inputs must match in size for dot product");
    }

    let s = 0.0;

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
        throw new Error("Input for subvec needs to be a vector");
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
        throw new Error("Input for subvec needs to be a row vector");
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

    const r = a.rows();
    const c = a.cols();

    // apply per col
    for (let j = 0; j < c; j++) {
        // compute v^T * A[:,j]

        const aj = col(out, j);
        let wj = aj.at(0);

        wj += dot(subvec(v, 1), subvec(aj, 1));

        const bwj = beta * wj;

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

    const r = a.rows();
    const c = a.cols();

    // apply per row
    for (let i = 0; i < r; i++) {
        const ai = transpose(row(out, i));

        let wi = ai.at(0);

        wi += dot(subvec(v, 1), subvec(ai, 1));

        const bwi = beta * wi;

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

/* eslint-disable no-unused-vars */
// TODO use
function computeQR(a, out) {
    const M = a.rows();
    const N = a.cols();

    const P = Math.min(M, N);
    if (out !== undefined) {
        insert(out, a);
        a = out;
    } else {
        a = copy(a);
    }


    for (let j = 0; j < N; j++) {
        let cj = subvec(col(a, j), j);
        let house = householderVector(cj);

        const blockj = block(a, j, j, M - j, N - j);
        applyHouseholderLeft(house.beta, house.v, blockj, blockj);

        insert(subvec(cj, 1), subvec(house.v, 1, M - j - 1));

    }

}
/* eslint-enable no-unused-vars */

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
    const M = a.rows();
    const N = a.cols();

    if (M < N) {
        throw new Error("Biadiagonal only implemented for M >= N");
    }
    if (out !== undefined) {
        insert(out, a);
        a = out;
    } else {
        a = copy(a);
    }

    for (let j = 0; j < N; j++) {
        let cj = subvec(col(a, j), j);
        let house = householderVector(cj);

        const blockj = block(a, j, j, M - j, N - j);
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
    const fac = TypedMatFactory.newFromMat(ubv);
    const B = fac.zeros(ubv.rows(), ubv.cols());

    insert(diag(B), diag(ubv));
    // superdiag
    for (let i = 0; i < N - 1; i++) {
        B.set(ubv.at(i, i + 1), i, i + 1);
    }

    const U = fac.zeros(M, M);
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

    const Vt = fac.id(N);
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
function isSymmetric(a, eps = 0.0) {
    const n = a.rows();
    // only square matrices are symmetric
    if (n !== a.cols()) {
        return false;
    }
    for (let j = 0; j < n; j++) {
        for (let i = 0; i < n; i++) {
            if (Math.abs(a.at(i, j) - a.at(j, i)) > eps) {
                return false;
            }
        }
    }

    return true;
}
//*****************************************************
/**
 * Adapted from the public domain JAMA library (JAMA)
 * @param {AbstractMat} V
 * @param {AbstractMat} d
 * @param {AbstractMat} e
 */
function tred2(V, d, e) {
    const n = V.rows();
    //  This is derived from the Algol procedures tred2 by
    //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
    //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
    //  Fortran subroutine in EISPACK.

    for (let j = 0; j < n; j++) {
        d[j] = V.at(n - 1, j);
    }

    // Householder reduction to tridiagonal form.

    for (let i = n - 1; i > 0; i--) {

        // Scale to avoid under/overflow.

        let scale = 0.0;
        let h = 0.0;
        for (let k = 0; k < i; k++) {
            scale = scale + Math.abs(d[k]);
        }
        if (scale === 0.0) {
            e[i] = d[i - 1];
            for (let j = 0; j < i; j++) {
                d[j] = V.at(i - 1, j);
                V.set(0.0, i, j);
                V.set(0.0, j, i);
            }
        } else {

            // Generate Householder vector.

            for (let k = 0; k < i; k++) {
                d[k] /= scale;
                h += d[k] * d[k];
            }
            let f = d[i - 1];
            let g = Math.sqrt(h);
            if (f > 0) {
                g = -g;
            }
            e[i] = scale * g;
            h = h - f * g;
            d[i - 1] = f - g;
            for (let j = 0; j < i; j++) {
                e[j] = 0.0;
            }

            // Apply similarity transformation to remaining columns.

            for (let j = 0; j < i; j++) {
                f = d[j];
                V.set(f, j, i);
                g = e[j] + V.at(j, j) * f;
                for (let k = j + 1; k <= i - 1; k++) {
                    g += V.at(k, j) * d[k];
                    e[k] += V.at(k, j) * f;
                }
                e[j] = g;
            }
            f = 0.0;
            for (let j = 0; j < i; j++) {
                e[j] /= h;
                f += e[j] * d[j];
            }
            const hh = f / (h + h);
            for (let j = 0; j < i; j++) {
                e[j] -= hh * d[j];
            }
            for (let j = 0; j < i; j++) {
                f = d[j];
                g = e[j];
                for (let k = j; k <= i - 1; k++) {
                    V.set(V.at(k, j) - (f * e[k] + g * d[k]), k, j);
                }
                d[j] = V.at(i - 1, j);
                V.set(0.0, i, j);
            }
        }
        d[i] = h;
    }

    // Accumulate transformations.

    for (let i = 0; i < n - 1; i++) {
        V.set(V.at(i, i), n - 1, i);
        V.set(1.0, i, i);
        const h = d[i + 1];
        if (h != 0.0) {
            for (let k = 0; k <= i; k++) {
                d[k] = V.at(k, i + 1) / h;
            }
            for (let j = 0; j <= i; j++) {
                let g = 0.0;
                for (let k = 0; k <= i; k++) {
                    g += V.at(k, i + 1) * V.at(k, j);
                }
                for (let k = 0; k <= i; k++) {
                    V.set(V.at(k, j) - g * d[k], k, j);
                }
            }
        }
        for (let k = 0; k <= i; k++) {
            V.set(0.0, k, i + 1);
        }
    }
    for (let j = 0; j < n; j++) {
        d[j] = V.at(n - 1, j);
        V.set(0.0, n - 1, j);
    }
    V.set(1.0, n - 1, n - 1);
    e[0] = 0.0;
}
//*****************************************************
/**
 * Adapted from the public domain JAMA library (JAMA)
 * @param {Number} xr Real part of x
 * @param {Number} xi Imaginary part of x
 * @param {Number} yr Real part of y
 * @param {Number} yi Imaginary part of y
 * @returns {Array} [r,i] The real and imaginary part of the division
 */
function cdiv(xr, xi, yr, yi) {
    let r, d;
    let cdivr, cdivi;
    if (Math.abs(yr) > Math.abs(yi)) {
        r = yi / yr;
        d = yr + r * yi;
        cdivr = (xr + r * xi) / d;
        cdivi = (xi - r * xr) / d;
    } else {
        r = yr / yi;
        d = yi + r * yr;
        cdivr = (r * xr + xi) / d;
        cdivi = (r * xi - xr) / d;
    }

    return [cdivr, cdivi];
}
//*****************************************************
/**
 * Adapted from the public domain JAMA library (JAMA)
 * @param {AbstractMat} V
 * @param {AbstractMat} H
 * @param {AbstractMat} d
 * @param {AbstractMat} e
 */
function hqr2(V, H, d, e) {
    const nn = V.rows();
    let n = nn - 1;
    const low = 0;
    const high = nn - 1;
    const eps = Math.pow(2.0, -52.0);
    let exshift = 0.0;
    let p = 0, q = 0, r = 0, s = 0, z = 0, t, w, x, y;

    // Store roots isolated by balanc and compute matrix norm

    let norm = 0.0;
    for (let i = 0; i < nn; i++) {
        if (i < low | i > high) {
            d[i] = H.at(i, i);
            e[i] = 0.0;
        }
        for (let j = Math.max(i - 1, 0); j < nn; j++) {
            norm = norm + Math.abs(H.at(i, j));
        }
    }

    // Outer loop over eigenvalue index

    let iter = 0;
    while (n >= low) {

        // Look for single small sub-diagonal element

        let l = n;
        while (l > low) {
            s = Math.abs(H.at(l - 1, l - 1)) + Math.abs(H.at(l, l));
            if (s === 0.0) {
                s = norm;
            }
            if (Math.abs(H.at(l, l - 1)) < eps * s) {
                break;
            }
            l--;
        }

        // Check for convergence
        // One root found

        if (l === n) {
            H.set(H.at(n, n) + exshift, n, n);
            d[n] = H.at(n, n);
            e[n] = 0.0;
            n--;
            iter = 0;

            // Two roots found

        } else if (l === n - 1) {
            w = H.at(n, n - 1) * H.at(n - 1, n);
            p = (H.at(n - 1, n - 1) - H.at(n, n)) / 2.0;
            q = p * p + w;
            z = Math.sqrt(Math.abs(q));
            H.set(H.at(n, n) + exshift, n, n);
            H.set(H.at(n - 1, n - 1) + exshift, n - 1, n - 1);
            x = H.at(n, n);

            // Real pair

            if (q >= 0) {
                if (p >= 0) {
                    z = p + z;
                } else {
                    z = p - z;
                }
                d[n - 1] = x + z;
                d[n] = d[n - 1];
                if (z !== 0.0) {
                    d[n] = x - w / z;
                }
                e[n - 1] = 0.0;
                e[n] = 0.0;
                x = H.at(n, n - 1);
                s = Math.abs(x) + Math.abs(z);
                p = x / s;
                q = z / s;
                r = Math.sqrt(p * p + q * q);
                p = p / r;
                q = q / r;

                // Row modification

                for (let j = n - 1; j < nn; j++) {
                    z = H.at(n - 1, j);
                    H.set(q * z + p * H.at(n, j), n - 1, j);
                    H.set(q * H.at(n, j) - p * z, n, j);
                }

                // Column modification

                for (let i = 0; i <= n; i++) {
                    z = H.at(i, n - 1);
                    H.set(q * z + p * H.at(i, n), i, n - 1);
                    H.set(q * H.at(i, n) - p * z, i, n);
                }

                // Accumulate transformations

                for (let i = low; i <= high; i++) {
                    z = V.at(i, n - 1);
                    V.set(q * z + p * V.at(i, n), i, n - 1);
                    V.set(q * V.at(i, n) - p * z, i, n);
                }

                // Complex pair

            } else {
                d[n - 1] = x + p;
                d[n] = x + p;
                e[n - 1] = z;
                e[n] = -z;
            }
            n = n - 2;
            iter = 0;

            // No convergence yet

        } else {

            // Form shift

            x = H.at(n, n);
            y = 0.0;
            w = 0.0;
            if (l < n) {
                y = H.at(n - 1, n - 1);
                w = H.at(n, n - 1) * H.at(n - 1, n);
            }

            // Wilkinson's original ad hoc shift

            if (iter === 10) {
                exshift += x;
                for (let i = low; i <= n; i++) {
                    H.set(H.at(i, i) - x, i, i);
                }
                s = Math.abs(H.at(n, n - 1)) + Math.abs(H.at(n - 1, n - 2));
                x = y = 0.75 * s;
                w = -0.4375 * s * s;
            }

            // MATLAB's new ad hoc shift

            if (iter === 30) {
                s = (y - x) / 2.0;
                s = s * s + w;
                if (s > 0) {
                    s = Math.sqrt(s);
                    if (y < x) {
                        s = -s;
                    }
                    s = x - w / ((y - x) / 2.0 + s);
                    for (let i = low; i <= n; i++) {
                        H.set(H.at(i, i) - s, i, i);
                    }
                    exshift += s;
                    x = y = w = 0.964;
                }
            }

            iter = iter + 1;   // (Could check iteration count here.)

            // Look for two consecutive small sub-diagonal elements

            let m = n - 2;
            while (m >= l) {
                z = H.at(m, m);
                r = x - z;
                s = y - z;
                p = (r * s - w) / H.at(m + 1, m) + H.at(m, m + 1);
                q = H.at(m + 1, m + 1) - z - r - s;
                r = H.at(m + 2, m + 1);
                s = Math.abs(p) + Math.abs(q) + Math.abs(r);
                p = p / s;
                q = q / s;
                r = r / s;
                if (m === l) {
                    break;
                }
                if (Math.abs(H.at(m, m - 1)) * (Math.abs(q) + Math.abs(r)) <
                    eps * (Math.abs(p) * (Math.abs(H.at(m - 1, m - 1)) + Math.abs(z) +
                        Math.abs(H.at(m + 1, m + 1))))) {
                    break;
                }
                m--;
            }

            for (let i = m + 2; i <= n; i++) {
                H.set(0.0, i, i - 2);
                if (i > m + 2) {
                    H.set(0.0, i, i - 3);
                }
            }

            // Double QR step involving rows l:n and columns m:n


            for (let k = m; k <= n - 1; k++) {
                const notlast = (k !== n - 1);
                if (k != m) {
                    p = H.at(k, k - 1);
                    q = H.at(k + 1, k - 1);
                    r = (notlast ? H.at(k + 2, k - 1) : 0.0);
                    x = Math.abs(p) + Math.abs(q) + Math.abs(r);
                    if (x === 0.0) {
                        continue;
                    }
                    p = p / x;
                    q = q / x;
                    r = r / x;
                }

                s = Math.sqrt(p * p + q * q + r * r);
                if (p < 0) {
                    s = -s;
                }
                if (s != 0) {
                    if (k != m) {
                        H.set(-s * x, k, k - 1);
                    } else if (l != m) {
                        H.set(-H.at(k, k - 1), k, k - 1);
                    }
                    p = p + s;
                    x = p / s;
                    y = q / s;
                    z = r / s;
                    q = q / p;
                    r = r / p;

                    // Row modification

                    for (let j = k; j < nn; j++) {
                        p = H.at(k, j) + q * H.at(k + 1, j);
                        if (notlast) {
                            p = p + r * H.at(k + 2, j);
                            H.set(H.at(k + 2, j) - p * z, k + 2, j);
                        }
                        H.set(H.at(k, j) - p * x, k, j);
                        H.set(H.at(k + 1, j) - p * y, k + 1, j);
                    }

                    // Column modification

                    for (let i = 0; i <= Math.min(n, k + 3); i++) {
                        p = x * H.at(i, k) + y * H.at(i, k + 1);
                        if (notlast) {
                            p = p + z * H.at(i, k + 2);
                            H.set(H.at(i, k + 2) - p * r, i, k + 2);
                        }
                        H.set(H.at(i, k) - p, i, k);
                        H.set(H.at(i, k + 1) - p * q, i, k + 1);
                    }

                    // Accumulate transformations

                    for (let i = low; i <= high; i++) {
                        p = x * V.at(i, k) + y * V.at(i, k + 1);
                        if (notlast) {
                            p = p + z * V.at(i, k + 2);
                            V.set(V.at(i, k + 2) - p * r, i, k + 2);
                        }
                        V.set(V.at(i, k) - p, i, k);
                        V.set(V.at(i, k + 1) - p * q, i, k + 1);
                    }
                }  // (s != 0)
            }  // k loop
        }  // check convergence
    }  // while (n >= low)

    // Backsubstitute to find vectors of upper triangular form

    if (norm === 0.0) {
        return;
    }

    for (n = nn - 1; n >= 0; n--) {
        p = d[n];
        q = e[n];

        // Real vector

        if (q === 0) {
            let l = n;
            H.set(1.0, n, n);
            for (let i = n - 1; i >= 0; i--) {
                w = H.at(i, i) - p;
                r = 0.0;
                for (let j = l; j <= n; j++) {
                    r = r + H.at(i, j) * H.at(j, n);
                }
                if (e[i] < 0.0) {
                    z = w;
                    s = r;
                } else {
                    l = i;
                    if (e[i] === 0.0) {
                        if (w !== 0.0) {
                            H.set(-r / w, i, n);
                        } else {
                            H.set(-r / (eps * norm), i, n);
                        }

                        // Solve real equations

                    } else {
                        x = H.at(i, i + 1);
                        y = H.at(i + 1, i);
                        q = (d[i] - p) * (d[i] - p) + e[i] * e[i];
                        t = (x * s - z * r) / q;
                        H.set(t, i, n);
                        if (Math.abs(x) > Math.abs(z)) {
                            H.set((-r - w * t) / x, i + 1, n);
                        } else {
                            H.set((-s - y * t) / z, i + 1, n);
                        }
                    }

                    // Overflow control

                    t = Math.abs(H.at(i, n));
                    if ((eps * t) * t > 1) {
                        for (let j = i; j <= n; j++) {
                            H.set(H[j][n] / t, j, n);
                        }
                    }
                }
            }

            // Complex vector

        } else if (q < 0) {
            let l = n - 1;

            // Last vector component imaginary so matrix is triangular

            if (Math.abs(H.at(n, n - 1)) > Math.abs(H.at(n - 1, n))) {
                H.set(q / H.at(n, n - 1), n - 1, n - 1);
                H.set(-(H.at(n, n) - p) / H.at(n, n - 1), n - 1, n);
            } else {
                const [cdivr, cdivi] = cdiv(0.0, -H.at(n - 1, n), H.at(n - 1, n - 1) - p, q);
                H.set(cdivr, n - 1, n - 1);
                H.set(cdivi, n - 1, n);
            }
            H.set(0.0, n, n - 1);
            H.set(1.0, n, n);
            for (let i = n - 2; i >= 0; i--) {
                let ra, sa, vr, vi;
                ra = 0.0;
                sa = 0.0;
                for (let j = l; j <= n; j++) {
                    ra = ra + H.at(i, j) * H.at(j, n - 1);
                    sa = sa + H.at(i, j) * H.at(j, n);
                }
                w = H.at(i, i) - p;

                if (e[i] < 0.0) {
                    z = w;
                    r = ra;
                    s = sa;
                } else {
                    l = i;
                    if (e[i] === 0) {
                        const [cdivr, cdivi] = cdiv(-ra, -sa, w, q);
                        H.set(cdivr, i, n - 1);
                        H.set(cdivi, i, n);
                    } else {

                        // Solve complex equations

                        x = H.at(i, i + 1);
                        y = H.at(i + 1, i);
                        vr = (d[i] - p) * (d[i] - p) + e[i] * e[i] - q * q;
                        vi = (d[i] - p) * 2.0 * q;
                        if (vr === 0.0 & vi === 0.0) {
                            vr = eps * norm * (Math.abs(w) + Math.abs(q) +
                                Math.abs(x) + Math.abs(y) + Math.abs(z));
                        }
                        const [cdivr, cdivi] = cdiv(x * r - z * ra + q * sa, x * s - z * sa - q * ra, vr, vi);
                        H.set(cdivr, i, n - 1);
                        H.set(cdivi, i, n);
                        if (Math.abs(x) > (Math.abs(z) + Math.abs(q))) {
                            H.set((-ra - w * H.at(i, n - 1) + q * H.at(i, n)) / x, i + 1, n - 1);
                            H.set((-sa - w * H.at(i, n) - q * H.at(i, n - 1)) / x, i + 1, n);
                        } else {
                            const [cdivr, cdivi] = cdiv(-r - y * H[i][n - 1], -s - y * H[i][n], z, q);
                            H.set(cdivr, i + 1, n - 1);
                            H.set(cdivi, i + 1, n);
                        }
                    }

                    // Overflow control

                    t = Math.max(Math.abs(H.at(i, n - 1)), Math.abs(H.at(i, n)));
                    if ((eps * t) * t > 1) {
                        for (let j = i; j <= n; j++) {
                            H.set(H.at(j, n - 1) / t, j, n - 1);
                            H.set(H.at(j, n) / t, j, n);
                        }
                    }
                }
            }
        }
    }

    // Vectors of isolated roots

    for (let i = 0; i < nn; i++) {
        if (i < low | i > high) {
            for (let j = i; j < nn; j++) {
                V.set(H.at(i, j), i, j);
            }
        }
    }

    // Back transformation to get eigenvectors of original matrix

    for (let j = nn - 1; j >= low; j--) {
        for (let i = low; i <= high; i++) {
            z = 0.0;
            for (let k = low; k <= Math.min(j, high); k++) {
                z = z + V.at(i, k) * H.at(k, j);
            }
            V.set(z, i, j);
        }
    }
}
//*****************************************************
/**
 * Adapted from the public domain JAMA library (JAMA)
 * @param {AbstractMat} V
 * @param {AbstractMat} d
 * @param {AbstractMat} e
 */
function tq12(V, d, e) {
    const n = V.rows();

    for (let i = 1; i < n; i++) {
        e[i - 1] = e[i];
    }
    e[n - 1] = 0.0;

    let f = 0.0;
    let tst1 = 0.0;
    let eps = Math.pow(2.0, -52.0);
    for (let l = 0; l < n; l++) {

        // Find small subdiagonal element

        tst1 = Math.max(tst1, Math.abs(d[l]) + Math.abs(e[l]));
        let m = l;
        while (m < n) {
            if (Math.abs(e[m]) <= eps * tst1) {
                break;
            }
            m++;
        }

        // If m == l, d[l] is an eigenvalue,
        // otherwise, iterate.

        if (m > l) {
            let iter = 0;
            do {
                iter = iter + 1;  // (Could check iteration count here.)

                // Compute implicit shift

                let g = d[l];
                let p = (d[l + 1] - g) / (2.0 * e[l]);
                let r = hypot(p, 1.0);
                if (p < 0) {
                    r = -r;
                }
                d[l] = e[l] / (p + r);
                d[l + 1] = e[l] * (p + r);
                const dl1 = d[l + 1];
                let h = g - d[l];
                for (let i = l + 2; i < n; i++) {
                    d[i] -= h;
                }
                f = f + h;

                // Implicit QL transformation.

                p = d[m];
                let c = 1.0;
                let c2 = c;
                let c3 = c;
                const el1 = e[l + 1];
                let s = 0.0;
                let s2 = 0.0;
                for (let i = m - 1; i >= l; i--) {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = hypot(p, e[i]);
                    e[i + 1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i + 1] = h + s * (c * g + s * d[i]);

                    // Accumulate transformation.

                    for (let k = 0; k < n; k++) {
                        h = V.at(k, i + 1);
                        V.set(s * V.at(k, i) + c * h, k, i + 1);
                        V.set(c * V.at(k, i) - s * h, k, i);
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;

                // Check for convergence.

            } while (Math.abs(e[l]) > eps * tst1);
        }
        d[l] = d[l] + f;
        e[l] = 0.0;
    }

    // Sort eigenvalues and corresponding vectors.

    for (let i = 0; i < n - 1; i++) {
        let k = i;
        let p = d[i];
        for (let j = i + 1; j < n; j++) {
            if (d[j] < p) {
                k = j;
                p = d[j];
            }
        }
        if (k != i) {
            d[k] = d[i];
            d[i] = p;
            for (let j = 0; j < n; j++) {
                p = V.at(j, i);
                V.set(V.at(j, k), j, i);
                V.set(p, j, k);
            }
        }
    }
}

//*****************************************************
/**
 * Adapted from the public domain JAMA library (JAMA)
 * @param {AbstractMat} V
 * @param {AbstractMat} H
 */
function orthes(V, H) {
    const n = V.rows();
    const low = 0;
    const high = n - 1;
    const type = V.type();
    const ort = new type(n);


    for (let m = low + 1; m <= high - 1; m++) {

        // Scale column.

        let scale = 0.0;
        for (let i = m; i <= high; i++) {
            scale = scale + Math.abs(H.at(i, m - 1));
        }
        if (scale !== 0.0) {

            // Compute Householder transformation.

            let h = 0.0;
            for (let i = high; i >= m; i--) {
                ort[i] = H.at(i, m - 1) / scale;
                h += ort[i] * ort[i];
            }
            let g = Math.sqrt(h);
            if (ort[m] > 0) {
                g = -g;
            }
            h = h - ort[m] * g;
            ort[m] = ort[m] - g;

            // Apply Householder similarity transformation
            // H = (I-u*u'/h)*H*(I-u*u')/h)

            for (let j = m; j < n; j++) {
                let f = 0.0;
                for (let i = high; i >= m; i--) {
                    f += ort[i] * H.at(i, j);
                }
                f = f / h;
                for (let i = m; i <= high; i++) {
                    H.set(H.at(i, j) - f * ort[i], i, j);
                }
            }

            for (let i = 0; i <= high; i++) {
                let f = 0.0;
                for (let j = high; j >= m; j--) {
                    f += ort[j] * H.at(i, j);
                }
                f = f / h;
                for (let j = m; j <= high; j++) {
                    H.set(H.at(i, j) - f * ort[j], i, j);
                }
            }
            ort[m] = scale * ort[m];
            H.set(scale * g, m, m - 1);
        }
    }

    // Accumulate transformations (Algol's ortran).

    for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
            V.set(i === j ? 1.0 : 0.0, i, j);
        }
    }

    for (let m = high - 1; m >= low + 1; m--) {
        if (H.at(m, m - 1) !== 0.0) {
            for (let i = m + 1; i <= high; i++) {
                ort[i] = H.at(i, m - 1);
            }
            for (let j = m; j <= high; j++) {
                let g = 0.0;
                for (let i = m; i <= high; i++) {
                    g += ort[i] * V.at(i, j);
                }
                // Double division avoids possible underflow
                g = (g / ort[m]) / H.at(m, m - 1);
                for (let i = m; i <= high; i++) {
                    V.set(V.at(i, j) + g * ort[i], i, j);
                }
            }
        }
    }
}

//*****************************************************
/**
 * Stores the Eigen Decomposition of a matrix: A*V = V*D
 * For a symmetric matrix, the eigenvalues are all real and D is orthogonal.
 * For a non-symmetric matrix, V is block-diagonal and V might be singular.
 */
class Eigen {

    /**
     * Construct a new Eigen object
     * @param {AbstractMat} V The Eigenvector matrix
     * @param {Array | TypedArray} d The real part of the eigenvalues
     * @param {Array | TypedArray} e The imaginary part of the eigenvalues
     */
    constructor(V, d, e) {
        this.V = V;
        this.d = d;
        this.e = e;
    }

    /**
     * Construct a new Eigen object
     * @param {AbstractMat} V The Eigenvector matrix
     * @param {Array | TypedArray} d The real part of the eigenvalues
     * @param {Array | TypedArray} e The imaginary part of the eigenvalues
     * @returns {Eigen} A new Eigen object
     */
    static new(V, d, e) {
        return new Eigen(V, d, e);
    }

    /**
     * Computes the eigenvalue decomposition of a matrix
     * @param {AbstractMat} a The input matrix
     * @returns {Eigen} The eigenvalue decomposition
     */
    static compute(a) {
        return computeEigen(a);
    }

    /**
     * @returns {AbstractMat} The eigenvector matrix
     */
    getV() {
        return this.V;
    }

    /**
     * The real part of the eigenvalues: real(lambda_i) for i=1:n
     * @returns {Array | TypedArray} The real eigenvalues
     */
    getRealEigenvalues() {
        return this.d;
    }
    /**
      * The imaginary part of the eigenvalues: imag(lambda_i) for i=1:n
      * @returns {Array | TypedArray} The imaginary eigenvalues
      */
    getImaginaryEigenvalues() {
        return this.e;
    }

    /**
     * Returns the eigenvalue matrix
     * For symmetric matrices, this will be diagonal, otherwise block-diagonal
     * @returns {AbstractMat} The eigenvalue matrix
     */
    getD() {
        const D = setZero(similar(this.V));
        const d = this.d;
        const e = this.e;
        const n = d.length;
        insert(diag(D), VecF32.from(d));
        for (let i = 0; i < n; i++) {
            if (e[i] > 0) {
                D.set(e[i], i, i + 1);
            } else if (e[i] < 0) {
                D.set(e[i], i, i - 1);
            }
        }
        return D;
    }
}
//*****************************************************

/**
 * Compute Eigenvalue decomposition of a matrix
 * This function will adapt to symmetric and non-symmetric matrices.
 * For more information about the result @see Eigen
 * Adapted from the public domain JAMA library (JAMA)
 * @returns {Eigen} The eigenvalue decomposition
 */
function computeEigen(a) {
    const n = a.cols();
    if (n !== a.rows()) {
        throw new Error("Eigenvalues only defined for square matrices");
    }

    const type = a.type();
    const d = new type(n);
    const e = new type(n);

    // check for symmetry
    if (isSymmetric(a)) {
        const V = copy(a);

        // tridiagonalize
        tred2(V, d, e);
        // diagonalize
        tq12(V, d, e);
        return Eigen.new(V, d, e);
    }
    else {
        const V = similar(a);
        const H = copy(a);

        orthes(V, H);

        hqr2(V, H, d, e);
        return Eigen.new(V, d, e);
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
     * Computes the singular value decomposition of a general matrix.
     *
     * A MxN matrix A can be decomposed into U * S * V^T
     *
     * For a MxN matrix with M >= N:
     * U is a unitary MxN matrix
     * S is a NxN diagonal matrix
     * V is a unitary NxN matrix
     * 
     * If M < N:
     * U is a unitary NxN matrix
     * S is a MxN diagonal matrix
     * V is a unitary MxN matrix
     *
     * The diagonal entries of S are the singular values of A and the columns of U and V are the corresponding left and right singular vectors.
     * The singular values and their singular vectors are sorted in ascending order. S is represented as a vector containing the singular values.
     *
     * Adapted from the JAMA library
     * @param {AbstractMat} m - The input matrix
     * @returns {SVD} The SVD
     */
    static compute(m) {
        return computeSVD(m);
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
        const U = this.U;
        const Vt = this.Vt;
        const S = this.S;

        const M = U.rows();
        // since V will be 
        const N = Vt.cols();
        out = out !== undefined ? out : similar(U, M, N);

        if (out.rows() !== M || out.cols() !== N) {
            throw new Error("Output matrix has wrong dimensions");
        }

        mult(U, mult(Diagonal.new(toArray(S), U.cols(), Vt.rows()), Vt), out);

        return out;

    }

    /**
     * Creates the pseudo inverse of this SVD
     * If the matrix has full rank, this will equal the inverse
     *
     * @returns {SVD} The SVD of the inverse
     */
    inv() {
        const U = copy(this.U);
        const V = copy(this.V);
        const S = copy(this.S);

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
 * For a MxN matrix with M >= N:
 * U is a unitary MxN matrix
 * S is a NxN diagonal matrix
 * V is a unitary NxN matrix
 * 
 * If M < N:
 * U is a unitary NxN matrix
 * S is a MxN diagonal matrix
 * V is a unitary MxN matrix
 *
 * The diagonal entries of S are the singular values of A and the columns of U and V are the corresponding left and right singular vectors.
 * The singular values and their singular vectors are sorted in ascending order. S is represented as a vector containing the singular values.
 *
 * Adapted from the JAMA library
 * @param {AbstractMat} m - The input matrix
 * @returns {SVD} The SVD
 */
function computeSVD(A) {

    // in case of m < n -> transpose temporarily
    let isTransposed = false;

    if (A.rows() < A.cols()) {
        A = transpose(A);
        isTransposed = true;
    }

    A = copy(A);

    const m = A.rows();
    const n = A.cols();

    const nu = Math.min(m, n);
    const type = A.type();
    const s = new type(Math.min(m + 1, n));
    const U = similar(A, m, nu);
    const V = similar(A, n, n);
    const e = new type(n);
    const work = new type(m);
    const wantu = true;
    const wantv = true;

    // Reduce A to bidiagonal form, storing the diagonal elements
    // in s and the super-diagonal elements in e.

    const nct = Math.min(m - 1, n);
    const nrt = Math.max(0, Math.min(n - 2, m));
    for (let k = 0; k < Math.max(nct, nrt); k++) {
        if (k < nct) {

            // Compute the transformation for the k-th column and
            // place the k-th diagonal in s[k].
            // Compute 2-norm of k-th column without under/overflow.
            s[k] = 0;
            for (let i = k; i < m; i++) {
                s[k] = hypot(s[k], A.at(i, k));
            }
            if (s[k] != 0.0) {
                if (A.at(k, k) < 0.0) {
                    s[k] = -s[k];
                }
                for (let i = k; i < m; i++) {
                    A.set(A.at(i, k) / s[k], i, k);
                }
                A.set(A.at(k, k) + 1.0, k, k);
            }
            s[k] = -s[k];
        }
        for (let j = k + 1; j < n; j++) {
            if ((k < nct) & (s[k] != 0.0)) {

                // Apply the transformation.

                let t = 0;
                for (let i = k; i < m; i++) {
                    t += A.at(i, k) * A.at(i, j);
                }
                t = -t / A.at(k, k);
                for (let i = k; i < m; i++) {
                    A.set(A.at(i, j) + t * A.at(i, k), i, j);
                }
            }

            // Place the k-th row of A into e for the
            // subsequent calculation of the row transformation.

            e[j] = A.at(k, j);
        }
        if (wantu & (k < nct)) {

            // Place the transformation in U for subsequent back
            // multiplication.

            for (let i = k; i < m; i++) {
                U.set(A.at(i, k), i, k);
            }
        }
        if (k < nrt) {

            // Compute the k-th row transformation and place the
            // k-th super-diagonal in e[k].
            // Compute 2-norm without under/overflow.
            e[k] = 0;
            for (let i = k + 1; i < n; i++) {
                e[k] = hypot(e[k], e[i]);
            }
            if (e[k] != 0.0) {
                if (e[k + 1] < 0.0) {
                    e[k] = -e[k];
                }
                for (let i = k + 1; i < n; i++) {
                    e[i] /= e[k];
                }
                e[k + 1] += 1.0;
            }
            e[k] = -e[k];
            if ((k + 1 < m) & (e[k] != 0.0)) {

                // Apply the transformation.

                for (let i = k + 1; i < m; i++) {
                    work[i] = 0.0;
                }
                for (let j = k + 1; j < n; j++) {
                    for (let i = k + 1; i < m; i++) {
                        work[i] += e[j] * A.at(i, j);
                    }
                }
                for (let j = k + 1; j < n; j++) {
                    const t = -e[j] / e[k + 1];
                    for (let i = k + 1; i < m; i++) {
                        A.set(A.at(i, j) + t * work[i], i, j);
                    }
                }
            }
            if (wantv) {

                // Place the transformation in V for subsequent
                // back multiplication.

                for (let i = k + 1; i < n; i++) {
                    V.set(e[i], i, k);
                }
            }
        }
    }

    // Set up the final bidiagonal matrix or order p.

    let p = Math.min(n, m + 1);
    if (nct < n) {
        s[nct] = A.at(nct, nct);
    }
    if (m < p) {
        s[p - 1] = 0.0;
    }
    if (nrt + 1 < p) {
        e[nrt] = A.at(nrt, p - 1);
    }
    e[p - 1] = 0.0;

    // If required, generate U.

    if (wantu) {
        for (let j = nct; j < nu; j++) {
            for (let i = 0; i < m; i++) {
                U.set(0.0, i, j);
            }
            U.set(1.0, j, j);
        }
        for (let k = nct - 1; k >= 0; k--) {
            if (s[k] != 0.0) {
                for (let j = k + 1; j < nu; j++) {
                    let t = 0;
                    for (let i = k; i < m; i++) {
                        t += U.at(i, k) * U.at(i, j);
                    }
                    t = -t / U.at(k, k);
                    for (let i = k; i < m; i++) {
                        U.set(U.at(i, j) + t * U.at(i, k), i, j);
                    }
                }
                for (let i = k; i < m; i++) {
                    U.set(-U.at(i, k), i, k);
                }
                U.set(1.0 + U.at(k, k), k, k);
                for (let i = 0; i < k - 1; i++) {
                    U.set(0.0, i, k);
                }
            } else {
                for (let i = 0; i < m; i++) {
                    U.set(0.0, i, k);
                }
                U.set(1.0, k, k);
            }
        }
    }

    // If required, generate V.

    if (wantv) {
        for (let k = n - 1; k >= 0; k--) {
            if ((k < nrt) & (e[k] != 0.0)) {
                for (let j = k + 1; j < nu; j++) {
                    let t = 0;
                    for (let i = k + 1; i < n; i++) {
                        t += V.at(i, k) * V.at(i, j);
                    }
                    t = -t / V.at(k + 1, k);
                    for (let i = k + 1; i < n; i++) {
                        V.set(V.at(i, j) + t * V.at(i, k), i, j);
                    }
                }
            }
            for (let i = 0; i < n; i++) {
                V.set(0.0, i, k);
            }
            V.set(1.0, k, k);
        }
    }

    // Main iteration loop for the singular values.

    const pp = p - 1;
    let iter = 0;
    const eps = Math.pow(2.0, -52.0);
    const tiny = Math.pow(2.0, -966.0);
    while (p > 0) {
        let k, kase;

        // Here is where a test for too many iterations would go.

        // This section of the program inspects for
        // negligible elements in the s and e arrays.  On
        // completion the variables kase and k are set as follows.

        // kase = 1     if s(p) and e[k-1] are negligible and k<p
        // kase = 2     if s(k) is negligible and k<p
        // kase = 3     if e[k-1] is negligible, k<p, and
        //              s(k), ..., s(p) are not negligible (qr step).
        // kase = 4     if e(p-1) is negligible (convergence).

        for (k = p - 2; k >= -1; k--) {
            if (k === -1) {
                break;
            }
            if (Math.abs(e[k]) <=
                tiny + eps * (Math.abs(s[k]) + Math.abs(s[k + 1]))) {
                e[k] = 0.0;
                break;
            }
        }
        if (k === p - 2) {
            kase = 4;
        } else {
            let ks;
            for (ks = p - 1; ks >= k; ks--) {
                if (ks === k) {
                    break;
                }
                const t = (ks != p ? Math.abs(e[ks]) : 0.) +
                    (ks != k + 1 ? Math.abs(e[ks - 1]) : 0.);
                if (Math.abs(s[ks]) <= tiny + eps * t) {
                    s[ks] = 0.0;
                    break;
                }
            }
            if (ks === k) {
                kase = 3;
            } else if (ks === p - 1) {
                kase = 1;
            } else {
                kase = 2;
                k = ks;
            }
        }
        k++;

        // Perform the task indicated by kase.

        switch (kase) {

            // Deflate negligible s(p).

            case 1: {
                let f = e[p - 2];
                e[p - 2] = 0.0;
                for (let j = p - 2; j >= k; j--) {
                    let t = hypot(s[j], f);
                    const cs = s[j] / t;
                    const sn = f / t;
                    s[j] = t;
                    if (j != k) {
                        f = -sn * e[j - 1];
                        e[j - 1] = cs * e[j - 1];
                    }
                    if (wantv) {
                        for (let i = 0; i < n; i++) {
                            t = cs * V.at(i, j) + sn * V.at(i, p - 1);
                            V.set(-sn * V.at(i, j) + cs * V.at(i, p - 1), i, p - 1);
                            V.set(t, i, j);
                        }
                    }
                }
            }
                break;

            // Split at negligible s(k).

            case 2: {
                let f = e[k - 1];
                e[k - 1] = 0.0;
                for (let j = k; j < p; j++) {
                    let t = hypot(s[j], f);
                    const cs = s[j] / t;
                    const sn = f / t;
                    s[j] = t;
                    f = -sn * e[j];
                    e[j] = cs * e[j];
                    if (wantu) {
                        for (let i = 0; i < m; i++) {
                            t = cs * U.at(i, j) + sn * U.at(i, k - 1);
                            U.set(-sn * U.at(i, j) + cs * U.at(i, k - 1), i, k - 1);
                            U.set(t, i, j);
                        }
                    }
                }
            }
                break;

            // Perform one qr step.

            case 3: {

                // Calculate the shift.

                const scale = Math.max(Math.max(Math.max(Math.max(
                    Math.abs(s[p - 1]), Math.abs(s[p - 2])), Math.abs(e[p - 2])),
                    Math.abs(s[k])), Math.abs(e[k]));
                const sp = s[p - 1] / scale;
                const spm1 = s[p - 2] / scale;
                const epm1 = e[p - 2] / scale;
                const sk = s[k] / scale;
                const ek = e[k] / scale;
                const b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2.0;
                const c = (sp * epm1) * (sp * epm1);
                let shift = 0.0;
                if ((b != 0.0) | (c != 0.0)) {
                    shift = Math.sqrt(b * b + c);
                    if (b < 0.0) {
                        shift = -shift;
                    }
                    shift = c / (b + shift);
                }
                let f = (sk + sp) * (sk - sp) + shift;
                let g = sk * ek;

                // Chase zeros.

                for (let j = k; j < p - 1; j++) {
                    let t = hypot(f, g);
                    let cs = f / t;
                    let sn = g / t;
                    if (j != k) {
                        e[j - 1] = t;
                    }
                    f = cs * s[j] + sn * e[j];
                    e[j] = cs * e[j] - sn * s[j];
                    g = sn * s[j + 1];
                    s[j + 1] = cs * s[j + 1];
                    if (wantv) {
                        for (let i = 0; i < n; i++) {
                            t = cs * V.at(i, j) + sn * V.at(i, j + 1);
                            V.set(-sn * V.at(i, j) + cs * V.at(i, j + 1), i, j + 1);
                            V.set(t, i, j);
                        }
                    }
                    t = hypot(f, g);
                    cs = f / t;
                    sn = g / t;
                    s[j] = t;
                    f = cs * e[j] + sn * s[j + 1];
                    s[j + 1] = -sn * e[j] + cs * s[j + 1];
                    g = sn * e[j + 1];
                    e[j + 1] = cs * e[j + 1];
                    if (wantu && (j < m - 1)) {
                        for (let i = 0; i < m; i++) {
                            t = cs * U.at(i, j) + sn * U.at(i, j + 1);
                            U.set(-sn * U.at(i, j) + cs * U.at(i, j + 1), i, j + 1);
                            U.set(t, i, j);
                        }
                    }
                }
                e[p - 2] = f;
                iter = iter + 1;
            }
                break;

            // Convergence.

            case 4: {

                // Make the singular values positive.

                if (s[k] <= 0.0) {
                    s[k] = (s[k] < 0.0 ? -s[k] : 0.0);
                    if (wantv) {
                        for (let i = 0; i <= pp; i++) {
                            V.set(-V.at(i, k), i, k);
                        }
                    }
                }

                // Order the singular values.

                while (k < pp) {
                    if (s[k] >= s[k + 1]) {
                        break;
                    }
                    let t = s[k];
                    s[k] = s[k + 1];
                    s[k + 1] = t;
                    if (wantv && (k < n - 1)) {
                        for (let i = 0; i < n; i++) {
                            t = V.at(i, k + 1);
                            V.set(V.at(i, k), i, k + 1);
                            V.set(t, i, k);
                        }
                    }
                    if (wantu && (k < m - 1)) {
                        for (let i = 0; i < m; i++) {
                            t = U.at(i, k + 1);
                            U.set(U.at(i, k), i, k + 1);
                            U.set(t, i, k);
                        }
                    }
                    k++;
                }
                iter = 0;
                p--;
            }
                break;
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

        return SVD.new(V, VecF32.from(s), U);

    } else {
        return SVD.new(U, VecF32.from(s), V);
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
    const mstr = map(m, x => x.toString(), MatAny.uninitialized(m.rows(), m.cols()));

    // find maximum length per col
    const ls = colreduce(mstr, x => reduce(x, (acc, v) => Math.max(acc, v.length), 0),
        TypedMatFactory.new(Int16Array).uninitialized(1, m.cols()));

    map(mstr, (x, i, j) => pad(x, ls.at(0, j)), mstr);
    const rows = rowreduce(mstr, x => toArray(x).join(" "), MatAny.uninitialized(m.rows(), 1));

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

function toJSON(a) {
    return {
        rows: a.rows(),
        cols: a.cols(),
        data: toArray(a)
    };
}


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
        throw new Error("Cross product only defined on vectors");
    }
    if (a.rows() !== 3 || b.rows() !== 3) {
        throw new Error("Cross product only defined in 3D");
    }


    const a0 = a.at(0);
    const a1 = a.at(1);
    const a2 = a.at(2);

    const b0 = b.at(0);
    const b1 = b.at(1);
    const b2 = b.at(2);

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
    const z = sub(eye, center);
    normalize(z, z);

    const x = cross(up, z);
    normalize(x, x);

    const y = cross(z, x);

    const V = setId(similar(eye, 4, 4));

    const R = block(V, 0, 0, 3, 3);

    insert(row(R, 0), transpose(x));
    insert(row(R, 1), transpose(y));
    insert(row(R, 2), transpose(z));

    const T = setId(similar(eye, 4, 4));
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
    const P = setId(mat(4, 4, f32));

    P.set(2.0 / (right - left), 0, 0);
    P.set(2.0 / (top - bottom), 1, 1);
    P.set(2.0 / (far - near), 2, 2);

    const c = col(P, 3);
    c.set(-(left + right) / (right - left), 0);
    c.set(-(bottom + top) / (top - bottom), 1);
    c.set(-(far + near) / (far - near), 2);

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
    const p = setId(mat(4, 4, f32));
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
    const P = setId(mat(4, 4, f32));

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

    const xsize = near * Math.tan(fov * 0.5);
    const ysize = xsize / aspect;

    const O = ortho(-xsize, xsize, -ysize, ysize, near, far);

    const rtol = setId(mat(4, 4, f32));
    diag(rtol).set(-1.0, 2);

    const persp = centralPerspective(near, far);

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
    const V = setId(mat(4, 4, f32));

    const d = diag(V);
    d.set(w / 2.0, 0);

    if (flipy) {
        d.set(-h / 2.0, 1);

    } else {
        d.set(h / 2.0, 1);

    }
    d.set(0.5, 2);

    const c = col(V, 3);

    c.set(w / 2 + x0, 0);
    if (flipy) {
        c.set(h / 2 - y0, 1);

    } else {

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
 * Creates a new homogeneous matrix from the given one.
 * If the given matrix is NxM, the result will be (N+1)x(M+1)
 * @param mat - The base matrix
 * @param [hcoord] - The homogeneous coordinate. If not given, will be equal to 1
 * @returns {AbstractMat} A homogeneous matrix
 */
function hmat(mat, hcoord = 1) {

    const result = setZero(similar(mat, mat.rows() + 1, mat.cols() + 1));
    insert(block(result, 0, 0, mat.rows(), mat.cols()), mat);
    result.set(hcoord, mat.rows(), mat.cols());

    return result;
}

/**
 * Creates a new homogeneous vector from the given one.
 * If the given vector is N-d, the result will be (N+1)-d
 * @param v - The base vector
 * @param [hcoord] - The homogeneous coordinate. If not given, will be equal to 1
 * @returns {AbstractMat} A homogeneous vector
 */
function hvec(v, hcoord = 1) {
    if (!isVec(v)) {
        throw new Error("Not a vector");
    }
    const result = setZero(similar(v, v.rows() + 1));
    insert(subvec(result, 0, v.rows()), v);
    result.set(hcoord, v.rows());

    return result;
}

/**
 * Computes a 3x3 cross product matrix K from a vector k such that for any vector v: Kv = cross(k,v)
 *
 * @param {AbstractMat} k - A 3D vector
 * @returns {Mat} The cross prodcut matrix
 */
function crossMatrix(k) {
    if (!isVec(k) || k.rows() !== 3) {
        throw new Error("k has to be a 3D vector to compute the cross product matrix");
    }
    const K = setZero(similar(k, 3, 3));

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
    const R = axisAngle(axis, angle);

    const R4 = setId(similar(R, 4, 4));
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
        throw new Error("Translation vector needs to be 3-dimensional");
    }
    const T = setId(similar(t, 4, 4));

    const subcol = subvec(col(T, 3), 0, 3);
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
        throw new Error("Scaling vector needs to be 3-dimensional");
    }
    const S = similar(s, 4, 4);

    const subdiag = subvec(diag(S), 0, 3);
    insert(subdiag, s);
    S.set(1, 3, 3);
    return S;
}


/**
 * Checks whether two matrices are approximately equal
 * @param a - The first matrix
 * @param b - The second matrix
 * @param [eps] - Optional epsilon value, defaults to 1E-7
 * @returns {boolean} True, if both matrices are approximately equal, false otherwise
 */
function approxEqual(a, b, eps = 1E-7) {

    return norm(sub(a, b)) < eps;
}


/**
 * Checks whether a matrix is approximately zero
 * @param a - The  matrix
 * @param [eps] - Optional epsilon value, defaults to 1E-7
 * @returns {boolean} True, if the matrix is approximately zero, false otherwise
 */
function approxZero(a, eps = 1E-7) {
    return norm(a) < eps;
}

/**
 * Constructs a cartesian vector from spherical coordinates
 * @param theta - The polar angle
 * @param phi - The azimuthal angle
 * @param r - The radial distance
 * @param [type] - The output type, will use the default type if not specified
 * @returns {AbstractMat} The cartesian vector
 */
function spherical(theta, phi, r, type = currentDefaultType) {
    const st = Math.sin(theta);

    return TypedVecFactory.new(type).from([r * st * Math.cos(phi), r * st * Math.sin(phi), r * Math.cos(theta)]);
}

/**
 * Computes the spherical representation of a cartesian vector
 * 
 * @param {AbstractMat} p A cartesian vector
 * @returns {{r: Number, theta:Number, phi:Number}} The spherical representation
 */
function cartesianToSpherical(p) {
    const [x, y, z] = [p.at(0), p.at(1), p.at(2)];
    const r = norm(p);
    const phi = Math.atan2(y, x);

    const theta = Math.acos(z / r);
    return { r, theta, phi };
}

/**
 * Permutes the rows of a matrix
 * @param a - The input matrix
 * @param indices - An array specifying the at each index, which other row should be put there. Needs to have the same length as a has rows
 * @param [out] - The output array. Will be created, if not given
 * @returns {AbstractMat}
 */
function permuteRows(a, indices, out) {
    if (indices.length !== a.rows()) {
        throw new Error("Row permutation list length is not equal to rows");
    }

    const tmp = similar(a);

    for (let i = 0; i < tmp.rows(); i++) {
        insert(row(tmp, i), row(a, indices[i]));
    }

    if (out !== undefined) {
        insert(out, tmp);
    } else {
        out = tmp;
    }

    return out;

}

/**
 * Returns the vector from a to b
 * @param a - The first point
 * @param b - The second point
 * @returns {AbstractMat} - The difference vector b-a
 */
function fromTo(a, b) {
    return sub(b, a);
}

/**
 * 3D x vector [1,0,0]
 */
const X = VecF32.from([1, 0, 0]);
/**
 * 3D y vector [0,1,0]
 */
const Y = VecF32.from([0, 1, 0]);
/**
 * 3D z vector [0,0,1]
 */
const Z = VecF32.from([0, 0, 1]);

/**
 * Homogenous 3D x vector [1,0,0,0]
 */
const Xh = VecF32.from([1, 0, 0, 0]);
/**
 * Homogenous 3D y vector [0,1,0,0]
 */
const Yh = VecF32.from([0, 1, 0, 0]);
/**
 * Homogenous 3D z vector [0,0,1,0]
 */
const Zh = VecF32.from([0, 0, 1, 0]);


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
    optional,
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
    fromTo,
    map,
    convert,
    reduce,
    sum,
    neg,
    abs,
    cwiseMult, cwiseDiv,
    cwiseMin, cwiseMax,
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
    swapCol,
    swapRow,
    det,
    inv,
    rank,
    cond,
    computePLUD,
    computeSVD,
    SVD,
    computeEigen,
    Eigen,
    computeUBVD,
    unpackUBV,
    solve,
    solvePLU,
    solveSVD,
    solveLowerTriagonal,
    solveUpperTriagonal,
    prettyprint,
    toString,
    toJSON,
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
    hvec, hmat,
    crossMatrix,
    axisAngle,
    axisAngle4,
    translation,
    scaling,
    approxEqual,
    approxZero,
    spherical,
    cartesianToSpherical,
    permuteRows,
    isSymmetric,
    X, Y, Z,
    Xh, Yh, Zh
};