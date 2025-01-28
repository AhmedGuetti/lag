#ifndef GMATRIX_H_
#define GMATRIX_H_

#include <stdexcept>
#include <vector>
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iostream>

template <class T>
class gmatrix {
public:
    /**
     * Constructs a matrix of size rows x cols and sets
     * all values to 0.
    */
    gmatrix();
    /**
     * Constructs a matrix of size rows x cols and sets
     * all values to 0.
     * @param nRow number of rows
     * @param nCol number of columns
    */
    gmatrix(int nRow, int nCol);
    gmatrix(int nRow, int nCol, const std::vector<T>& data);
    /**
     * Copy constructor.
     * @param mat Matrix from which to copy
    */
    template <class U>
    gmatrix(const gmatrix<U>& matrix);
    gmatrix(const gmatrix<T>& matrix);

    /**
     * Move constructor
     * @param mat
     */
    gmatrix(gmatrix<T>&& mat);

    
    ~gmatrix();

    // Configuration
    /**
     * Resize the matrix dimetion.
     * @param nRow number of the row
     * @param nCol numnber of the column 
     * @return true if we resize with sucess, else otherwise
    */
    bool resize(int nRow, int nCol);
    /**
     * we get the index from 2D --> 1D 
     * @param nRow number of the row
     * @param nCol numnber of the column 
     * @return the Linear index of the element in the position (nRow, nCols)
    */
    int linear_index(int nRow, int nCol) const;
    /**
     * Check if the two matrices given have the same dimetion.
     * @param lhs pointer to the left hand side matrix 
     * @param rhs pointer to the right hand side matrix
     * @return true if lhs and rhs have the same number of rows and column
    */
    inline static bool equalDimension(const gmatrix<T>& lhs, const gmatrix<T>& rhs);

    /**
     * @param nRow number of the row
     * @param nCol numnber of the column 
     * @return return the element in the position (nRow, nCol)
    */
    T get(int nRow, int nCol) const;

    /**
     * @return the number of rows in the matrix
    */
    int rows() const;

    /**
     * @return the number of columns in the matrix
    */
    int cols() const;

    /** the size of a matrix is the number of element inside of it
     * @return the size of the matrix
    */
    int size() const;


    /** Data is the matrix holding the acctual serialised version of the array
     * @return the data array of a matrix 
    */
    inline T*       data();
    
    /** Data is the matrix holding the acctual serialised version of the array
     * @return the data array of a matrix 
    */
    inline const T* data() const;
    
    /** Copying the data from one matrix to another with different data type 
     * @param src the matrix we want to copy from
     * @param dst the matrix we want to copy to
    */
    template <class U>
    static void copyMat(const gmatrix<T>& src, gmatrix<U>& dst);

    /** Copying the data from one matrix to another with the same data type 
     * @param src the matrix we want to copy from
     * @param dst the matrix we want to copy to
    */
    static void copyMat(const gmatrix<T>& src, gmatrix<T>& dst);

    /** Copying the data from a vector to a matrix
     * @param src the vector we want to copy from
     * @param dst the matrix we want to copy to
    */
    static void copyVec(std::vector<T>& vec, const gmatrix<T>& mat);
    
    /**Allocate space for n number of element in the matrix where 
     * n is the size of the matrix that can be found 
    */
    void allocate();

    /**Delete the table holding the data and free the resources 
    */
    void release();

    /**
     * Set the element in the (nRow, nCol) to the value given.
     * @param nRow number of the row
     * @param nCol numnber of the column 
     * @param value Value
     * @return true if we can set element false in case of error 
     */
    bool set(int nRow, int nCol, T value);

    /**
     * Sets each element to the value val.
     * @param val Value
     */
    inline void fill(T e);

    /**
     * Returns the row at position r
     * @param r
     * @return row-vector of size 1xr
     */
    std::vector<T> row(int r) const;

    /**
     * Returns the row at position c
     * @param c
     * @return column-vector of size cx1
     */
    std::vector<T> column(int c) const;

    /**
     * Remove the r th row 
     * @param r Index of row to remove.
     */
    void delRow(int r);

    /**
     * Remove the c th row
     * @param c Index of column to remove.
    */
    void delCol(int c);

    /**
     * Print the n first number of rows and by default if noting is passed 5
     * will be the default value 
    */
    void head(int r = 5);

    /**
     * Print the n last number of rows and by default if noting is passed 5
     * will be the default value 
    */
    void tail(int r = 5);

    /**
     * Print element in the cell (i, j)
     * @param row the row index
     * @param col the column index 
    */
    void view(int row, int col);

    /**
     * Print a sub matrix starting from the cell (start_row, start_col)
     * to the cell (end_row, end_col)
     * @param start_row the starting position on the row index
     * @param end_row the end position on the row index
     * @param start_col the starting position on the column index
     * @param end_col the end position on the column index
    */
    void view(int start_row, int end_row, int start_col, int end_col);


    /**
     * Slice the matric to get a sub matrix starting from the cell 
     * (start_row, start_col) to the cell (end_row, end_col)
     * @param start_row the starting position on the row index
     * @param end_row the end position on the row index
     * @param start_col the starting position on the column index
     * @param end_col the end position on the column index
     * @note the matrix is not changed
    */
    gmatrix<T>& slice(int start_row, int end_row, int start_col, int end_col);

    /**
     * Transpose the matrix
     * @return the transposed matrix
     * @note the matrix is not changed
    */
    gmatrix<T>& transpose();

    /**
     * Printing operator overrloading for matricies 
     * @param os outout stream
     * @param matrix the matrix 
     * @return a reference to the output stream
    */
    friend std::ostream &operator<<(std::ostream &os, const gmatrix<T> &matrix);

    /**
     * Asignment operator: overwrite content of this matrix.
     * @param mat
     * @return
     */
    gmatrix<T>& operator=(const gmatrix<T>& mat);
    /**
     * Asignment operator: this is used with std::move.
     * @param mat
     * @return
     */
    gmatrix<T>& operator=(gmatrix<T>&& mat);

    /**
     * Check operator: true if the two matrices have the same value in the same index.
     * @param rhs the matrix we compare to it the "this" matrix 
     * @return true if equale false otherwise
     */
    bool operator==(const gmatrix<T>& rhs) const;

    /**
     * Plus operator
     * @param rhs right matrix
     * @param lhs left matrix
     * @return a gmatrix with the result of the plus operation
    */
    template <class U>
    friend gmatrix<U> operator+(const gmatrix<U>& lhs, const gmatrix<U>& rhs);
    /**
     * Plus operator
     * @param rhs right scalar
     * @param lhs left matrix
     * @return a gmatrix with the result of the plus operation
    */
    template <class U>
    friend gmatrix<U> operator+(const U& lhs, const gmatrix<U>& rhs);
    /**
     * Plus operator
     * @param rhs right matrix
     * @param lhs left scalar
     * @return a gmatrix with the result of the plus operation
    */
    template <class U>
    friend gmatrix<U> operator+(const gmatrix<U>& lhs, const U& rhs);

    /**
     * Substraction operator
     * @param rhs right matrix
     * @param lhs left matrix
     * @return a gmatrix with the result of the Substraction operation
    */
    template <class U>
    friend gmatrix<U> operator-(const gmatrix<U>& lhs, const gmatrix<U>& rhs);
    /**
     * Substraction operator
     * @param rhs right scalar
     * @param lhs left matrix
     * @return a gmatrix with the result of the Substraction operation
    */
    template <class U>
    friend gmatrix<U> operator-(const U& lhs, const gmatrix<U>& rhs);
    /**
     * Substraction operator
     * @param rhs right matrix
     * @param lhs left scalar
     * @return a gmatrix with the result of the Substraction operation
    */
    template <class U>
    friend gmatrix<U> operator-(const gmatrix<U>& lhs, const U& rhs);

    /**
     * multiplication operator
     * @param rhs right matrix
     * @param lhs left matrix
     * @return a gmatrix with the result of the multiplication operation
    */
    template <class U>
    friend gmatrix<U> operator*(const gmatrix<U>& lhs, const gmatrix<U>& rhs);
    /**
     * multiplication operator
     * @param rhs right scalar
     * @param lhs left matrix
     * @return a gmatrix with the result of the multiplication operation
    */
    template <class U>
    friend gmatrix<U> operator*(const U& lhs, const gmatrix<U>& rhs);
    /**
     * multiplication operator
     * @param rhs right matrix
     * @param lhs left scalar
     * @return a gmatrix with the result of the multiplication operation
    */
    template <class U>
    friend gmatrix<U> operator*(const gmatrix<U>& lhs, const U& rhs);



    
private:
    T* gmatrix_data;
    int rows_, cols_, size_;
};

// Constructors
template <class T>
gmatrix<T>::gmatrix() 
: rows_(1), cols_(1), size_(1) {
    gmatrix_data = allocate();
    gmatrix_data[0] = 0;
}

template <class T>
gmatrix<T>::gmatrix(int nRow, int nCol) 
: rows_(nRow), cols_(nCol), size_(nRow * nCol)  {
    allocate();
    gmatrix_data[0] = 0;
}


template <class T>
gmatrix<T>::gmatrix(int nRow, int nCol, const std::vector<T>& data) 
: rows_(nRow), cols_(nCol),  size_(nRow * nCol)  {
    allocate();
    copyVec(data, gmatrix_data);
}

template <class T>
template <class U>
gmatrix<T>::gmatrix(const gmatrix<U>& matrix) 
: gmatrix<T>(matrix.rows_, matrix.cols_m, matrix.size_) {
    allocate();
    copyMat(matrix, *this);
}

template <class T>
gmatrix<T>::gmatrix(const gmatrix<T>& matrix) 
: gmatrix<T>(matrix.rows_, matrix.cols_m, matrix.size_) {
    allocate();
    copyMat(matrix, *this);
}
template <class T>
gmatrix<T>::gmatrix(gmatrix<T>&& mat){
    this->gmatrix_data = std::move(mat.gmatrix_data);
    this->rows_ = mat.rows();
    this->cols_ = mat.cols();
    this->size_ = rows() * cols();
}

template <class T>
template <class U>
void gmatrix<T>::copyMat(const gmatrix<T>& src, gmatrix<U>& dst){
    std::copy(src, src + src.size(), dst);
}

template <class T>
void gmatrix<T>::copyMat(const gmatrix<T>& src, gmatrix<T>& dst){
    std::memcpy(dst, src, dst.size() * sizeof(T));
}

template <class T>
void gmatrix<T>::copyVec(std::vector<T>& vec, const gmatrix<T>& mat){
    if (vec.size() != mat.rows()){
        std::cout << "Error in copyVec" << std::end;
        // Needs work on error handling
        return;
    }    
    std::copy(vec.begin(), vec.end(), mat);
}

// Destructor
template <class T>
gmatrix<T>::~gmatrix() {
    release();
}

// Configuration functions
template <class T>
void gmatrix<T>::allocate() {
    gmatrix_data = new T[size_]();
}

template <class T>
void gmatrix<T>::release() {
    delete[] gmatrix_data;
    gmatrix_data = nullptr;
    size_ = rows_ = cols_ = 0;
}

template <class T>
bool gmatrix<T>::resize(int nRow, int nCol) {
    release();
    rows_ = nRow;
    cols_ = nCol;
    size_ = nRow * nCol;
    gmatrix_data = allocate();
    return gmatrix_data != nullptr;
}

template <class T>
int gmatrix<T>::linear_index(int nRow, int nCol) const {
    if (nRow < 0 || nRow >= rows()) {
        throw std::out_of_range("Invalid row index: " + std::to_string(nRow)+ " total rows are " + std::to_string(rows()));
    }
    if (nCol < 0 || nCol >= cols()) {
        throw std::out_of_range("Invalid column index: " + std::to_string(nCol));
    }
    return nRow * cols() + nCol;
}

template <class T>
inline bool gmatrix<T>::equalDimension(const gmatrix<T>& lhs, const gmatrix<T>& rhs) {
    return lhs.rows_ == rhs.rows_ && lhs.cols_ == rhs.cols_;
}

// Access Data
template <class T>
T gmatrix<T>::get(int nRow, int nCol) const {
    return gmatrix_data[linear_index(nRow, nCol)];
}

template <class T>
bool gmatrix<T>::set(int nRow, int nCol, T value) {
    int lindex = linear_index(nRow, nCol);
    gmatrix_data[lindex] = value;
    return true;
}

template<class T>
inline void gmatrix<T>::fill(T e){
    for (int i = 0; i < size(); ++i) {
        gmatrix_data[i] = e;
    }
}

template <class T>
std::vector<T> gmatrix<T>::row(int r) const{
    if (r < 0 || r >= rows())
        throw std::out_of_range("Row index out of range");
    
    std::vector<T> vec(gmatrix_data + r * cols(), gmatrix_data + (r + 1) * cols());
    return vec; 
}

template <class T>
std::vector<T> gmatrix<T>::column(int c) const{
    if (c < 0 || c >= cols())
        throw std::out_of_range("Column index out of range");
    std::vector<T> vec;
    for (int i = c; i < size(); i+=cols){
        /*
            we add cols to get to the same element in the next coloumn
            it took me a while to figure it out :>
        */
        vec.push_back(gmatrix_data[i]);
    }
    return vec;
}
template <class T>
void gmatrix<T>::delRow(int r) {
    if (r < 0 || r >= rows())
        throw std::out_of_range("Row index out of range");
    T* new_data = new T[(rows() - 1) * cols()]();
    int k = 0;
    for (int i = 0; i < rows(); ++i) {
        if (i == r) continue;
        for (int j = 0; j < cols(); ++j) 
            new_data[k++] = gmatrix_data[linear_index(i, j)];
    }
    // release(); this will change the row_ and col_ is better to jest desallocate here with delete[]
    delete[] gmatrix_data;
    --rows_;
    size_ = rows() * cols();
    gmatrix_data = new_data;
}

template <class T>
void gmatrix<T>::delCol(int c) {
    if (c < 0 || c >= cols())
        throw std::out_of_range("Column index out of range");
    T* new_data = new T[rows() * (cols() - 1)];
    int k = 0;
    for (int i = 0; i < rows(); ++i) {
        for (int j = 0; j < cols(); ++j) {
            if (j == c) continue;
            new_data[k++] = gmatrix_data[linear_index(i, j)];
        }
    }
    // release(); The same problem here 
    delete[] gmatrix_data;
    --cols_;
    size_ = rows() * cols();
    gmatrix_data = new_data;
}


template <class T>
void gmatrix<T>::head(int r) {
    int row = rows_() < r ? rows_() : r;
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < cols_; j++)
            std::cout << this.get(i, j) << "\t";
        std::cout << std::endl;
    }
}

template <class T>
void gmatrix<T>::tail(int r) {
    int row = rows_() < r ? rows_() : r;
    for (int i = rows_() - row; i < rows_(); i++) {
        for (int j = 0; j < cols_; j++)
            std::cout << this->get(i, j) << "\t";
        std::cout << std::endl;
    }
}

template <class T>
void gmatrix<T>::view(int row, int col){
    std::cout << this->get(row, col) << std::endl;
}

template <class T>
void gmatrix<T>::view(int start_row, int end_row, int start_col, int end_col) {
    if (start_row < 0 || 
        start_col < 0 || 
        end_row > this->rows() || 
        end_col > this->cols() || 
        start_row >= end_row || 
        start_col >= end_col) {
        throw std::out_of_range("The slicing parameters are out of bounds of the matrix size.");
    }

    for (int i = start_row; i < end_row; i++) {
        for (int j = start_col; j < end_col; j++) {
            std::cout << this->get(i, j) << "\t";
        }
        std::cout << std::endl;
    }
}

template <class T>
gmatrix<T>& gmatrix<T>::slice(int start_row, int end_row, int start_col, int end_col){
    if (start_row < 0 || 
        start_col < 0 || 
        end_row > this->rows() || 
        end_col > this->cols() || 
        start_row >= end_row || 
        start_col >= end_col) {
            throw std::out_of_range("The slicing parameters are out of bounds of the matrix size.");
    }
    
    gmatrix<T> mat(end_row - start_row, end_col - start_col);

    for (int i = start_row; i < end_row; i++) {
        for (int j = start_col; j < end_col; j++) {
            mat.set(i - start_row, j - start_col, this->get(i, j));
        }
    }
    return mat;
}

template <class T>
gmatrix<T>& gmatrix<T>::transpose(){
    gmatrix<T> mat(this->cols(), this->rows());

    for (int i = 0; i < this->cols(); i++) {
        for (int j = 0; j < this->rows(); j++) {
            mat.set(i, j, this->get(j, i));
        }
    }
    return mat;
}




template <class T>
int gmatrix<T>::rows() const {return rows_;}

template <class T>
int gmatrix<T>::cols() const {return cols_;}

template <class T>
int gmatrix<T>::size() const {return size_;}

template <class T>
inline T* gmatrix<T>::data(){ return gmatrix_data; }
template <class T>
inline const T* gmatrix<T>::data() const{ return gmatrix_data;}



template <typename T>
std::ostream& operator<<(std::ostream& os, const gmatrix<T>& matrix) {
    for (int i = 0; i < matrix.rows_; ++i) {
        for (int j = 0; j < matrix.cols_; ++j) {
            os << matrix.get(i, j) << ' ';
        }
        os << '\n';
    }
    return os;
}

template <class T>
gmatrix<T>& gmatrix<T>::operator=(const gmatrix<T>& mat){
    if (this != &mat) { // we don't have mat = mat self assign
        if (size_ != mat.size())
            gmatrix_data.resize(mat.rows(), mat.cols());

        rows_ = mat.rows();
        cols_ = mat.cols();
        size_ = mat.size();

        copyMat(mat, *this);
    }
    return *this;
}

template <class T>
gmatrix<T>& gmatrix<T>::operator=(gmatrix<T>&& mat){
    if (this != &mat) { // we don't have mat = mat self assign 
        this->gmatrix_data = std::move(mat.gmatrix_data);
        rows_ = mat.rows();
        cols_ = mat.cols();
        size_ = mat.size();
    }
    return *this;
}



// Operator overloader ==
template <class T>
bool gmatrix<T>::operator==(const gmatrix<T>& rhs) const {
    if (!equalDimension(*this, rhs))
        return false;
    for (int i = 0; i < size_; i++) {
        if (gmatrix_data[i] != rhs.gmatrix_data[i])
            return false;
    }
    return true;
}

// + Operator
template <class U>
gmatrix<U> operator+(const gmatrix<U>& lhs, const gmatrix<U>& rhs) {
    if (!gmatrix<U>::equalDimension(lhs, rhs))
        throw std::out_of_range("Dimensions don't match");
    gmatrix<U> result(lhs.rows_, lhs.cols_);
    for (int i = 0; i < lhs.size_; i++) {
        result.gmatrix_data[i] = lhs.gmatrix_data[i] + rhs.gmatrix_data[i];
    }
    return result;
}

// - Operator
template <class U>
gmatrix<U> operator-(const gmatrix<U>& lhs, const gmatrix<U>& rhs) {
    if (!gmatrix<U>::equalDimension(lhs, rhs))
        throw std::out_of_range("Dimensions don't match");
    gmatrix<U> result(lhs.rows_, lhs.cols_);
    for (int i = 0; i < lhs.size_; i++) {
        result.gmatrix_data[i] = lhs.gmatrix_data[i] - rhs.gmatrix_data[i];
    }
    return result;
}

// * Operator (Matrix Multiplication)
template <class U>
gmatrix<U> operator*(const gmatrix<U>& lhs, const gmatrix<U>& rhs) {
    if (lhs.cols_ != rhs.rows_)
        throw std::out_of_range("Dimensions do not match for matrix multiplication");
    gmatrix<U> result(lhs.rows_, rhs.cols_);
    for (int i = 0; i < lhs.rows_; ++i) {
        for (int j = 0; j < rhs.cols_; ++j) {
            U sum = 0;
            for (int k = 0; k < lhs.cols_; ++k) {
                sum += lhs.get(i, k) * rhs.get(k, j);
            }
            result.set(i, j, sum);
        }
    }
    return result;
}


#endif // GMATRIX_H_
