#ifndef GMATRIX_H_
#define GMATRIX_H_

#include <stdexcept>

template <class T>
class gmatrix {
public:
    // Constructors
    gmatrix();
    gmatrix(int nRow, int nCol);
    gmatrix(int nRow, int nCol, const T* data);
    gmatrix(const gmatrix<T>& matrix);

    // Destructor
    ~gmatrix();

    // Configuration
    bool resize(int nRow, int nCol);
    int linear_index(int nRow, int nCol) const;
    static bool equalDimension(const gmatrix<T>& lhs, const gmatrix<T>& rhs);

    // Access Data
    T get(int nRow, int nCol) const;
    bool set(int nRow, int nCol, T value);
    int rows() const;
    int cols() const;
    int size() const;

    // Operator overloader
    bool operator==(const gmatrix<T>& rhs) const;

    // Overload +, -, and *
    template <class U>
    friend gmatrix<U> operator+(const gmatrix<U>& lhs, const gmatrix<U>& rhs);
    template <class U>
    friend gmatrix<U> operator+(const U& lhs, const gmatrix<U>& rhs);
    template <class U>
    friend gmatrix<U> operator+(const gmatrix<U>& lhs, const U& rhs);

    template <class U>
    friend gmatrix<U> operator-(const gmatrix<U>& lhs, const gmatrix<U>& rhs);
    template <class U>
    friend gmatrix<U> operator-(const U& lhs, const gmatrix<U>& rhs);
    template <class U>
    friend gmatrix<U> operator-(const gmatrix<U>& lhs, const U& rhs);

    template <class U>
    friend gmatrix<U> operator*(const gmatrix<U>& lhs, const gmatrix<U>& rhs);
    template <class U>
    friend gmatrix<U> operator*(const U& lhs, const gmatrix<U>& rhs);
    template <class U>
    friend gmatrix<U> operator*(const gmatrix<U>& lhs, const U& rhs);

private:
    T* gmatrix_data;
    int rows_, cols_, size_;
};

// Constructors
template <class T>
gmatrix<T>::gmatrix() : rows_(1), cols_(1), size_(1) {
    gmatrix_data = new T[size_];
    gmatrix_data[0] = 0;
}

template <class T>
gmatrix<T>::gmatrix(int nRow, int nCol) : rows_(nRow), cols_(nCol), size_(nRow * nCol) {
    gmatrix_data = new T[size_]();
}

template <class T>
gmatrix<T>::gmatrix(int nRow, int nCol, const T* data) : rows_(nRow), cols_(nCol), size_(nRow * nCol) {
    gmatrix_data = new T[size_];
    for (int i = 0; i < size_; i++)
        gmatrix_data[i] = data[i];
}

template <class T>
gmatrix<T>::gmatrix(const gmatrix<T>& matrix) : rows_(matrix.rows_), cols_(matrix.cols_), size_(matrix.size_) {
    gmatrix_data = new T[size_];
    for (int i = 0; i < size_; i++)
        gmatrix_data[i] = matrix.gmatrix_data[i];
}

// Destructor
template <class T>
gmatrix<T>::~gmatrix() {
    delete[] gmatrix_data;
}

// Configuration functions
template <class T>
bool gmatrix<T>::resize(int nRow, int nCol) {
    delete[] gmatrix_data;
    rows_ = nRow;
    cols_ = nCol;
    size_ = nRow * nCol;
    gmatrix_data = new T[size_]();
    return gmatrix_data != nullptr;
}

template <class T>
int gmatrix<T>::linear_index(int nRow, int nCol) const {
    if (nRow >= 0 && nRow < rows_ && nCol >= 0 && nCol < cols_) {
        return nRow * cols_ + nCol;
    }
    throw std::out_of_range("Invalid row or column index");
}

template <class T>
bool gmatrix<T>::equalDimension(const gmatrix<T>& lhs, const gmatrix<T>& rhs) {
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

template <class T>
int gmatrix<T>::rows() const {
    return rows_;
}

template <class T>
int gmatrix<T>::cols() const {
    return cols_;
}

template <class T>
int gmatrix<T>::size() const {
    return size_;
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
