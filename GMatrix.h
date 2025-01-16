#if !defined(GMATRIX_H_)
#define GMATRIX_H_

#include <stdexcept>

template <class T>
class gmatrix
{
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
        int linear_index(int nRow, int nCol);
        bool equalDimension(const gmatrix<T> lhs, const gmatrix<T> rhs) const;

        // Access Data
        T get(int nRow, int nCol);
        bool set(int nRow, int nCol, T value);
        int rows() const;
        int cols() const;

        // Operator overloader
        bool operator == (const gmatrix<T>& rhs);

        // Overload +, - and *
        template <class U> 
        friend gmatrix<U> operator+ (const gmatrix<U> lhs, const gmatrix<U> rhs);
        template <class U> 
        friend gmatrix<U> operator+ (const U& lhs, const gmatrix<U> rhs);
        template <class U> 
        friend gmatrix<U> operator+ (const gmatrix<U> lhs, const U& rhs);

        template <class U> 
        friend gmatrix<U> operator- (const gmatrix<U> lhs, const gmatrix<U> rhs);
        template <class U> 
        friend gmatrix<U> operator- (const U& lhs, const gmatrix<U> rhs);
        template <class U> 
        friend gmatrix<U> operator- (const gmatrix<U> lhs, const U& rhs);
        
        template <class U> 
        friend gmatrix<U> operator* (const gmatrix<U> lhs, const gmatrix<U> rhs);
        template <class U> 
        friend gmatrix<U> operator* (const U& lhs, const gmatrix<U> rhs);
        template <class U> 
        friend gmatrix<U> operator* (const gmatrix<U> lhs, const U& rhs);

    private:
        T* gmatrix_data;
        int rows_, cols_, count;
};




// Constroctors
template<class T>
gmatrix<T>::gmatrix() 
: rows_(1) , cols_(1), count(1)
{
    gmatrix_data = new T[count];
    gmatrix_data[0] = 0.0;
}
template<class T>
gmatrix<T>::gmatrix(int nRow, int nCol)
: rows_(nRow), cols_(nCol), count(nRow * nCol)
{
    gmatrix_data = new T[count];
    for (int i = 0; i < count; i++)
        gmatrix_data[0] = 0.0;
}
template<class T>
gmatrix<T>::gmatrix(int nRow, int nCol, const T* data)
: rows_(nRow), cols_(nCol), count(nRow * nCol)
{
    gmatrix_data = new T[count];
    for (int i = 0; i < count; i++)
        gmatrix_data[0] = data[i];
}

template<class T>
gmatrix<T>::gmatrix(const gmatrix<T>& matrix)
: rows_(matrix.rows_), cols_(matrix.cols_), count(matrix.count)
{
    gmatrix_data = new T[count];
    for (int i = 0; i < count; i++)
        gmatrix_data[0] = matrix.gmatrix_data[i];
}


// Destructor 
template<class T>
gmatrix<T>::~gmatrix(){
    if(gmatrix_data != nullptr)
        delete[] gmatrix_data;
    
}

// Configuration functions
template <class T>
bool gmatrix<T>::resize(int nRow, int nCol){
    rows_ = nRow;
    cols_ = nCol;
    count = nRow * nCol;
    delete[] gmatrix_data;
    if(gmatrix_data != nullptr)
        for (int i = 0; i < count; i++){
            gmatrix_data[i] = 0.0;
            return true;
        }
    return false;
}
template <class T>
int linear_index (int nRow, int nCol)
{
    if ((rows_ < nRow) && (rows_ > 0) && (cols_ < nCol) && (cols_ > 0))
        return (rows_ * nCol) +  cols_;
    return -1;
}
template <class T>
bool gmatrix<T>::equalDimension(const gmatrix<T> lhs, const gmatrix<T> rhs) const
{
    return lhs.cols() == rhs.cols() && lhs.rows() == rhs.rows();
}

// Access Function
template <class T>
T gmatrix<T>::get(int nRow, int nCol){
    int lindex = linear_index(nRow, nCol);
    if(lindex >=0)
        return gmatrix_data[lindex];
    throw std::out_of_range("Invalid row or column index");
}


template <class T>
bool gmatrix<T>::set(int nRow, int nCol, T value){

    int lindex = linear_index(nRow, nCol);
    if(lindex >= 0) {
        gmatrix_data[lindex] = value;
        return true;
    }
    return false;
}
template <class T>
int gmatrix<T>::rows() const {return rows_;}

template <class T>
int gmatrix<T>::cols() const {return cols_;}


// + Operator
template <class U> 
gmatrix<U> operator+ (const gmatrix<U> lhs, const gmatrix<U> rhs){
    if(!equalDimension(lhs, rhs))
        throw std::out_of_range("Dimetion are  dosn't match");
    int numberCount = lhs.rows_ * lhs.cols_;
    T* temp = new T[numberCount];
    // c_ij = a_ij + b_ij
    for (int i = 0; i < count; i++)
        temp[i] = lhs.gmatrix_data[i] + rhs.gmatrix_data[i];
    gmatrix<T> result(lhs.rows_, lhs.cols_, temp);
    delete[] temp;
    return  result;
}


template <class U> 
gmatrix<U> operator+ (const U& lhs, const gmatrix<U> rhs){
    int numberCount = lhs.rows_ * lhs.cols_;
    T* temp = new T[numberCount];
    // (c_ij = scalar + b_ij) we do the sum element wise
    for (int i = 0; i < count; i++)
        temp[i] = lhs + rhs.gmatrix_data[i];
    gmatrix<T> result(lhs.rows_, lhs.cols_, temp);
    delete[] temp;
    return  result;
}


template <class U> 
gmatrix<U> operator+ (const gmatrix<U> lhs, const U& rhs){
    int numberCount = lhs.rows_ * lhs.cols_;
    T* temp = new T[numberCount];
    // (c_ij = a_ij + scalar) we do the sum element wise
    for (int i = 0; i < count; i++)
        temp[i] = lhs.gmatrix_data[i] + rhs;
    gmatrix<T> result(lhs.rows_, lhs.cols_, temp);
    delete[] temp;
    return  result;
}


// - Operator
template <class U> 
gmatrix<U> operator- (const gmatrix<U> lhs, const gmatrix<U> rhs){
    if(!equalDimension(lhs, rhs))
        throw std::out_of_range("Dimetion are  dosn't match");
    int numberCount = lhs.rows_ * lhs.cols_;
    T* temp = new T[numberCount];
    for (int i = 0; i < count; i++)
        temp[i] = lhs.gmatrix_data[i] - rhs.gmatrix_data[i];
    gmatrix<T> result(lhs.rows_, lhs.cols_, temp);
    delete[] temp;
    return  result;
}

template <class U> 
gmatrix<U> operator- (const U& lhs, const gmatrix<U> rhs){
    int numberCount = lhs.rows_ * lhs.cols_;
    T* temp = new T[numberCount];
    for (int i = 0; i < count; i++)
        temp[i] = lhs - rhs.gmatrix_data[i];
    gmatrix<T> result(lhs.rows_, lhs.cols_, temp);
    delete[] temp;
    return  result;
}

template <class U> 
gmatrix<U> operator- (const gmatrix<U> lhs, const U& rhs){
    int numberCount = lhs.rows_ * lhs.cols_;
    T* temp = new T[numberCount];
    for (int i = 0; i < count; i++)
        temp[i] = lhs.gmatrix_data[i] - rhs;
    gmatrix<T> result(lhs.rows_, lhs.cols_, temp);
    delete[] temp;
    return  result;
}



// * Operator (Matrix Multiplication)
template <class U> 
gmatrix<U> operator* (const gmatrix<U> lhs, const gmatrix<U> rhs){

    // have to be worked on 
}


template <class U> 
gmatrix<U> operator* (const U& lhs, const gmatrix<U> rhs){
    int numberCount = lhs.rows_ * lhs.cols_;
    T* temp = new T[numberCount];
    for (int i = 0; i < count; i++)
        temp[i] = lhs * rhs.gmatrix_data[i];
    gmatrix<T> result(lhs.rows_, lhs.cols_, temp);
    delete[] temp;
    return  result;
}

template <class U> 
gmatrix<U> operator* (const gmatrix<U> lhs, const U& rhs){
    int numberCount = lhs.rows_ * lhs.cols_;
    T* temp = new T[numberCount];
    for (int i = 0; i < count; i++)
        temp[i] = lhs.gmatrix_data[i] * rhs;
    gmatrix<T> result(lhs.rows_, lhs.cols_, temp);
    delete[] temp;
    return  result;
}





#endif // GMATRIX_H_
