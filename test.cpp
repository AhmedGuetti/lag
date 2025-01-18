#include "gmatrix.h"
#include <iostream>
#include <vector>
#include <cassert>

int main() {
    gmatrix<int> mat(3, 3);
    mat.fill(1);
    assert(mat.get(1, 1) == 1);

    mat.set(1, 1, 42);
    assert(mat.get(1, 1) == 42);

    std::vector<int> row = mat.row(1);
    assert(row[1] == 42);



    // for (int i = 0; i < 3; i++){
    //     for (int j = 0; j < 3; j++)
    //         std::cout << mat.get(i,j) << "\t";
       
    //     std::cout<<std::endl;
    // }
    
    mat.delRow(1);

    
    // for (int i = 0; i < 2; i++){
    //     for (int j = 0; j < 3; j++)
    //         std::cout << mat.get(i,j) << "\t";
    //     std::cout<<std::endl;
    // }
    
    assert(mat.rows() == 2);

    mat.delCol(1);
    assert(mat.cols() == 2);

    gmatrix<int> mat2(2, 2);
    mat2.fill(5);
    gmatrix<int> sum = mat + mat2;
    assert(sum.get(0, 0) == 6);

    gmatrix<int> diff = mat2 - mat;
    assert(diff.get(0, 0) == 4);

    gmatrix<int> mat3(2, 3);
    gmatrix<int> mat4(3, 2);
    mat3.fill(1);
    mat4.fill(2);
    gmatrix<int> product = mat3 * mat4;
    assert(product.rows() == 2 && product.cols() == 2);

    std::cout << "All tests passed!" << std::endl;
    return 0;
}
