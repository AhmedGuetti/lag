#include <iostream>
#include "./gmatrix.h"

int main() {
    try {
        // Create matrices
        int data1[] = {1, 2, 3, 4, 5, 6};
        int data2[] = {7, 8, 9, 10, 11, 12};

        // Matrix A (2x3)
        gmatrix<int> mat1(2, 3, data1);

        // Matrix B (3x2)
        gmatrix<int> mat2(3, 2, data2);

        // Perform Matrix Multiplication
        gmatrix<int> result = mat1 * mat2;

        // Print Matrix A
        std::cout << "Matrix A (2x3):" << std::endl;
        for (int i = 0; i < mat1.rows(); ++i) {
            for (int j = 0; j < mat1.cols(); ++j) {
                std::cout << mat1.get(i, j) << " ";
            }
            std::cout << std::endl;
        }

        // Print Matrix B
        std::cout << "Matrix B (3x2):" << std::endl;
        for (int i = 0; i < mat2.rows(); ++i) {
            for (int j = 0; j < mat2.cols(); ++j) {
                std::cout << mat2.get(i, j) << " ";
            }
            std::cout << std::endl;
        }

        // Print Result Matrix
        std::cout << "Result Matrix (2x2):" << std::endl;
        for (int i = 0; i < result.rows(); ++i) {
            for (int j = 0; j < result.cols(); ++j) {
                std::cout << result.get(i, j) << " ";
            }
            std::cout << std::endl;
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}
