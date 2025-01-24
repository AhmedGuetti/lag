#include "gmatrix.h"
#include "gcomplex.h"
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

    std::cout << "All tests passed for gmatrix!" << std::endl;



    // // Test default constructor
    // gcomplex<double> c1;
    // assert(c1.real() == 0);
    // assert(c1.img() == 0);

    // // Test parameterized constructor
    // gcomplex<double> c2(3.0, 4.0);
    // assert(c2.real() == 3.0);
    // assert(c2.img() == 4.0);

    // // Test magnitude
    // assert(fabs(c2.magnitude() - 5.0) < 1e-9);

    // // Test argument
    // double expected_argument = atan2(4.0, 3.0);
    // assert(fabs(c2.argument() - expected_argument) < 1e-9);

    // // Test conjugate
    // gcomplex<double> c3 = c2.conjugate();
    // assert(c3.real() == 3.0);
    // assert(c3.img() == -4.0);

    // // Test addition
    // gcomplex<double> c4 = c2 + gcomplex<double>(1.0, 1.0);
    // assert(c4.real() == 4.0);
    // assert(c4.img() == 5.0);

    // // Test subtraction
    // gcomplex<double> c5 = c2 - gcomplex<double>(1.0, 1.0);
    // assert(c5.real() == 2.0);
    // assert(c5.img() == 3.0);

    // // Test multiplication
    // gcomplex<double> c6 = c2 * gcomplex<double>(2.0, 0.0);
    // assert(c6.real() == 6.0);
    // assert(c6.img() == 8.0);

    // // Test division
    // gcomplex<double> c7 = c2 / gcomplex<double>(1.0, 1.0);
    // assert(fabs(c7.real() - 3.5) < 1e-9);
    // assert(fabs(c7.img() - 0.5) < 1e-9);

    // // Test power
    // gcomplex<double> c8 = c2.power(2);
    // assert(fabs(c8.real() - (-7.0)) < 1e-9);
    // assert(fabs(c8.img() - (24.0)) < 1e-9);

    // // Test assignment operator
    // gcomplex<double> c9;
    // c9 = c2;
    // assert(c9.real() == 3.0);
    // assert(c9.img() == 4.0);

    // std::cout << "All tests passed for ggcomplex!" << std::endl;



    // Test 1: Constructor and Output
    gcomplex<double> c1(3.0, 4.0); // 3 + 4i
    gcomplex<double> c2(1.0, 2.0); // 1 + 2i
    std::cout << "c1: " << c1 << std::endl;
    std::cout << "c2: " << c2 << std::endl;

    // Test 2: Addition
    gcomplex<double> c3 = c1 + c2;
    std::cout << "c1 + c2 = " << c3 << std::endl;

    // Test 3: Subtraction
    gcomplex<double> c4 = c1 - c2;
    std::cout << "c1 - c2 = " << c4 << std::endl;

    // Test 4: Multiplication
    gcomplex<double> c5 = c1 * c2;
    std::cout << "c1 * c2 = " << c5 << std::endl;

    // Test 5: Division
    try {
        gcomplex<double> c6 = c1 / c2;
        std::cout << "c1 / c2 = " << c6 << std::endl;
    } catch (const std::exception& e) {
        std::cout << "Error during division: " << e.what() << std::endl;
    }

    // Test 6: Compound Addition
    c1 += c2;
    std::cout << "After c1 += c2, c1 = " << c1 << std::endl;

    // Test 7: Compound Subtraction
    c1 -= c2;
    std::cout << "After c1 -= c2, c1 = " << c1 << std::endl;

    // Test 8: Compound Multiplication
    c1 *= c2;
    std::cout << "After c1 *= c2, c1 = " << c1 << std::endl;

    // Test 9: Compound Division
    try {
        c1 /= c2;
        std::cout << "After c1 /= c2, c1 = " << c1 << std::endl;
    } catch (const std::exception& e) {
        std::cout << "Error during compound division: " << e.what() << std::endl;
    }

    // Test 10: Default Constructor and Copy Constructor
    gcomplex<double> c7; // Default constructor
    gcomplex<double> c8(c3); // Copy constructor
    std::cout << "Default constructor c7 = " << c7 << std::endl;
    std::cout << "Copy constructor c8 (copy of c3) = " << c8 << std::endl;

    // Test 11: Division by Zero
    try {
        gcomplex<double> c9(0.0, 0.0);
        gcomplex<double> result = c1 / c9;
        std::cout << "c1 / c9 = " << result << std::endl;
    } catch (const std::exception& e) {
        std::cout << "Error during division by zero: " << e.what() << std::endl;
    }




        // Create complex numbers
    gcomplex<double> z1(1, 2);  // 1 + 2i
    gcomplex<double> z2(3, 4);  // 3 + 4i

    std::cout << "Test Complex Numbers:" << std::endl;
    std::cout << "z1 = " << z1 << std::endl;
    std::cout << "z2 = " << z2 << std::endl;

    // Test addition, subtraction, multiplication, and division
    gcomplex<double> add_result = z1 + z2;
    gcomplex<double> sub_result = z1 - z2;
    gcomplex<double> mul_result = z1 * z2;
    gcomplex<double> div_result = z1 / z2;

    std::cout << "\nArithmetic operations:" << std::endl;
    std::cout << "z1 + z2 = " << add_result << std::endl;
    std::cout << "z1 - z2 = " << sub_result << std::endl;
    std::cout << "z1 * z2 = " << mul_result << std::endl;
    std::cout << "z1 / z2 = " << div_result << std::endl;

    // Test magnitude and argument (angle)
    std::cout << "\nMagnitude and argument:" << std::endl;
    std::cout << "Magnitude of z1 = " << z1.magnitude() << std::endl;
    std::cout << "Argument of z1 = " << z1.argument() << " radians" << std::endl;

    // Test trigonometric and exponential functions
    std::cout << "\nTrigonometric and exponential functions:" << std::endl;
    std::cout << "sin(z1) = " << z1.sin() << std::endl;
    std::cout << "cos(z1) = " << z1.cos() << std::endl;
    std::cout << "exp(z1) = " << z1.exp() << std::endl;

    // Test powers and logarithms
    std::cout << "\nPowers and logarithms:" << std::endl;
    std::cout << "z1^2 = " << z1.power(2) << std::endl;
    std::cout << "log(z1) = " << z1.log() << std::endl;

    // Test the conjugate function
    std::cout << "\nConjugate of z1 = " << z1.conjugate() << std::endl;

    // Test the square root
    std::cout << "\nSquare root of z1 = " << z1.sqrt() << std::endl;

    // Check assignment operator
    gcomplex<double> z3 = z1;
    std::cout << "\nAfter assigning z1 to z3: z3 = " << z3 << std::endl;

    // Check equality and inequality
    std::cout << "\nEquality and Inequality:" << std::endl;
    std::cout << "z1 == z3: " << (z1 == z3 ? "True" : "False") << std::endl;
    std::cout << "z1 != z2: " << (z1 != z2 ? "True" : "False") << std::endl;
    return 0;
}
