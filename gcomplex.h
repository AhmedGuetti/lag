#if !defined(GCOMLEX_H_)
#define GCOMLEX_H_

#include <iostream>
#include <math.h>

/*
* z = m * e^(i*phi)
* z = m*cos(phi)+i*m*sin(phi)
*
*
*
**/
template <class T>
class gcomplex {
    private:
        // Complex c(re, im) ====> re + i * im;
        T re; // The real size
        T im; // The imaginary part
    
    public:
        /*
        * The default constructor to make the number 
        */
        gcomplex();

        /**
         * Constructs a Complex number 
         * with the real and imaginary number set to 0.
         * @param re real part 
         * @param im imaginary part
        */
        gcomplex(const T re, const T im);
        /**
         * Copy constructor
         * with the real and imaginary set to the part of the z Complex number.
         * @param z Complex number  
        */
        template <class U>
        gcomplex(const gcomplex<U>& z);
        gcomplex(const gcomplex<T>& z);

        ~gcomplex();

        /**
         * Rewrite the complex number.
         * @return complex number from the mag and phase
        */
        gcomplex<T> get_r_i(const double m, const double phi);

        /**
         * Calculates the conjugate of the Complex number.
         * @return conjugate Complex number
         */
        gcomplex<T> conjugate() const;

        /**
         * Calculates the magnitude of the Complex number.
         * @return magnitude as a double
         */
        double magnitude() const;

        /**
         * Calculates the argument (angle) of the Complex number.
         * @return argument as a double in radians
         */
        double argument() const;

        /**
         * Raises the Complex number to the power of n.
         * @param n exponent
         * @return resulting Complex number
         */
        gcomplex<T> power(const int n) const;

        /**
         * Raises the Complex number to the power of n.
         * @param n exponent
         * @return resulting Complex number
         */
        T get_re() const ;

        /**
         * Raises the Complex number to the power of n.
         * @param n exponent
         * @return resulting Complex number
         */
        T get_im() const;

        /**
         * Use magniture and phase to get the complex number as z = a + bi.
         * @param m the magniture of the complex number
         * @param phi the phase of the complex number also called the argument
         * @return Complex number represenrtaion
         */
        gcomplex<T> toRI(const double m, const double phi) ;

        /**
         * Compute the square root of the complex number.
         * @return The square root as a complex number.
        */
        gcomplex<T> sqrt();

        /**
         * Compute the natural logarithm of the complex number.
         * @return The natural logarithm as a complex number.
         */
        gcomplex<T> log();

        /**
         * Compute the exponential of the complex number.
         * @return The exponential as a complex number.
         */
        gcomplex<T> exp();

        /**
         * Raise the complex number to the power of c.
         * @param c The exponent.
         * @return The result of the exponentiation as a complex number.
         */
        gcomplex<T> pow(double c);

        /**
         * Compute the sine of the complex number.
         * @return The sine as a complex number.
         */
        gcomplex<T> sin();

        /**
         * Compute the cosine of the complex number.
         * @return The cosine as a complex number.
         */
        gcomplex<T> cos();

        /**
         * Compute the tangent of the complex number.
         * @return The tangent as a complex number.
         */
        gcomplex<T> tan();

        /**
         * Compute the secant of the complex number.
         * @return The secant as a complex number.
         */
        gcomplex<T> sec();

        /**
         * Compute the cosecant of the complex number.
         * @return The cosecant as a complex number.
         */
        gcomplex<T> csc();

        /**
         * Compute the cotangent of the complex number.
         * @return The cotangent as a complex number.
         */
        gcomplex<T> cot();

        /**
         * Compute the hyperbolic sine of the complex number.
         * @return The hyperbolic sine as a complex number.
         */
        gcomplex<T> sinh();

        /**
         * Compute the hyperbolic cosine of the complex number.
         * @return The hyperbolic cosine as a complex number.
         */
        gcomplex<T> cosh();

        /**
         * Compute the hyperbolic tangent of the complex number.
         * @return The hyperbolic tangent as a complex number.
         */
        gcomplex<T> tanh();

        /**
         * Compute the hyperbolic secant of the complex number.
         * @return The hyperbolic secant as a complex number.
         */
        gcomplex<T> sech();

        /**
         * Compute the hyperbolic cosecant of the complex number.
         * @return The hyperbolic cosecant as a complex number.
         */
        gcomplex<T> csch();

        /**
         * Compute the hyperbolic cotangent of the complex number.
         * @return The hyperbolic cotangent as a complex number.
         */
        gcomplex<T> coth();

        /**
         * Compute the inverse sine (arcsine) of the complex number.
         * @return The inverse sine as a complex number.
         */
        gcomplex<T> asin();

        /**
         * Compute the inverse cosine (arccosine) of the complex number.
         * @return The inverse cosine as a complex number.
         */
        gcomplex<T> acos();

        /**
         * Compute the inverse tangent (arctangent) of the complex number.
         * @return The inverse tangent as a complex number.
         */
        gcomplex<T> atan();

        /**
         * Compute the inverse hyperbolic sine of the complex number.
         * @return The inverse hyperbolic sine as a complex number.
         */
        gcomplex<T> asinh();

        /**
         * Compute the inverse hyperbolic cosine of the complex number.
         * @return The inverse hyperbolic cosine as a complex number.
         */
        gcomplex<T> acosh();

        /**
         * Compute the inverse hyperbolic tangent of the complex number.
         * @return The inverse hyperbolic tangent as a complex number.
         */
        gcomplex<T> atanh();


        /**
         * Assign a real number to this complex number, setting the imaginary part to zero.
         * @param re The real number to assign.
         * @return Reference to this complex number after assignment.
        */
        gcomplex<T>& operator=(const T re);
        
        /**
         * Copy assignment operator.
         * @param z The complex number to copy.
         * @return Reference to this complex number after assignment.
        */
        gcomplex<T>& operator=(const gcomplex<T>& z);

        /**
         * Add a real number to this complex number.
         * @param n The real number to add.
         * @return Reference to this complex number after addition.
        */
        gcomplex<T>& operator+=(const T n);

        /**
         * Subtract a real number from this complex number.
         * @param n The real number to subtract.
         * @return Reference to this complex number after subtraction.
        */
        gcomplex<T>& operator-=(const T n);

        /**
         * Multiply this complex number by a real number.
         * @param n The real number to multiply by.
         * @return Reference to this complex number after multiplication.
        */
        gcomplex<T>& operator*=(const T n);

        /**
         * Divide this complex number by a real number.
         * @param n The real number to divide by.
         * @return Reference to this complex number after division.
        */
        gcomplex<T>& operator/=(const T n);

        // Arithmetic operators with complex numbers
        /**
         * Add another complex number to this complex number.
         * @param z The complex number to add.
         * @return Reference to this complex number after addition.
        */
        gcomplex<T>& operator+=(const gcomplex<T>& z);

        /**
         * Subtract another complex number from this complex number.
         * @param z The complex number to subtract.
         * @return Reference to this complex number after subtraction.
        */
        gcomplex<T>& operator-=(const gcomplex<T>& z);

        /**
         * Multiply this complex number by another complex number.
         * @param z The complex number to multiply by.
         * @return Reference to this complex number after multiplication.
        */
        gcomplex<T>& operator*=(const gcomplex<T>& z);

        /**
         * Divide this complex number by another complex number.
         * @param z The complex number to divide by.
         * @return Reference to this complex number after division.
        */
        gcomplex<T>& operator/=(const gcomplex<T>& z);

        /**
         * Check equality between this complex number and another.
         * @param z The complex number to compare with.
         * @return True if both real and imaginary parts are equal, false otherwise.
         */
        bool operator==(const gcomplex<T>& z) const;

        /**
         * Check inequality between this complex number and another.
         * @param z The complex number to compare with.
         * @return True if real or imaginary parts are not equal, false otherwise.
         */
        bool operator!=(const gcomplex<T>& z) const;


        /**
         * Add a real number to a complex number.
         * @param lhs The real number to add.
         * @param rhs The complex number to add.
         * @return The sum as a new complex number.
         */
        friend gcomplex<T> operator+(const T& lhs, const gcomplex<T>& rhs) {
            return gcomplex<T>(lhs + rhs.re, rhs.im);
        }

        /**
         * Add a complex number to a real number.
         * @param lhs The complex number to add.
         * @param rhs The real number to add.
         * @return The sum as a new complex number.
         */
        friend gcomplex<T> operator+(const gcomplex<T>& lhs, const T& rhs) {
            return gcomplex<T>(lhs.re + rhs, lhs.im);
        }

        /**
         * Add two complex numbers.
         * @param lhs The first complex number to add.
         * @param rhs The second complex number to add.
         * @return The sum as a new complex number.
         */
        friend gcomplex<T> operator+(const gcomplex<T>& lhs, const gcomplex<T>& rhs) {
            return gcomplex<T>(lhs.re + rhs.re, lhs.im + rhs.im);
        }

        /**
         * Subtract a complex number from a real number.
         * @param lhs The real number to subtract from.
         * @param rhs The complex number to subtract.
         * @return The result as a new complex number.
         */
        friend gcomplex<T> operator-(const T& lhs, const gcomplex<T>& rhs) {
            return gcomplex<T>(lhs - rhs.re, -rhs.im);
        }

        /**
         * Subtract a real number from a complex number.
         * @param lhs The complex number to subtract from.
         * @param rhs The real number to subtract.
         * @return The result as a new complex number.
         */
        friend gcomplex<T> operator-(const gcomplex<T>& lhs, const T& rhs) {
            return gcomplex<T>(lhs.re - rhs, lhs.im);
        }

        /**
         * Subtract two complex numbers.
         * @param lhs The complex number to subtract from.
         * @param rhs The complex number to subtract.
         * @return The result as a new complex number.
        */
        friend gcomplex<T> operator-(const gcomplex<T>& lhs, const gcomplex<T>& rhs) {
            return gcomplex<T>(lhs.re - rhs.re, lhs.im - rhs.im);
        }

        /**
         * Multiply a real number by a complex number.
         * @param lhs The real number to multiply.
         * @param rhs The complex number to multiply.
         * @return The product as a new complex number.
        */
        friend gcomplex<T> operator*(const T& lhs, const gcomplex<T>& rhs) {
            return gcomplex<T>(lhs * rhs.re, lhs * rhs.im);
        }

        /**
         * Multiply a complex number by a real number.
         * @param lhs The complex number to multiply.
         * @param rhs The real number to multiply.
         * @return The product as a new complex number.
        */
        friend gcomplex<T> operator*(const gcomplex<T>& lhs, const T& rhs) {
            return gcomplex<T>(lhs.re * rhs, lhs.im * rhs);
        }

        /**
         * Multiply two complex numbers.
         * @param lhs The first complex number to multiply.
         * @param rhs The second complex number to multiply.
         * @return The product as a new complex number.
        */
        friend gcomplex<T> operator*(const gcomplex<T>& lhs, const gcomplex<T>& rhs) {
            T realPart = lhs.re * rhs.re - lhs.im * rhs.im;
            T imaginaryPart = lhs.re * rhs.im + lhs.im * rhs.re;
            return gcomplex<T>(realPart, imaginaryPart);
        }

        /**
         * Divide a real number by a complex number.
         * @param lhs The real number (dividend).
         * @param rhs The complex number (divisor).
         * @return The quotient as a new complex number.
         * @throws std::runtime_error if division by zero occurs.
        */
        friend gcomplex<T> operator/(const T& lhs, const gcomplex<T>& rhs) {
            T denominator = rhs.re * rhs.re + rhs.im * rhs.im;
            if (denominator == 0) throw std::runtime_error("Division by zero");
            T realPart = (lhs * rhs.re) / denominator;
            T imaginaryPart = -(lhs * rhs.im) / denominator;
            return gcomplex<T>(realPart, imaginaryPart);
        }

        /**
         * Divide a complex number by a real number.
         * @param lhs The complex number (dividend).
         * @param rhs The real number (divisor).
         * @return The quotient as a new complex number.
         * @throws std::runtime_error if division by zero occurs.
        */
        friend gcomplex<T> operator/(const gcomplex<T>& lhs, const T& rhs) {
            if (rhs == 0) throw std::runtime_error("Division by zero");
            return gcomplex<T>(lhs.re / rhs, lhs.im / rhs);
        }

        /**
         * Divide one complex number by another.
         * @param lhs The complex number (dividend).
         * @param rhs The complex number (divisor).
         * @return The quotient as a new complex number.
         * @throws std::runtime_error if division by zero occurs.
        */
        friend gcomplex<T> operator/(const gcomplex<T>& lhs, const gcomplex<T>& rhs) {
            T denominator = rhs.re * rhs.re + rhs.im * rhs.im;
            if (denominator == 0) throw std::runtime_error("Division by zero");
            T realPart = (lhs.re * rhs.re + lhs.im * rhs.im) / denominator;
            T imaginaryPart = (lhs.im * rhs.re - lhs.re * rhs.im) / denominator;
            return gcomplex<T>(realPart, imaginaryPart);
        }

        /**
         * Output the complex number to a stream in the form "a + bi" or "a - bi".
         * @param out The output stream to write to.
         * @param z The complex number to output.
         * @return The output stream after writing the complex number.
        */
        friend std::ostream& operator<<(std::ostream& out, const gcomplex<T>& z) {
            if (z.im >= 0) {
                out << z.re << " + " << z.im << "i";
            } else {
                out << z.re << " - " << -z.im << "i";
            }
            return out;
        }
};

// Contractors
template <class T>
gcomplex<T>::gcomplex() : re(0), im(0) {}

template <class T>
gcomplex<T>::gcomplex(const T re, const T im) : re(re), im(im) {}

template <class T>
gcomplex<T>::gcomplex(const gcomplex<T>& z) : re(z.re), im(z.im) {}

template <class T>
template <class U>
gcomplex<T>::gcomplex(const gcomplex<U>& z) : re(static_cast<T>(z.re)), im(static_cast<T>(z.im)) {}

template <class T>
gcomplex<T>::~gcomplex() {}

template <class T>
gcomplex<T> gcomplex<T>::toRI(const double m, const double phi) {
    return gcomplex<T>(m * std::cos(phi), m * std::sin(phi));
}

template <class T>
gcomplex<T> gcomplex<T>::conjugate() const {
    return gcomplex<T>(re, -im);
}

template <class T>
double gcomplex<T>::magnitude() const {
    return std::sqrt(re * re + im * im);
}

template <class T>
double gcomplex<T>::argument() const {
    if (re == 0 && im == 0) throw std::runtime_error("Undefined argument for zero complex number");
    return std::atan2(im, re);
}

template <class T>
gcomplex<T> gcomplex<T>::power(const int n) const {
    double r = magnitude();
    double a = argument();
    double rn = std::pow(r, n);
    return gcomplex<T>(rn * std::cos(n * a), rn * std::sin(n * a));
}

template <class T>
T gcomplex<T>::get_re() const {
    return this->re;
}

template <class T>
T gcomplex<T>::get_im() const {
    return this->im;
}




template <class T>
gcomplex<T> gcomplex<T>::sqrt()
{
    // no need to calculate those twice
    double inverse_sqrt2 = (1/std::sqrt(2));
    double a_b_2 =  std::sqrt(this->re * this->re + this->im * this->im);

    double realPart = inverse_sqrt2 * std::sqrt(a_b_2 + this->re);
    double imaginePart;

    if(this->im == 0){
        if(this->re <0) throw std::runtime_error("sqrt error");
        else return gcomplex<T>(realPart, 0);
    }
    imaginePart = this->im > 0 ? inverse_sqrt2 * std::sqrt(a_b_2 - this->re) : (-1) * inverse_sqrt2 * std::sqrt(a_b_2 - this->re);
    return gcomplex<T>(realPart, imaginePart);
  
}

template <class T>
gcomplex<T> gcomplex<T>::log()
{
  return *this;
}

template <class T>
gcomplex<T> gcomplex<T>::exp()
{
    return *this;
}

template <class T>
gcomplex<T> gcomplex<T>::pow(double c)
{
    return *this;
}

template <class T>
gcomplex<T> gcomplex<T>::sin()
{
return *this;
}

template <class T>
gcomplex<T> gcomplex<T>::cos()
{
return *this;
}


template <class T>
gcomplex<T> gcomplex<T>::tan() { return *this;}

template <class T>
gcomplex<T> gcomplex<T>::sec() {return *this;}

template <class T>
gcomplex<T> gcomplex<T>::csc(){return *this;}

template <class T>
gcomplex<T> gcomplex<T>::cot() {return *this;}

template <class T>
gcomplex<T> gcomplex<T>::sinh() {return *this;}

template <class T>
gcomplex<T> gcomplex<T>::cosh() {return *this;}

template <class T>
gcomplex<T> gcomplex<T>::tanh() {return *this;}

template <class T>
gcomplex<T> gcomplex<T>::sech() {return *this;}

template <class T>
gcomplex<T> gcomplex<T>::csch(){return *this;}

template <class T>
gcomplex<T> gcomplex<T>::coth() {return *this;}
template <class T>
gcomplex<T> gcomplex<T>::asin() {return *this;}

template <class T>
gcomplex<T> gcomplex<T>::acos() {return *this;}

template <class T>
gcomplex<T> gcomplex<T>::atan(){return *this;}

template <class T>
gcomplex<T> gcomplex<T>::asinh(){return *this;}

template <class T>
gcomplex<T> gcomplex<T>::acosh() {return *this;}

template <class T>
gcomplex<T> gcomplex<T>::atanh() {return *this;}


template <class T>
gcomplex<T>& gcomplex<T>::operator=(const T re) {
    this->re = re;
    this->im = 0;
    return *this;
}

template <class T>
gcomplex<T>& gcomplex<T>::operator=(const gcomplex<T>& z) {
    this->re = z.re;
    this->im = z.im;
    return *this;
}

template <class T>
gcomplex<T>& gcomplex<T>::operator+=(const T n) {
    re += n;
    return *this;
}

template <class T>
gcomplex<T>& gcomplex<T>::operator-=(const T n) {
    re -= n;
    return *this;
}

template <class T>
gcomplex<T>& gcomplex<T>::operator*=(const T n) {
    re *= n;
    im *= n;
    return *this;
}

template <class T>
gcomplex<T>& gcomplex<T>::operator/=(const T n) {
    if (n == 0) throw std::runtime_error("Division by zero");
    re /= n;
    im /= n;
    return *this;
}

template <class T>
gcomplex<T>& gcomplex<T>::operator+=(const gcomplex<T>& z) {
    re += z.re;
    im += z.im;
    return *this;
}

template <class T>
gcomplex<T>& gcomplex<T>::operator-=(const gcomplex<T>& z) {
    re -= z.re;
    im -= z.im;
    return *this;
}

template <class T>
gcomplex<T>& gcomplex<T>::operator*=(const gcomplex<T>& z) {
    T realPart = re * z.re - im * z.im;
    T imaginaryPart = re * z.im + im * z.re;
    re = realPart;
    im = imaginaryPart;
    return *this;
}

template <class T>
gcomplex<T>& gcomplex<T>::operator/=(const gcomplex<T>& z) {
    T denominator = z.re * z.re + z.im * z.im;
    if (denominator == 0) throw std::runtime_error("Division by zero");
    T realPart = (re * z.re + im * z.im) / denominator;
    T imaginaryPart = (im * z.re - re * z.im) / denominator;
    re = realPart;
    im = imaginaryPart;
    return *this;
}

template <class T>
bool gcomplex<T>::operator==(const gcomplex<T>& z) const {
    return re == z.re && im == z.im;
}

template <class T>
bool gcomplex<T>::operator!=(const gcomplex<T>& z) const {
    return !(*this == z);
}

#endif // GCOMPLEX_H_
