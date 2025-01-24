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

        gcomplex<T> sqrt();

        gcomplex<T> log();

        gcomplex<T> exp();

        gcomplex<T> pow(double c);

        gcomplex<T> sin();
        
        gcomplex<T> cos();

        gcomplex<T> tan();

        gcomplex<T> sec();

        gcomplex<T> csc();

        gcomplex<T> cot();

        gcomplex<T> sinh();

        gcomplex<T> cosh();

        gcomplex<T> tanh();

        gcomplex<T> sech();

        gcomplex<T> csch();

        gcomplex<T> coth();

        gcomplex<T> asin();

        gcomplex<T> acos();

        gcomplex<T> atan();

        gcomplex<T> asinh();

        gcomplex<T> acosh();

        gcomplex<T> atanh();

        // Assignment operators
        gcomplex<T>& operator=(const T re);
        gcomplex<T>& operator=(const gcomplex<T>& z);

        // Arithmetic operators with real numbers
        gcomplex<T>& operator+=(const T n);
        gcomplex<T>& operator-=(const T n);
        gcomplex<T>& operator*=(const T n);
        gcomplex<T>& operator/=(const T n);

        // Arithmetic operators with complex numbers
        gcomplex<T>& operator+=(const gcomplex<T>& z);
        gcomplex<T>& operator-=(const gcomplex<T>& z);
        gcomplex<T>& operator*=(const gcomplex<T>& z);
        gcomplex<T>& operator/=(const gcomplex<T>& z);

        // Equality operators
        bool operator==(const gcomplex<T>& z) const;
        bool operator!=(const gcomplex<T>& z) const;

        // Friend functions for arithmetic operators
        friend gcomplex<T> operator+(const T& lhs, const gcomplex<T>& rhs) {
            return gcomplex<T>(lhs + rhs.re, rhs.im);
        }

        friend gcomplex<T> operator+(const gcomplex<T>& lhs, const T& rhs) {
            return gcomplex<T>(lhs.re + rhs, lhs.im);
        }

        friend gcomplex<T> operator+(const gcomplex<T>& lhs, const gcomplex<T>& rhs) {
            return gcomplex<T>(lhs.re + rhs.re, lhs.im + rhs.im);
        }

        friend gcomplex<T> operator-(const T& lhs, const gcomplex<T>& rhs) {
            return gcomplex<T>(lhs - rhs.re, -rhs.im);
        }

        friend gcomplex<T> operator-(const gcomplex<T>& lhs, const T& rhs) {
            return gcomplex<T>(lhs.re - rhs, lhs.im);
        }

        friend gcomplex<T> operator-(const gcomplex<T>& lhs, const gcomplex<T>& rhs) {
            return gcomplex<T>(lhs.re - rhs.re, lhs.im - rhs.im);
        }

        friend gcomplex<T> operator*(const T& lhs, const gcomplex<T>& rhs) {
            return gcomplex<T>(lhs * rhs.re, lhs * rhs.im);
        }

        friend gcomplex<T> operator*(const gcomplex<T>& lhs, const T& rhs) {
            return gcomplex<T>(lhs.re * rhs, lhs.im * rhs);
        }

        friend gcomplex<T> operator*(const gcomplex<T>& lhs, const gcomplex<T>& rhs) {
            T realPart = lhs.re * rhs.re - lhs.im * rhs.im;
            T imaginaryPart = lhs.re * rhs.im + lhs.im * rhs.re;
            return gcomplex<T>(realPart, imaginaryPart);
        }

        friend gcomplex<T> operator/(const T& lhs, const gcomplex<T>& rhs) {
            T denominator = rhs.re * rhs.re + rhs.im * rhs.im;
            if (denominator == 0) throw std::runtime_error("Division by zero");
            T realPart = (lhs * rhs.re) / denominator;
            T imaginaryPart = -(lhs * rhs.im) / denominator;
            return gcomplex<T>(realPart, imaginaryPart);
        }

        friend gcomplex<T> operator/(const gcomplex<T>& lhs, const T& rhs) {
            if (rhs == 0) throw std::runtime_error("Division by zero");
            return gcomplex<T>(lhs.re / rhs, lhs.im / rhs);
        }

        friend gcomplex<T> operator/(const gcomplex<T>& lhs, const gcomplex<T>& rhs) {
            T denominator = rhs.re * rhs.re + rhs.im * rhs.im;
            if (denominator == 0) throw std::runtime_error("Division by zero");
            T realPart = (lhs.re * rhs.re + lhs.im * rhs.im) / denominator;
            T imaginaryPart = (lhs.im * rhs.re - lhs.re * rhs.im) / denominator;
            return gcomplex<T>(realPart, imaginaryPart);
        }

        // Output stream overload
        friend std::ostream& operator<<(std::ostream& out, const gcomplex<T>& z) {
            if (z.im >= 0) {
                out << z.re << " + " << z.im << "i";
            } else {
                out << z.re << " - " << -z.im << "i";
            }
            return out;
        }
};

// Implementations

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
