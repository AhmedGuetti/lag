#ifndef GCOMPLEX_H_
#define GCOMPLEX_H_

#include <iostream>
#include <cmath>
#include <stdexcept>

template <class T>
class gcomplex {
private:
    T re; // Real part
    T im; // Imaginary part

public:
    // Default constructor
    gcomplex();

    // Constructor with real and imaginary parts
    gcomplex(const T re, const T im);

    // Copy constructors
    gcomplex(const gcomplex<T>& z);
    template <class U>
    gcomplex(const gcomplex<U>& z);

    // Destructor
    ~gcomplex();

    // Create a complex number from magnitude and phase
    static gcomplex<T> toRI(const double m, const double phi);

    // Conjugate of the complex number
    gcomplex<T> conjugate() const;

    // Magnitude of the complex number
    double magnitude() const;

    // Argument (angle) of the complex number in radians
    double argument() const;

    // Raise the complex number to the power of n
    gcomplex<T> power(const int n) const;

    // Getters for real and imaginary parts
    T real() const;
    T img() const;

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
T gcomplex<T>::real() const {
    return re;
}

template <class T>
T gcomplex<T>::img() const {
    return im;
}

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
