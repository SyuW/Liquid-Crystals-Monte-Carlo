#ifndef VECTORHEADERDEF
#define VECTORHEADERDEF

#include <initializer_list>

class Vector
{
private:
    double* m_data {};
    int m_size {};

public:
    // constructors
    Vector(const Vector& otherVector);
    Vector(int size);
    ~Vector();

    // accessors
    int GetSize() const;

    // overloaded operators
    double& operator[] (int i); // zero-based indexing
    double& operator() (int i); // one-based indexing
    Vector operator+() const; // unary + operator
    Vector operator-() const; // unary - operator
    Vector& operator=(const Vector& otherVector); // overwrite with another vector
    Vector operator+(const Vector& v1) const; // add with another vector
    Vector operator-(const Vector& v1) const; // subtract with another vector
    Vector operator*(double a) const; // scalar multiplication

    // linear algebra
    double CalculateNorm(int p=2) const; // calculate the norm
    double dot(const Vector& v1) const;

    // miscellaneous functions
    double read(int i) const; // read only version of zero-based indexing
    friend int length(const Vector& v); // declare length() as a friend function
};

// function prototype for length()
int length(const Vector& v);

#endif