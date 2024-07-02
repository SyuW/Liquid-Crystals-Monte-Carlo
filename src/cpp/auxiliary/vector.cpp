#include <cmath>
#include <iostream>
#include <cassert>
#include "vector.hpp"

Vector::Vector(const Vector& otherVector)
{
    m_size = otherVector.GetSize();
    m_data = new double [m_size];
    for (int i=0; i<m_size; ++i)
    {
        m_data[i] = otherVector.m_data[i];
    }
}

Vector::Vector(int size)
{
    assert(size > 0);
    m_size = size;
    m_data = new double [m_size];
    for (int i=0; i<m_size; ++i)
    {
        m_data[i] = 0.0;
    }
}

Vector::~Vector()
{
    delete[] m_data;
}

int Vector::GetSize() const
{
    return m_size;
}

double& Vector::operator[] (int i)
{
    assert(i > -1);
    assert(i < m_size);
    return m_data[i];
}

double Vector::read(int i) const
{
    assert(i > -1);
    assert(i < m_size);
    return m_data[i];
}

double& Vector::operator() (int i)
{
    assert(i > 0);
    assert(i < m_size+1);
    return m_data[i-1];
}

Vector Vector::operator+() const
{
    Vector v(m_size);
    for (int i=0; i<m_size; ++i)
    {
        v[i] = m_data[i];
    }
    return v;
}

Vector Vector::operator-() const
{
    Vector v(m_size);
    for (int i=0; i<m_size; ++i)
    {
        v[i] = -m_data[i];
    }
    return v;
}

Vector& Vector::operator=(const Vector& otherVector)
{
    assert(m_size == otherVector.m_size);
    for (int i=0; i<m_size; ++i)
    {
        m_data[i] = otherVector.m_data[i];
    }
    return *this;
}

Vector Vector::operator+(const Vector& v1) const
{
    assert(m_size == v1.m_size);
    Vector v(m_size);
    for (int i=0; i<m_size; ++i)
    {
        v[i] = m_data[i] + v1.m_data[i];
    }
    return v;
}

Vector Vector::operator-(const Vector& v1) const
{
    assert(m_size == v1.m_size);
    Vector v(m_size);
    for (int i=0; i<m_size; ++i)
    {
        v[i] = m_data[i] - v1.m_data[i];
    }
    return v;
}

Vector Vector::operator*(double a) const
{
    Vector v(m_size);
    for (int i=0; i<m_size; ++i)
    {
        v[i] = a * m_data[i];
    }
    return v;
}

// dot product
double Vector::dot(const Vector& v1) const
{
    assert(m_size == v1.m_size);
    double dot_product;
    for (int i=0; i<m_size; ++i)
    {
        dot_product += m_data[i] * v1.m_data[i];
    }
    return dot_product;
}

// p-norm
double Vector::CalculateNorm(int p) const
{
    double norm_val;
    double sum = 0.0;
    for (int i=0; i<m_size; ++i)
    {
        sum += pow(fabs(m_data[i]), p);
    }
    norm_val = pow(sum, 1.0/(static_cast<double>(p)));
    return norm_val;
}

int length(const Vector& v)
{
    return v.m_size;
}