#include <cmath>
#include <iostream>
#include <cassert>
#include "matrix.hpp"

// copy constructor
Matrix::Matrix(const Matrix& otherMatrix)
{
    m_num_rows = otherMatrix.getNumberOfRows();
    m_num_cols = otherMatrix.getNumberOfColumns();
    m_num_elts = m_num_rows * m_num_cols;
    m_data = new (std::nothrow) double [m_num_elts];
    for (int i=0; i<m_num_elts; ++i)
    {
        m_data[i] = otherMatrix.m_data[i];
    }
}

// constructor (with zero'd out entries)
Matrix::Matrix(int numRows, int numCols)
{
    assert(numRows > 0 && numCols > 0);
    m_num_rows = numRows;
    m_num_cols = numCols;
    m_num_elts = m_num_rows * m_num_cols;
    m_data = new (std::nothrow) double [m_num_elts];

    if (!m_data)
    {
        std::cerr << "Could not allocate new memory\n";
    }

    for (int i=0; i<m_num_elts; ++i)
    {
        m_data[i] = 0.0;
    }
}

// destructor
Matrix::~Matrix()
{
    delete[] m_data;
}

// accessor for number of rows
int Matrix::getNumberOfRows() const
{
    return m_num_rows;
}

// accessor for number of columns
int Matrix::getNumberOfColumns() const
{
    return m_num_cols;
}

// accessor for number of elements
int Matrix::getNumberOfElements() const
{
    return m_num_elts;
}

// print method
void Matrix::print() const
{
    for (int i = 0; i < m_num_rows; ++i)
    {
        for (int j = 0; j < m_num_cols; ++j)
        {
            std::cout << m_data[m_num_cols*i + j] << " ";
        }
        std::cout << "\n";
    }
}

// 
void Matrix::setsize(int m, int n)
{
    
}

// zero-based indexing
double& Matrix::operator() (int i, int j)
{
    assert(i >= 0 && i < m_num_rows);
    assert(j >= 0 && j < m_num_cols);

    return m_data[m_num_cols*i + j];
}

// zero-based indexing with flattened index
double& Matrix::operator[] (int i)
{
    assert(i >= 0 && i < m_num_elts);

    return m_data[i];
}

// unary identity
Matrix Matrix::operator+() const
{
    Matrix M(m_num_rows, m_num_cols);
    for (int i=0; i<m_num_elts; ++i)
    {
        M.m_data[i] = m_data[i];
    }
    return M;
}

// unary negation
Matrix Matrix::operator-() const
{
    Matrix M(m_num_rows, m_num_cols);
    for (int i=0; i<m_num_elts; ++i)
    {
        M.m_data[i] = -m_data[i];
    }
    return M;
}

// copy assignment
Matrix& Matrix::operator=(const Matrix& otherMatrix)
{
    assert(m_num_rows == otherMatrix.m_num_rows);
    assert(m_num_cols == otherMatrix.m_num_cols);
    for (int i=0; i<m_num_elts; ++i)
    {
        m_data[i] = otherMatrix.m_data[i];
    }
    return *this;
}

// binary addition
Matrix Matrix::operator+(const Matrix& otherMatrix) const
{
    assert(m_num_rows == otherMatrix.m_num_rows);
    assert(m_num_cols == otherMatrix.m_num_cols);
    Matrix M(m_num_rows, m_num_cols);
    for (int i=0; i<m_num_elts; ++i)
    {
        M.m_data[i] = m_data[i] + otherMatrix.m_data[i];
    }
    return M;
}

// binary subtraction
Matrix Matrix::operator-(const Matrix& otherMatrix) const
{
    assert(m_num_rows == otherMatrix.m_num_rows);
    assert(m_num_cols == otherMatrix.m_num_cols);
    Matrix M(m_num_rows, m_num_cols);
    for (int i=0; i<m_num_elts; ++i)
    {
        M.m_data[i] = m_data[i] - otherMatrix.m_data[i];
    }
    return M;
}

// scalar multiplication
Matrix Matrix::operator*(double a) const
{
    Matrix M(m_num_rows, m_num_cols);
    for (int i=0; i<m_num_elts; ++i)
    {
        M.m_data[i] = a * m_data[i];
    }
    return M;
}
