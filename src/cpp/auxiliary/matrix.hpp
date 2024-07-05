#ifndef MATRIXHEADERDEF
#define MATRIXHEADERDEF

class Matrix
{
public:
    // constructors
    Matrix(const Matrix& otherMatrix);
    Matrix(int numRows, int numCols);
    ~Matrix();

    // accessors
    int getNumberOfRows() const;
    int getNumberOfColumns() const;

    // printing
    void print() const;

    // overloaded operators
    double& operator() (int i, int j); // zero-based indexing
    Matrix operator+() const; // unary + operator
    Matrix operator-() const; // unary - operator
    Matrix& operator=(const Matrix& otherMatrix); // overwrite with another matrix
    Matrix operator+(const Matrix& otherMatrix) const; // add with another vector
    Matrix operator-(const Matrix& otherMatrix) const; // subtract with another vector
    Matrix operator*(double a) const; // scalar multiplication

private:
    double* m_data {};
    int m_num_rows {};
    int m_num_cols {};
    int m_num_elts {};

};

#endif