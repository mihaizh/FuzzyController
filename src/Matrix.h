#ifndef MATRIX_H
#define MATRIX_H

class Matrix
{
private:
    float** m_data = nullptr;

    int m_rows;
    int m_cols;

public:
    Matrix(int rows, int cols);
    Matrix(const Matrix& other);
    Matrix& operator=(const Matrix& other);
    ~Matrix();

    Matrix operator*(const Matrix& other);
    Matrix operator+(const Matrix& other);
    Matrix operator*(float val);

    int rows() const
    {
        return m_rows;
    }

    int cols() const
    {
        return m_cols;
    }

    float& at(int row, int col)
    {
        return m_data[row][col];
    }

    const float& at(int row, int col) const
    {
        return m_data[row][col];
    }
};

#endif
