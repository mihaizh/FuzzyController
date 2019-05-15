#include "Matrix.h"

#include <new>
#include <assert.h>

Matrix::Matrix(int rows, int cols)
    : m_rows(rows),
    m_cols(cols)
{
    m_data = new float*[rows];
    for (int i = 0; i < rows; ++i)
    {
        m_data[i] = new float[cols];
    }
}

Matrix::Matrix(const Matrix& other)
    : m_rows(other.m_rows),
    m_cols(other.m_cols)
{
    m_data = new float*[m_rows];
    for (int i = 0; i < m_rows; ++i)
    {
        m_data[i] = new float[m_cols];
        for (int j = 0; j < m_cols; ++j)
        {
            m_data[i][j] = other.m_data[i][j];
        }
    }
}

Matrix& Matrix::operator=(const Matrix& other)
{
    if (this != &other)
    {
        this->~Matrix();
        new (this) Matrix(other);
    }

    return *this;
}

Matrix::~Matrix()
{
    for (int i = 0; i < m_rows; ++i)
    {
        delete[] m_data[i];
    }

    delete[] m_data;
}

Matrix Matrix::operator*(const Matrix& rhs)
{
    assert(m_cols == rhs.m_rows);

    Matrix res(m_rows, rhs.m_cols);
    for (int i = 0; i < m_rows; ++i)
    {
        for (int j = 0; j < rhs.m_cols; ++j)
        {
            res.m_data[i][j] = 0.F;
            for (int k = 0; k < m_cols; ++k)
            {
                res.m_data[i][j] += m_data[i][k] * rhs.m_data[k][j];
            }
        }
    }

    return res;
}

Matrix Matrix::operator+(const Matrix& rhs)
{
    assert(m_cols == rhs.m_cols);
    assert(m_rows == rhs.m_rows);

    Matrix res(*this);
    for (int i = 0; i < m_rows; ++i)
    {
        for (int j = 0; j < m_cols; ++j)
        {
            res.m_data[i][j] += rhs.m_data[i][j];
        }
    }

    return res;
}

Matrix Matrix::operator*(float val)
{
    Matrix res(*this);
    for (int i = 0; i < m_rows; ++i)
    {
        for (int j = 0; j < m_cols; ++j)
        {
            res.m_data[i][j] *= val;
        }
    }

    return res;
}

