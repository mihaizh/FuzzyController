#include "TfRunner.h"

TfRunner::TfRunner(
    const Matrix& A,
    const Matrix& B,
    const Matrix& C,
    const Matrix& D)
    : m_A(A), m_B(B), m_C(C), m_D(D),
    m_x(m_B.rows(), m_B.cols())
{
    for (int i = 0; i < m_x.rows(); ++i)
    {
        for (int j = 0; j < m_x.cols(); ++j)
        {
            m_x.at(i, j) = 0.F;
        }
    }
}

TfRunner::TfRunner(
    const Matrix& A,
    const Matrix& B,
    const Matrix& C,
    const Matrix& D,
    const Matrix& x0)
    : m_A(A), m_B(B), m_C(C), m_D(D), m_x(x0)
{
}

float TfRunner::step(float u)
{
    const Matrix y = (m_C * m_x) + (m_D * u);
    m_x = (m_A * m_x) + (m_B * u);

    return y.at(0, 0);
}

