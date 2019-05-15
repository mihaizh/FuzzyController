#ifndef TF_RUNNER_H
#define TF_RUNNER_H

#include "Matrix.h"

class TfRunner
{
private:
    Matrix m_A;
    Matrix m_B;
    Matrix m_C;
    Matrix m_D;

    Matrix m_x;
public:
    TfRunner(
        const Matrix& A,
        const Matrix& B,
        const Matrix& C,
        const Matrix& D);

    TfRunner(
        const Matrix& A,
        const Matrix& B,
        const Matrix& C,
        const Matrix& D,
        const Matrix& x0);

    float step(float u);
};

#endif
