#include "FuzzyController.h"

#include <math.h>
#include <assert.h>

Range::Range()
    : m_min(0.F),
    m_max(0.F),
    m_step(0.F),
    m_size(0),
    m_data(nullptr)
{
}

Range::Range(float min, float max, float step)
    : m_min(min),
    m_max(max),
    m_step(step),
    m_size(0)
{
    assert((min < max) && "min should be smaller than max");

    m_size = (int)floorf((fabsf(min) + fabs(max)) / step);
    m_data = new float[m_size];

    for (int i = 0; i < m_size; ++i)
    {
        m_data[i] = min + (i * step);
    }
}

Range::Range(const Range& other)
    : m_min(other.m_min),
    m_max(other.m_max),
    m_step(other.m_step),
    m_size(other.m_size)
{
    m_data = new float[m_size];
    for (int i = 0; i < m_size; ++i)
    {
        m_data[i] = other.m_data[i];
    }
}

Range& Range::operator=(const Range& other)
{
    if (this != &other)
    {
        m_min = other.m_min;
        m_max = other.m_max;
        m_step = other.m_step;
        m_size = other.m_size;

        delete[] m_data;
        m_data = new float[m_size];

        for (int i = 0; i < m_size; ++i)
        {
            m_data[i] = other.m_data[i];
        }
    }

    return *this;
}

int Range::pre(float val)
{
    if (val < m_min)
    {
        return 0;
    }

    if (val > m_max)
    {
        return (m_size - 1);
    }

    for (int i = 1; i < m_size; ++i)
    {
        if ((m_data[i - 1] <= val) && (m_data[i] >= val))
        {
            return (i - 1);
        }
    }

    assert(false && "Range for value could not be found");
    return 0;
}

float FuzzyController::centroid(int i, int j)
{
    // Assumption: hits at most 2 MFs
    float h[N][2] = { { 1.F, 1.F }, { 1.F, 1.F }, { 1.F, 1.F } };
    int mf[N][2];

    memset(mf, 0, sizeof(mf));

    // heights
    for (int l = 0; l < (N - 1); ++l)
    {
        int idx = 0;
        for (int k = 0; k < m_mfs[l].size(); ++k)
        {
            const float x = m_ranges[l].at((l == 0) ? i : j);
            const float ret = getH(m_mfs[l][k], x); // ??
            if (ret > std::numeric_limits<float>::epsilon())
            {
                h[l][idx] = ret;
                mf[l][idx] = k;
                ++idx;
            }
        }
    }

    // min
    float minh[N][2];
    for (int l = 0; l < (N - 1); ++l)
    {
        for (int m = 0; m < 2; ++m)
        {
            if (h[0][m] < h[1][l])
            {
                minh[l][m] = h[0][m];
            }
            else
            {
                minh[l][m] = h[1][l];
            }
        }
    }

    // centers
    float num = 0.F;
    float den = 0.F;
    for (int l = 0; l < (N - 1); ++l)
    {
        for (int m = 0; m < 2; ++m)
        {
            const float c = getC(m_mfs[2][m_rules[mf[0][m]][mf[1][l]]], minh[l][m]);
            num += c * minh[l][m];
            den += minh[l][m];
        }
    }

    return num / den;
}

void FuzzyController::build()
{
    for (int i = 0; i < N; ++i)
    {
        assert((m_ranges[i].size() > 0) && "At least one of the ranges was not set.");
        assert((m_mfs[i].size() > 0U) && "MFs for at least one variable were not set.");
        assert((m_rules[i].size() > 0U) && "Rules for at least one variable were not set.");
    }

    m_table = new float*[m_ranges[0].size()];
    for (int i = 0; i < m_ranges[0].size(); ++i)
    {
        m_table[i] = new float[m_ranges[1].size()];
    }

    for (int i = 0; i < m_ranges[0].size(); ++i)
    {
        for (int j = 0; j < m_ranges[1].size(); ++j)
        {
            m_table[i][j] = centroid(i, j);
        }
    }
}

float FuzzyController::calculate(float v1, float v2)
{
    const int v1Idx = m_ranges[0].pre(v1);
    const int v2Idx = m_ranges[1].pre(v2);

    const float den = fabsf(m_ranges[0].at(v1Idx) - m_ranges[0].at(v1Idx + 1)) *
        fabsf(m_ranges[1].at(v2Idx) - m_ranges[1].at(v2Idx + 1));

    const float w11 = (fabsf(v1 - m_ranges[0].at(v1Idx + 1)) * fabsf(v2 - m_ranges[1].at(v2Idx + 1))) / den;
    const float w12 = (fabsf(v1 - m_ranges[0].at(v1Idx)) * fabsf(v2 - m_ranges[1].at(v2Idx + 1))) / den;
    const float w21 = (fabsf(v1 - m_ranges[0].at(v1Idx + 1)) * fabsf(v2 - m_ranges[1].at(v2Idx))) / den;
    const float w22 = (fabsf(v1 - m_ranges[0].at(v1Idx)) * fabsf(v2 - m_ranges[1].at(v2Idx))) / den;

    const float x11 = m_table[v1Idx][v2Idx];
    const float x12 = m_table[v1Idx][v2Idx + 1];
    const float x21 = m_table[v1Idx + 1][v2Idx];
    const float x22 = m_table[v1Idx + 1][v2Idx + 1];

    return (w11 * x11) + (w12 * x12) + (w21 * x21) + (w22 * x22);
}
