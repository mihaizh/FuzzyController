#include "FuzzyController.h"

#include <math.h>
#include <assert.h>

Range::Range()
    : m_min(0.F),
    m_max(0.F),
    m_step(0.F)
{
}

Range::Range(float min, float max, float step)
    : m_min(min),
    m_max(max),
    m_step(step)
{
    assert((min < max) && "min should be smaller than max");

    const int size = (int)floorf((fabsf(min) + fabs(max)) / step);
    m_data.resize(static_cast<size_t>(size));

    for (size_t i = 0U; i < m_data.size(); ++i)
    {
        m_data[i] = min + (i * step);
    }
}

Range::Range(float min, float max, size_t N)
    : m_min(min),
    m_max(max),
    m_step(0.F)
{
    assert((min < max) && "min should be smaller than max");

    m_step = ceilf((fabsf(min) + fabs(max)) / N);
    m_data.resize(N);

    for (size_t i = 0; i < m_data.size(); ++i)
    {
        m_data[i] = min + (i * m_step);
    }
}

Range::Range(const Range& other)
    : m_min(other.m_min),
    m_max(other.m_max),
    m_step(other.m_step),
    m_data(other.m_data)
{
}

Range& Range::operator=(const Range& other)
{
    if (this != &other)
    {
        m_min = other.m_min;
        m_max = other.m_max;
        m_step = other.m_step;

        std::copy(
            other.m_data.begin(),
            other.m_data.end(),
            std::back_inserter(m_data));
    }

    return *this;
}

size_t Range::pre(float val)
{
    if (val < m_min)
    {
        return 0;
    }

    if (val > m_max)
    {
        return (m_data.size() - 1U);
    }

    const float pos = floorf((val - m_min) / m_step);
    return static_cast<int>(pos);
}

template <typename RuleType_t>
FuzzyController<RuleType_t>::FuzzyController()
    : m_coeffs(N)
{
    std::fill(m_coeffs.begin(), m_coeffs.end(), 1.F);
}

template <>
float FuzzyController<size_t>::center( // Mamdani
    const Table<float>& minHeight,
    const Table<size_t>& memFuncIdx)
{
    float num = 0.F;
    float den = 0.F;
    for (size_t k = 0; k < minHeight.size(); ++k)
    {
        for (size_t l = 0; l < minHeight[k].size(); ++l)
        {
            const size_t mf1 = memFuncIdx[0][k];
            const size_t mf2 = memFuncIdx[1][l];
            const size_t ruleIdx = m_rules[mf1][mf2];

            const float c = getC(m_mfs[2][ruleIdx], minHeight[k][l]);

            num += c * minHeight[k][l];
            den += minHeight[k][l];
        }
    }

    return num / den;
}

template <>
float FuzzyController<float>::center( // Sugeno-Takagi
    const Table<float>& minHeight,
    const Table<size_t>& memFuncIdx)
{
    float num = 0.F;
    float den = 0.F;
    for (size_t k = 0; k < minHeight.size(); ++k)
    {
        for (size_t l = 0; l < minHeight[k].size(); ++l)
        {
            const size_t mf1 = memFuncIdx[0][k];
            const size_t mf2 = memFuncIdx[1][l];

            const float c = m_rules[mf1][mf2];

            num += c * minHeight[k][l];
            den += minHeight[k][l];
        }
    }

    return num / den;
}

template <typename RuleType_t>
float FuzzyController<RuleType_t>::centroid(int i, int j)
{
    const size_t INPUTS = N - 1;

    Table<float> height(INPUTS);
    Table<size_t> memFuncIdx(INPUTS);

    // compute heights
    for (size_t k = 0; k < INPUTS; ++k)
    {
        const size_t ruleIdx = (k == 0) ? i : j;
        for (size_t l = 0; l < m_mfs[k].size(); ++l)
        {
            const float v = m_ranges[k].at(ruleIdx);
            const float h = getH(m_mfs[k][l], v);

            if (h > std::numeric_limits<float>::epsilon())
            {
                height[k].push_back(h);
                memFuncIdx[k].push_back(l);
            }
        }
    }

    // apply "min" operator
    Table<float> minHeight(height[0].size());
    for (size_t k = 0; k < height[0].size(); ++k)
    {
        minHeight[k].resize(height[1].size());
        for (size_t l = 0; l < height[1].size(); ++l)
        {
            if (height[0][k] < height[1][l])
            {
                minHeight[k][l] = height[0][k];
            }
            else
            {
                minHeight[k][l] = height[1][l];
            }
        }
    }

    // calculate the value using the centroid formula
    // sum(center * height) / sum(heights)
    return center(minHeight, memFuncIdx);
}

template <typename RuleType_t>
void FuzzyController<RuleType_t>::build()
{
    for (int i = 0; i < N; ++i)
    {
        assert((m_ranges[i].size() > 0) && "At least one of the ranges was not set.");
        assert((m_mfs[i].size() > 0) && "MFs for at least one variable were not set.");

        if (std::is_same<RuleType_t, size_t>::value) // Mamdani
        {
            assert((m_rules[i].size() > 0) && "Rules for at least one variable were not set.");
        }
    }

    m_table.resize(m_ranges[0].size());
    for (auto& v : m_table)
    {
        v.resize(m_ranges[1].size());
    }

    for (int i = 0; i < m_ranges[0].size(); ++i)
    {
        for (int j = 0; j < m_ranges[1].size(); ++j)
        {
            m_table[i][j] = centroid(i, j);
        }
    }
}

template <typename RuleType_t>
float FuzzyController<RuleType_t>::calculate(float v1, float v2)
{
    // scale inputs
    v1 *= m_coeffs[0];
    v2 *= m_coeffs[1];

    // find lower bound index for value
    // in the range
    const size_t v1Idx = m_ranges[0].pre(v1);
    const size_t v2Idx = m_ranges[1].pre(v2);

    // compute the weight denominator
    // this could be pre-computed for uniform ranges
    const float den = fabsf(m_ranges[0].at(v1Idx) - m_ranges[0].at(v1Idx + 1)) *
        fabsf(m_ranges[1].at(v2Idx) - m_ranges[1].at(v2Idx + 1));

    // compute weights
    const float w11 = (fabsf(v1 - m_ranges[0].at(v1Idx + 1)) * fabsf(v2 - m_ranges[1].at(v2Idx + 1))) / den;
    const float w12 = (fabsf(v1 - m_ranges[0].at(v1Idx))     * fabsf(v2 - m_ranges[1].at(v2Idx + 1))) / den;
    const float w21 = (fabsf(v1 - m_ranges[0].at(v1Idx + 1)) * fabsf(v2 - m_ranges[1].at(v2Idx)))     / den;
    const float w22 = (fabsf(v1 - m_ranges[0].at(v1Idx))     * fabsf(v2 - m_ranges[1].at(v2Idx)))     / den;

    // values from the Mamdani table
    const float x11 = m_table[v1Idx][v2Idx];
    const float x12 = m_table[v1Idx][v2Idx + 1];
    const float x21 = m_table[v1Idx + 1][v2Idx];
    const float x22 = m_table[v1Idx + 1][v2Idx + 1];

    // weighted mean, then scaled
    return ((w11 * x11) + (w12 * x12) + (w21 * x21) + (w22 * x22)) * m_coeffs[2];
}

template class FuzzyController<size_t>; // Mamdani implementation
template class FuzzyController<float>;  // Sugeno implementation

