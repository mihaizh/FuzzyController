#ifndef FUZZY_C_H
#define FUZZY_C_H

#include "fuzzy_mfs.h"
#include <vector>

using MFVector = std::vector<sBasicMF>;
using Rule = std::vector<int>;
using RuleVector = std::vector<Rule>;

class Range
{
private:
    float m_min;
    float m_max;
    float m_step;

    float* m_data;
    int m_size;

public:
    Range();
    Range(float min, float max, float step);
    Range(const Range& other);
    Range& operator=(const Range& other);

    float at(int idx) const
    {
        return m_data[idx];
    }

    int size() const
    {
        return m_size;
    }

    int pre(float val);
};

class FuzzyController
{
private:
    static const int N = 3; // 2 input + 1 output

    Range m_ranges[N];
    MFVector m_mfs[N];
    RuleVector m_rules;

    float** m_table;
    float m_coeffs[N];

    float centroid(int i, int j);

public:
    FuzzyController();

    void setRange(int varIdx, const Range& range)
    {
        m_ranges[varIdx] = range;
    }

    void setMFs(int varIdx, const MFVector& mfs)
    {
        m_mfs[varIdx] = mfs;
    }

    void setRules(const RuleVector& rules)
    {
        m_rules = rules;
    }

    void setCoeffs(const float coeffs[N])
    {
        memcpy(m_coeffs, coeffs, sizeof(m_coeffs));
    }

    void build();
    float calculate(float v1, float v2);
};

#endif // FUZZY_C_H
