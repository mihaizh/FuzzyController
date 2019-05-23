#ifndef FUZZY_C_H
#define FUZZY_C_H

#include "fuzzy_mfs.h"
#include <vector>

using MFVector = std::vector<sBasicMF>;

template <typename RuleType_t>
using Rule = std::vector<RuleType_t>;

template <typename RuleType_t>
using RuleVector = std::vector<Rule<RuleType_t>>;

template <typename TableType_t>
using Table = std::vector<std::vector<TableType_t>>;

class Range
{
private:
    float m_min;
    float m_max;
    float m_step;

    std::vector<float> m_data;

public:
    Range();
    Range(float min, float max, float step);
    Range(float min, float max, size_t N);
    Range(const Range& other);
    Range(Range&& other) = default;
    Range& operator=(const Range& other);
    Range& operator=(Range&& other) = default;

    float at(size_t idx) const
    {
        return m_data[idx];
    }

    size_t size() const
    {
        return m_data.size();
    }

    size_t pre(float val);
};

template <typename RuleType_t>
class FuzzyController
{
private:
    static const size_t N = 3U; // 2 input + 1 output

    Range m_ranges[N];
    MFVector m_mfs[N];
    RuleVector<RuleType_t> m_rules;
    Table<float> m_table;

    std::vector<float> m_coeffs;

    float centroid(int i, int j);
    float center(
        const Table<float>& minHeight,
        const Table<size_t>& memFuncIdx);

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

    void setRules(const RuleVector<RuleType_t>& rules)
    {
        m_rules = rules;
    }

    void setCoeffs(const std::vector<float>& coeffs)
    {
        assert((coeffs.size() <= N) &&
            "Number of scaling coefficients has to be lower than number of variables.");

        std::copy(coeffs.begin(), coeffs.end(), m_coeffs.begin());
    }

    void build();
    float calculate(float v1, float v2);
};

using MamdaniFuzzyController = FuzzyController<size_t>;
using SugenoFuzzyController = FuzzyController<float>;

using MamdaniRuleVector = RuleVector<size_t>;
using SugenoRuleVector = RuleVector<float>;

#endif // FUZZY_C_H
