#ifndef FUZZY_MFS_H
#define FUZZY_MFS_H

#include <assert.h>

enum MF_TYPE
{
    MF_NONE,

    MF_TRI,
    MF_TRAP,

    MF_NUM
};

struct sBasicMF // trimf, trapmf
{
    sBasicMF(MF_TYPE t, float _a, float _b, float _c, float _d)
        : type(t),
        a(_a),
        b(_b),
        c(_c),
        d(_d)
    {
    }

    MF_TYPE type;

    float a;
    float b;
    float c;
    float d;
};

inline float getH(const sBasicMF& mf, float x)
{
    float val = 0.F;

    switch (mf.type)
    {
        case MF_TRI:
        {
            if ((x >= mf.a) && (x <= mf.c))
            {
                if (x <= mf.b)
                {
                    val = (x - mf.a) / (mf.b - mf.a);
                }
                else
                {
                    val = (mf.c - x) / (mf.c - mf.b);
                }
            }
        }
        break;
        case MF_TRAP:
        {
            if ((x >= mf.a) && (x <= mf.d))
            {
                if (x <= mf.b)
                {
                    val = (x - mf.a) / (mf.b - mf.a);
                }
                else if (x <= mf.c)
                {
                    val = 1.F;
                }
                else
                {
                    val = (mf.d - x) / (mf.d - mf.c);
                }
            }
        }
        break;
        default:
        {
            assert(false && "Unknown MF type!");
        }
        break;
    }

    return val;
}

inline float getC(const sBasicMF& mf, float h)
{
    float val = 0.F;

    switch (mf.type)
    {
        case MF_TRI:
        {
            const float v1 = (-h * (mf.c - mf.b)) + mf.c;
            const float v2 = (h * (mf.b - mf.a)) + mf.a;
            val = (v1 + v2) * 0.5F;
        }
        break;
        case MF_TRAP: // ??
        {
            const float v1 = (-h * (mf.d - mf.c)) + mf.d;
            const float v2 = (h * (mf.b - mf.a)) + mf.a;
            val = (v1 + v2) * 0.5F;
        }
        break;
        default:
        {
            assert(false && "Unknown MF type!");
        }
        break;
    }

    return val;
}

#endif // FUZZY_MFS_H
