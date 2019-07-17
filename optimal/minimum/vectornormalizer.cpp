#include "vectornormalizer.h"
#include <math.h>
#include <float.h>

auto IVectorNormalizer::EuclideanNorm(const DoubleVector &v) -> double
{
    return v.EuclideanNorm();
}

auto IVectorNormalizer::EuclideanNormalize(DoubleVector &v) -> void
{
    v.EuclideanNormalize();
}

auto IVectorNormalizer::L1Norm(const DoubleVector &v) -> double
{
    return v.L1Norm();
}

auto IVectorNormalizer::L1Normalize(DoubleVector &v) -> void
{
    v.L1Normalize();
}

auto IVectorNormalizer::L2Norm(const DoubleVector &v) -> double
{
    return v.L2Norm();
}

auto IVectorNormalizer::L2Normalize(DoubleVector &v) -> void
{
    v.L2Normalize();
}

auto IVectorNormalizer::L2FuncNorm(const DoubleVector &v, double h) -> double
{
    if (v.length() < 2) return 0.0;

    double norm = 0.0;
    auto length = v.length();

    norm = v.mData[0] * v.mData[0];
    for (unsigned int i=1; i < length-1; i++)
    {
        norm += 2.0 * v.mData[i] * v.mData[i];
    }
    norm += v.mData[length-1] * v.mData[length-1];
    norm *= 0.5*h;

    return sqrt(norm);
}

auto IVectorNormalizer::L2FuncNormalize(DoubleVector &v, double h) -> void
{
    if (v.length() < 2) return;
    auto length = v.length();
    auto norm = L2FuncNorm(v, h);
    if (norm > DBL_EPSILON) for (unsigned int i=0; i<length; i++) v.mData[i] /= norm;
}

auto IVectorNormalizer::LpNorm(const DoubleVector &v, double p) -> double
{
    if (v.empty()) return 0.0;
    double norm = 0.0;
    auto length = v.length();
    for (unsigned int i=0; i < length; i++)
    {
        norm += pow(v.mData[i], p);
    }
    return pow(norm, 1.0/p);
}

auto IVectorNormalizer::LpNormalize(DoubleVector &v, double p) -> void
{
    if (v.empty()) return;
    auto length = v.length();
    auto norm = LpNorm(v, p);
    if (norm > DBL_EPSILON) for (unsigned int i=0; i<length; i++) v.mData[i] /= norm;
}

auto IVectorNormalizer::LInfNorm(const DoubleVector &v) -> double
{
    if (v.empty()) return 0.0;
    auto length = v.length();
    auto norm = fabs(v.mData[0]);
    for (unsigned int i=1; i < length; i++)
    {
        if (norm < fabs(v.mData[i])) norm = fabs(v.mData[i]);
    }
    return norm;
}

auto IVectorNormalizer::LInfNormalize(DoubleVector &v) -> void
{
    if (v.empty()) return;
    auto length = v.length();
    auto norm = LInfNorm(v);
    if (norm > DBL_EPSILON) for (unsigned int i=0; i<length; i++) v.mData[i] /= norm;
}
