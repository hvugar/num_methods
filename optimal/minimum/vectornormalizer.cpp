#include "vectornormalizer.h"
#include <math.h>

auto IVectorNormalizer::EuclideanNorm(const DoubleVector &v) -> double
{
    if (v.empty()) return 0.0;
    auto norm = 0.0;
    auto length = v.length();
    for (unsigned int i=0; i < length; i++)
    {
        double item = v.mData[i];
        norm += item * item;
    }
    return sqrt(norm);
}

auto IVectorNormalizer::EuclideanNormalize(DoubleVector &v) -> void
{
    if (v.empty()) return;
    auto length = v.length();
    auto norm = EuclideanNorm(v);
    for (unsigned int i=0; i<length; i++) v.mData[i] /= norm;
}

auto IVectorNormalizer::L1Norm(const DoubleVector &v) -> double
{
    if (v.empty()) return 0.0;
    auto norm = 0.0;
    auto length = v.length();
    for (unsigned int i=0; i < length; i++)
    {
        norm += fabs(v.mData[i]);
    }
    return norm;
}

auto IVectorNormalizer::L1Normalize(DoubleVector &v) -> void
{
    if (v.empty()) return;
    auto length = v.length();
    auto norm = L1Norm(v);
    for (unsigned int i=0; i<length; i++) v.mData[i] /= norm;
}

auto IVectorNormalizer::L2Norm(const DoubleVector &v) -> double
{
    if (v.empty()) return 0.0;
    double norm = 0.0;
    auto length = v.length();
    for (unsigned int i=0; i < length; i++)
    {
        norm += v.mData[i] * v.mData[i];
    }
    return sqrt(norm);
}

auto IVectorNormalizer::L2Normalize(DoubleVector &v) -> void
{
    if (v.empty()) return;
    auto length = v.length();
    auto norm = L2Norm(v);
    for (unsigned int i=0; i<length; i++) v.mData[i] /= norm;
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
    for (unsigned int i=0; i<length; i++) v.mData[i] /= norm;
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
    for (unsigned int i=0; i<length; i++) v.mData[i] /= norm;
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
    for (unsigned int i=0; i<length; i++) v.mData[i] /= norm;
}
