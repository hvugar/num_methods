#ifndef VECTORNORMALIZER_H
#define VECTORNORMALIZER_H

#include "global.h"
#include "vector2d.h"

class MINIMUMSHARED_EXPORT IVectorNormalizer
{
public:
    enum Norm
    {
        EUCLIDEAN_NORM,
        L1_NORM,
        L2_NORM,
        L2_FUNC_NORM,
        Lp_NORM,
        LInf_NORM,
        NOT_NORM
    };

public:
    virtual auto norm(const DoubleVector &v) const -> double = 0;
    virtual auto normalize(DoubleVector &v) const -> void = 0;

public:
    static auto EuclideanNorm(const DoubleVector &v) -> double;
    static auto EuclideanNormalize(DoubleVector &v) -> void;
    static auto L1Norm(const DoubleVector &v) -> double;
    static auto L1Normalize(DoubleVector &v) -> void;
    static auto L2Norm(const DoubleVector &v) -> double;
    static auto L2Normalize(DoubleVector &v) -> void;
    static auto L2FuncNorm(const DoubleVector &v, double h) -> double;
    static auto L2FuncNormalize(DoubleVector &v, double h) -> void;
    static auto LpNorm(const DoubleVector &v, double p) -> double;
    static auto LpNormalize(DoubleVector &v, double p) -> void;
    static auto LInfNorm(const DoubleVector &v) -> double;
    static auto LInfNormalize(DoubleVector &v) -> void;

    friend class DoubleVector;
};

#endif // VECTORNORMALIZER_H
