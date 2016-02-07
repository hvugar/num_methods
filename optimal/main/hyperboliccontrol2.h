#ifndef HYPERBOLICCONTROL2_H
#define HYPERBOLICCONTROL2_H

#include <function.h>
#include <hyperbolicequation.h>
#include <printer.h>

class HyperbolicControl2 : public RnFunction, public IGradient, public IHyperbolicEquation, public IBackwardHyperbolicEquation, public IPrinter
{
public:
    HyperbolicControl2();
    ~HyperbolicControl2() {}
};

#endif // HYPERBOLICCONTROL2_H
