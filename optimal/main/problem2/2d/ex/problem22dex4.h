#ifndef PROBLEM22DEX4_H
#define PROBLEM22DEX4_H

#include "../abs/abstractproblem22d.h"
#include <imaging.h>

//---------------------------------------------------------------------------------------------------------------//
// Məqalə üçün hazırlanmış test məsələ
// İdarə nöqtələrinin sayı 2.
// Ölçü cihazlarının sayı 3.
//---------------------------------------------------------------------------------------------------------------//

class Problem22DEx4 : public AbstactProblem22D
{
public:
    static void Main(int argc, char* argv[]);

    Problem22DEx4();
    virtual ~Problem22DEx4();

    void compareNandAGradients();
};

#endif // PROBLEM22DEX4_H
