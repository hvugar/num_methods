#include <iostream>

using namespace std;

#include <grid/bvp.h>
#include <ode/lode1o.h>


class FirstOrderLinearODEY01 : public FirstOrderLinearODE
{};

class FirstOrderLinearODEY02 : public FirstOrderLinearODE
{};

class FirstOrderLinearODEY11 : public FirstOrderLinearODE
{};

class FirstOrderLinearODEY12 : public FirstOrderLinearODE
{};



int main()
{
    FirstOrderLinearODEY01 f1;

    cout << "Hello World!" << endl;
    return 0;
}
