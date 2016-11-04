#ifndef BORDERTEST2_H
#define BORDERTEST2_H


class BorderTest2
{
public:
    BorderTest2();
    virtual ~BorderTest2() {}

    double f(int i, int j) const;
    double phi1(int i, int j) const;
    double phi2(int i, int j) const;

    double initial(int i) const;
    double border(int type, int j) const;
};

#endif // BORDERTEST2_H
