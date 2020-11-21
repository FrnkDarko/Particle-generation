#ifndef RESONANCETYPE_H
#define RESONANCETYPE_H
#include "ParticleType.h"
class ResonanceType : public ParticleType
{
public:
    ResonanceType(const char *Name, const double Mass, const int Charge, const double Width);
    double GetWidth();
    void Print() const;

private:
    const double fWidth;
};
#endif