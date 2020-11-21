#include <iostream>
#include "ResonanceType.h"
ResonanceType::ResonanceType(const char *Name, const double Mass, const int Charge, const double Width)
    : ParticleType(Name, Mass, Charge), fWidth(Width) {}
double ResonanceType::GetWidth() { return fWidth; }
void ResonanceType::Print() const
{
    ParticleType::Print();
    std::cout << "Width: " << fWidth << '\n'
              << '\n';
}
