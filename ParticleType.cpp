#include <iostream>
#include "ParticleType.h"
ParticleType::ParticleType(const char *Name, const double Mass, const int Charge)
    : fName(Name), fMass(Mass), fCharge(Charge) {}
const char *ParticleType::GetName() const { return fName; }
double ParticleType::GetMass() const { return fMass; }
int ParticleType::GetCharge() const { return fCharge; }
int ParticleType::GetWidth() const { return 0; }
void ParticleType::Print() const
{
    std::cout << "Name: " << fName << '\n';
    std::cout << "Mass: " << fMass << '\n';
    std::cout << "Charge: " << fCharge << '\n' << '\n';
}
