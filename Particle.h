#ifndef PARTICLE_H
#define PARTICLE_H
#include "ParticleType.h"
#include "ResonanceType.h"

class Particle
{
public:
    Particle();
    Particle(const char *Name, double Px, double Py, double Pz);
    int GetIndex() const;
    double GetPx() const;
    double GetPy() const;
    double GetPz() const;
    void SetIndex(int const Index);
    void SetIndex(const char *name);
    void SetP(double Px, double Py, double Pz);
    double GetMass() const;
    double GetEnergy() const;
    double InvMass(Particle &p) const;
    static void PrintIndex();
    void PrintParticle();
    static void AddParticleType(const char *name, const double mass, const int charge, const double width);
    int Decay2Body(Particle &dau1, Particle &dau2) const;

private:
    const char *fName;
    static int const fMaxNum = 10;
    static ParticleType *fParticleType[fMaxNum];
    static int fNParticleType;
    int fIndex;
    double fPx;
    double fPy;
    double fPz;
    static int FindParticle(const char *name);
    void Boost(double bx, double by, double bz);
};
#endif