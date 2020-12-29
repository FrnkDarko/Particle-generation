#include "Particle.h"
#include <iostream>
#include <cmath>
#include <cstdlib>

Particle::Particle() : fIndex{0}, fPx{0.}, fPy{0.}, fPz{0.} {}

Particle::Particle(const char *Name, double Px, double Py, double Pz)
    : fName(Name), fPx(Px), fPy(Py), fPz(Pz) { fIndex = FindParticle(Name); }
int Particle::fNParticleType = 0;
ParticleType *Particle::fParticleType[fMaxNum];
int Particle::GetIndex() const { return fIndex; }
double Particle::GetPx() const { return fPx; }
double Particle::GetPy() const { return fPy; }
double Particle::GetPz() const { return fPz; }
void Particle::SetIndex(int const Index) { fIndex = Index; }
void Particle::SetIndex(const char *name)
{
    fIndex = FindParticle(name);
}

void Particle::SetP(double Px, double Py, double Pz)
{
    fPx = Px;
    fPy = Py;
    fPz = Pz;
}

double Particle::GetMass() const
{
    return (fParticleType[fIndex]->GetMass());
}

double Particle::GetEnergy() const
{
    double m = pow(GetMass(), 2);
    double p = pow(fPx, 2) + pow(fPy, 2) + pow(fPz, 2);
    return sqrt(m + p);
}

double Particle::InvMass(Particle &p) const
{
    double E2 = (GetEnergy() + p.GetEnergy()) * (GetEnergy() + p.GetEnergy());
    double Px = fPx + p.GetPx();
    double Py = fPy + p.GetPy();
    double Pz = fPz + p.GetPz();
    double P2 = pow(Px, 2) + pow(Py, 2) + pow(Pz, 2);
    return sqrt(E2 - P2);
}

void Particle::PrintIndex()
{
    for (int i = 0; i < fNParticleType; i++)
    {
        fParticleType[i]->Print();
    }
}

void Particle::PrintParticle()
{
    std::cout << "Type: " << FindParticle(fName) << '\n';
    std::cout << "Name: " << fName << '\n';
    std::cout << "Momentum: (" << fPx << ", " << fPy << ", " << fPz << ")" << '\n';
}

int Particle::FindParticle(const char *name)
{
    int f;
    int result;
    int i = 0;
    for (int j = 0; j < fNParticleType; j++)
    {
        if (fParticleType[j]->GetName() == name)
        {
            f = 1;
            break;
        }
        else
        {
            i++;
        }
    }
    if (f == 1)
    {
        result = i;
    }
    else
    {
        result = fMaxNum + 1;
    }
    return result;
}

void Particle::AddParticleType(const char *name, const double mass, const int charge, const double width)
{
    int FindIndex = FindParticle(name);
    if (FindIndex == fMaxNum + 1)
    {
        if (fNParticleType < fMaxNum)
        {
            if (width == 0)
            {
                ParticleType *p = new ParticleType(name, mass, charge);
                fParticleType[fNParticleType] = p;
                ++fNParticleType;
            }
            else
            {
                ResonanceType *r = new ResonanceType(name, mass, charge, width);
                fParticleType[fNParticleType] = r;
                ++fNParticleType;
            }
        }
        else
        {
            std::cout << "Cannot add the particle, the index is full" << '\n';
        }
    }
    else
    {
        std::cout << "The particle type is already listed in the index" << '\n';
    }
}

void Particle::Boost(double bx, double by, double bz)
{

    double energy = GetEnergy();

    double b2 = bx * bx + by * by + bz * bz;
    double gamma = 1.0 / sqrt(1.0 - b2);
    double bp = bx * fPx + by * fPy + bz * fPz;
    double gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;

    fPx += gamma2 * bp * bx + gamma * bx * energy;
    fPy += gamma2 * bp * by + gamma * by * energy;
    fPz += gamma2 * bp * bz + gamma * bz * energy;
}

int Particle::Decay2Body(Particle &dau1, Particle &dau2) const
{
    if (GetMass() == 0.0)
    {
        std::cout << "Decayment cannot be preformed if mass is zero" << '\n';
        return 1;
    }

    double massMot = GetMass();
    double massDau1 = dau1.GetMass();
    double massDau2 = dau2.GetMass();

    if (fIndex > -1)
    {
        float x1, x2, w, y1, y2;

        double invnum = 1. / RAND_MAX;
        do
        {
            x1 = 2.0 * rand() * invnum - 1.0;
            x2 = 2.0 * rand() * invnum - 1.0;
            w = x1 * x1 + x2 * x2;
        } while (w >= 1.0);

        w = sqrt((-2.0 * log(w)) / w);
        y1 = x1 * w;
        y2 = x2 * w;

        massMot += fParticleType[fIndex]->GetWidth() * y1;
    }

    if (massMot < massDau1 + massDau2)
    {
        std::cout << "Decayment cannot be preformed because mass is too low in this channel" << '\n';
        return 2;
    }

    double pout = sqrt((massMot * massMot - (massDau1 + massDau2) * (massDau1 + massDau2)) * (massMot * massMot - (massDau1 - massDau2) * (massDau1 - massDau2))) / massMot * 0.5;

    double norm = 2 * M_PI / RAND_MAX;

    double phi = rand() * norm;
    double theta = rand() * norm * 0.5 - M_PI / 2.;
    dau1.SetP(pout * sin(theta) * cos(phi), pout * sin(theta) * sin(phi), pout * cos(theta));
    dau2.SetP(-pout * sin(theta) * cos(phi), -pout * sin(theta) * sin(phi), -pout * cos(theta));

    double energy = sqrt(fPx * fPx + fPy * fPy + fPz * fPz + massMot * massMot);

    double bx = fPx / energy;
    double by = fPy / energy;
    double bz = fPz / energy;

    dau1.Boost(bx, by, bz);
    dau2.Boost(bx, by, bz);

    return 0;
}
