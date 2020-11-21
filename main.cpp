#include "ParticleType.h"
#include "ResonanceType.h"
#include "Particle.h"
#include <iostream>
#include <cmath>
#include <TRandom.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TFile.h>

int main()
{
    Particle::AddParticleType("Pi+", 0.13957, +1, 0);
    Particle::AddParticleType("Pi-", 0.13957, -1, 0);
    Particle::AddParticleType("K+", 0.49367, +1, 0);
    Particle::AddParticleType("K-", 0.49367, -1, 0);
    Particle::AddParticleType("P+", 0.93827, +1, 0);
    Particle::AddParticleType("P-", 0.93827, -1, 0);
    Particle::AddParticleType("K*", 0.89166, 0, 0.050);

    TH1F *particle_type = new TH1F("particle_type", "Types of particles generated", 7, 0., 7.);
    TH1F *phi_angle = new TH1F("phi_angle", "Azimutal angle", 100, 0., 2 * M_PI);
    TH1F *theta_angle = new TH1F("theta_angle", "Polar angle", 100, 0., M_PI);
    TH1F *impulse = new TH1F("impulse", "Impulse", 100, 0., 7.);
    TH1F *t_impulse = new TH1F("t_impulse", "Trasversal impulse", 100, 0., 5.);
    TH1F *energy = new TH1F("energy", "Energy", 100, 0., 7.);
    TH1F *invmass_tot = new TH1F("invmass_tot", "Total invariant mass", 100, 0., 2.);
    TH1F *invmass_disc = new TH1F("invmass_disc", "Discord invariant mass", 80, 0., 2.);
    TH1F *invmass_conc = new TH1F("invmass_conc", "Concord invariant mass", 80, 0., 2.);
    TH1F *invmass_pi_k_disc = new TH1F("invmass_pi_k_disc", "Discord Pion and Kaon invariant mass", 80, 0., 2.);
    TH1F *invmass_pi_k_conc = new TH1F("invmass_pi_k_conc", "Concord Pion and Kaon invariant mass", 80, 0., 2.);
    TH1F *invmass_k = new TH1F("invmass_k", "Invariant mass of decayed particles", 80, 0., 2.);

    gRandom->SetSeed();
    int N = 120;
    Particle particle[N];
    for (int i = 0; i != 1E5; i++)
    {
        for (int j = 0; j != 100; j++)
        {
            Particle part;
            double phi = gRandom->Uniform(0, 2 * M_PI);
            phi_angle->Fill(phi);
            double theta = gRandom->Uniform(0, M_PI);
            theta_angle->Fill(theta);
            double P = gRandom->Exp(1);
            impulse->Fill(P);
            double px = P * sin(theta) * cos(phi);
            double py = P * sin(theta) * sin(phi);
            double pz = P * cos(theta);
            double pt = sqrt(pow(px, 2) + pow(px, 2));
            t_impulse->Fill(pt);
            part.SetP(px, py, pz);
            int k = 0;
            double x = gRandom->Uniform(0, 1);
            if (x < 0.8)
            {
                double y = gRandom->Uniform(0, 1);
                if (y < 0.5)
                {
                    part.SetIndex("Pi+");
                }
                else
                {
                    part.SetIndex("Pi-");
                }
            }
            if (x > 0.8 && x < 0.9)
            {
                double y = gRandom->Uniform(0, 1);
                if (y < 0.5)
                {
                    part.SetIndex("K+");
                }
                else
                {
                    part.SetIndex("K-");
                }
            }
            if (x > 0.9 && x < 0.99)
            {
                double y = gRandom->Uniform(0, 1);
                if (y < 0.5)
                {
                    part.SetIndex("P+");
                }
                else
                {
                    part.SetIndex("P-");
                }
            }
            else if (x > 0.99)
            {
                part.SetIndex("K*");
                double y = gRandom->Uniform(0, 1);
                Particle dau1;
                Particle dau2;
                if (y < 0.5)
                {
                    dau1.SetIndex("Pi+");
                    dau2.SetIndex("K-");
                }
                else
                {
                    dau1.SetIndex("Pi-");
                    dau2.SetIndex("K+");
                }
                part.Decay2Body(dau1, dau2);
                invmass_k->Fill(dau1.InvMass(dau2));
                particle[100 + k] = dau1;
                k++;
                particle[100 + k] = dau2;
                k++;
            }
            particle_type->Fill(part.GetIndex());
            particle[j] = part;
            energy->Fill(part.GetEnergy());

            if (j > 0)
            {
                invmass_tot->Fill(particle[0].InvMass(particle[j]));
                if (particle[0].GetIndex() % 2 == 0)
                {
                    if (particle[j].GetIndex() % 2 == 0)
                    {
                        invmass_conc->Fill(particle[0].InvMass(particle[j]));
                    }
                    else
                    {
                        invmass_disc->Fill(particle[0].InvMass(particle[j]));
                    }
                }
                if (particle[0].GetIndex() % 2 != 0)
                {
                    if (particle[j].GetIndex() % 2 != 0)
                    {
                        invmass_conc->Fill(particle[0].InvMass(particle[j]));
                    }
                    else
                    {
                        invmass_disc->Fill(particle[0].InvMass(particle[j]));
                    }
                }
            }
            for (int m = 0; m != j; m++)
            {
                if (particle[m].GetIndex() == 0 && particle[j].GetIndex() == 3)
                {
                    invmass_pi_k_disc->Fill(particle[m].InvMass(particle[j]));
                }
                else if (particle[m].GetIndex() == 1 && particle[j].GetIndex() == 2)
                {
                    invmass_pi_k_disc->Fill(particle[m].InvMass(particle[j]));
                }
                else if (particle[m].GetIndex() == 0 && particle[j].GetIndex() == 2)
                {
                    invmass_pi_k_conc->Fill(particle[m].InvMass(particle[j]));
                }
                else if (particle[m].GetIndex() == 1 && particle[j].GetIndex() == 3)
                {
                    invmass_pi_k_conc->Fill(particle[m].InvMass(particle[j]));
                }
            }
        }
    }
    TCanvas *c1 = new TCanvas("c1", "Particle types");
    c1->Divide(3, 2);
    c1->cd(1);
    particle_type->Draw();
    c1->cd(2);
    phi_angle->Draw();
    c1->cd(3);
    theta_angle->Draw();
    c1->cd(4);
    impulse->Draw();
    c1->cd(5);
    t_impulse->Draw();
    c1->cd(6);
    energy->Draw();

    TCanvas *c2 = new TCanvas("c2", "Invariant mass");
    c2->Divide(3, 2);
    c2->cd(1);
    invmass_tot->Draw();
    c2->cd(2);
    invmass_disc->Draw();
    c2->cd(3);
    invmass_conc->Draw();
    c2->cd(4);
    invmass_pi_k_disc->Draw();
    c2->cd(5);
    invmass_pi_k_conc->Draw();
    c2->cd(6);
    invmass_k->Draw();

    TFile *file = new TFile("particles.root", "RECREATE");
    particle_type->Write();
    phi_angle->Write();
    theta_angle->Write();
    impulse->Write();
    energy->Write();
    invmass_tot->Write();
    invmass_disc->Write();
    invmass_conc->Write();
    invmass_pi_k_disc->Write();
    invmass_pi_k_conc->Write();
    invmass_k->Write();
    file->Close();
}