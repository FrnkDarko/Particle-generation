#include <TFile.h>
#include <TH1F.h>
#include <iostream>
#include <TCanvas.h>

void analysis(){
    TCanvas *c1 = new TCanvas("particles", "Particles generated");
    TFile *histo = new TFile("Particles.root");
    TH1F *particle_type = (TH1F*) histo->Get("particle_type");
    c1->Divide(3,2);
    c1->cd(1);
    particle_type->Draw();
    
    std::cout << "Number of particles generated: " << particle_type->GetEntries() << '\n' << '\n';
    double pions = particle_type->GetBinContent(1) + particle_type->GetBinContent(2);
    double e_pions = particle_type->GetBinError(1) + particle_type->GetBinError(2);
    double kaons = particle_type->GetBinContent(3) + particle_type->GetBinContent(4);
    double e_kaons = particle_type->GetBinError(3) + particle_type->GetBinError(4);
    double protons = particle_type->GetBinContent(5) + particle_type->GetBinContent(6);
    double e_protons = particle_type->GetBinError(5) + particle_type->GetBinError(6);
    double resonance = particle_type->GetBinContent(7);
    double e_resonance = particle_type->GetBinError(7);
    std::cout << "Pions: " << pions << " ± " << e_pions << '\n';
    std::cout << "Kaons: " << kaons << " ± " << e_kaons << '\n';
    std::cout << "Protons: " << protons << " ± " << e_protons << '\n';
    std::cout << "K* resonance: " << resonance << " ± " << e_resonance << '\n' << '\n';
    std::cout << "Percentage of particle types:" << '\n';
    std::cout << "Pions: " << 100 * pions / particle_type->GetEntries() << "%" << '\n';
    std::cout << "Kaons: " << 100 * kaons / particle_type->GetEntries() << "%" << '\n';
    std::cout << "Protons: " << 100 * protons / particle_type->GetEntries() << "%" << '\n';
    std::cout << "K* resonance: " << 100 * resonance / particle_type->GetEntries() << "%" << '\n' <<'\n';

    TH1F *theta_angle = (TH1F *) histo->Get("theta_angle");
    TH1F *phi_angle = (TH1F *) histo->Get("phi_angle");
    c1->cd(2);
    theta_angle->Fit("pol0", "Q");
    TF1 *theta = theta_angle->GetFunction("pol0");
    std::cout << "Polar angle fit parameter: " << theta->GetParameter(0) << " ± " << theta->GetParError(0) <<'\n';
    theta_angle->Draw();
    c1->cd(3);
    phi_angle->Fit("pol0", "Q");
    TF1 *phi = phi_angle->GetFunction("pol0");
    std::cout << "Azimutal angle fit parameter: " << phi->GetParameter(0) << " ± " << phi->GetParError(0) <<'\n';
    phi_angle->Draw();

    TH1F *impulse = (TH1F *) histo->Get("impulse");
    c1->cd(4);
    impulse->Fit("expo", "Q");
    std::cout << "Impulse mean: " << impulse->GetMean() << " ± " << impulse->GetMeanError() << '\n';
    impulse->Draw();
    
    TH1F *invmass_disc = (TH1F *) histo->Get("invmass_disc");
    TH1F *invmass_conc = (TH1F *) histo->Get("invmass_conc");
    TH1F *invmass_pi_k_disc = (TH1F *) histo->Get("invmass_pi_k_disc");
    TH1F *invmass_pi_k_conc = (TH1F *) histo->Get("invmass_pi_k_conc");
    TH1F *invmass_k = (TH1F *) histo->Get("invmass_k");

    c1->cd(5);
    TH1F *diff_pi_k = new TH1F("diff_pi_k", "Difference between Pi and K of same and opposite sign", 80, 0, 2);
    diff_pi_k->Add(invmass_pi_k_disc, invmass_pi_k_conc, 1, -1);
    diff_pi_k->Fit("gaus", "Q");
    diff_pi_k->Draw();
    TF1 *fit = diff_pi_k->GetFunction("gaus");
    std::cout << "K* mass mean: " << fit->GetParameter(1) << " ± " << fit->GetParError(1) << '\n';
    std::cout << "K* mass width: " << fit->GetParameter(2) << " ± " << fit->GetParError(2) << '\n';
    c1->cd(6);
    TH1F *diff_tot = new TH1F("diff_tot", "Difference between all particles of opposite sign", 80, 0, 2);
    diff_tot->Add(invmass_disc, invmass_conc, 1, -1);
    diff_tot->Fit("gaus", "Q");
    diff_tot->Draw();
    TF1 *fit_tot = diff_tot->GetFunction("gaus");
    std::cout << "K* mass mean: " << fit_tot->GetParameter(1) << " ± " << fit_tot->GetParError(1) << '\n';
    std::cout << "K* mass width: " << fit_tot->GetParameter(2) << " ± " << fit_tot->GetParError(2) << '\n';

}