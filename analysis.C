#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TROOT.h>
#include <iostream>
#include <TCanvas.h>
#include <cmath>

void analysis(){

    TFile *histo = new TFile("Particles.root");
    TH1F *particle_type = (TH1F*) histo->Get("particle_type");
    particle_type->GetXaxis()->SetTitle("Particle type");
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
    theta_angle->GetXaxis()->SetTitle("Angle");
    theta_angle->Fit("pol0", "Q", 0, M_PI);
    TF1 *theta = theta_angle->GetFunction("pol0");
    std::cout << "Polar angle fit parameter: " << theta->GetParameter(0) << " ± " << theta->GetParError(0) <<'\n' <<'\n';
   
    TH1F *phi_angle = (TH1F *) histo->Get("phi_angle");
    theta_angle->GetXaxis()->SetTitle("Angle");
    phi_angle->Fit("pol0", "Q", 0, 2 * M_PI);
    TF1 *phi = phi_angle->GetFunction("pol0");
    std::cout << "Azimutal angle fit parameter: " << phi->GetParameter(0) << " ± " << phi->GetParError(0) <<'\n' <<'\n';

    TH1F *impulse = (TH1F *) histo->Get("impulse"); 
    impulse->GetXaxis()->SetTitle("Impulse (GeV)");
    impulse->Fit("expo", "Q");
    std::cout << "Impulse mean: " << impulse->GetMean() << " ± " << impulse->GetMeanError() << '\n' <<'\n';

    TH1F *invmass_k = (TH1F *) histo->Get("invmass_k");
    invmass_k->GetXaxis()->SetTitle("Invariant mass (GeV)");
    invmass_k->Fit("gaus", "Q");
    TF1 *fit_k = invmass_k->GetFunction("gaus");
    std::cout << "Invariant mass of particles produced after a decay:" << '\n';
    std::cout << "Mean: " << fit_k->GetParameter(1) << " ± " << fit_k->GetParError(1) << '\n';
    std::cout << "Width: " << fit_k->GetParameter(2) << " ± " << fit_k->GetParError(2) << '\n' << '\n';

    TH1F *invmass_pi_k_disc = (TH1F *) histo->Get("invmass_pi_k_disc");
    TH1F *invmass_pi_k_conc = (TH1F *) histo->Get("invmass_pi_k_conc");
    TH1F *diff_pi_k = new TH1F("diff_pi_k", "Difference between Pi and K of same and opposite sign", 80, 0, 2);
    diff_pi_k->Add(invmass_pi_k_disc, invmass_pi_k_conc, 1, -1);
    diff_pi_k->GetXaxis()->SetTitle("Invariant mass (GeV)");
    diff_pi_k->Fit("gaus", "Q");
    TF1 *fit = diff_pi_k->GetFunction("gaus");
    std::cout << "Invariant mass of pions and kaons:" <<'\n';
    std::cout << "Mean: " << fit->GetParameter(1) << " ± " << fit->GetParError(1) << '\n';
    std::cout << "Width: " << fit->GetParameter(2) << " ± " << fit->GetParError(2) << '\n' << '\n';
    
    TH1F *invmass_disc = (TH1F *) histo->Get("invmass_disc");
    TH1F *invmass_conc = (TH1F *) histo->Get("invmass_conc");
    TH1F *diff_tot = new TH1F("diff_tot", "Difference between all particles", 80, 0, 2);
    diff_tot->Add(invmass_disc, invmass_conc, 1, -1);
    diff_tot->GetXaxis()->SetTitle("Invariant mass (GeV)");
    diff_tot->Fit("gaus", "Q");
    TF1 *fit_tot = diff_tot->GetFunction("gaus");
    std::cout << "Total invariant mass:" << '\n';
    std::cout << "Mean: " << fit_tot->GetParameter(1) << " ± " << fit_tot->GetParError(1) << '\n';
    std::cout << "Width: " << fit_tot->GetParameter(2) << " ± " << fit_tot->GetParError(2) << '\n' << '\n';

    TCanvas *c1 = new TCanvas("particles", "Particles generated, angles, impulse");
    c1->Divide(2,2);
    c1->cd(1);
    particle_type->Draw("E, H, SAME");
    c1->cd(2);
    theta_angle->Draw("E, H, SAME");
    c1->cd(3);
    phi_angle->Draw("E, H, SAME");
    c1->cd(4);
    impulse->Draw("E, H, SAME");

    TCanvas *c2 = new TCanvas("invmass", "Invariant mass");
    c2->Divide(1, 3);
    c2->cd(1);
    invmass_k->Draw("E, H, SAME");
    c2->cd(2);
    diff_pi_k->Draw("E, H, SAME");
    c2->cd(3);
    diff_tot->Draw("E, H, SAME");

    TFile *analysis = new TFile("analysis_result.root", "RECREATE");
    particle_type->Write();
    theta_angle->Write();
    phi_angle->Write();
    impulse->Write();
    invmass_k->Write();
    diff_pi_k->Write();
    diff_tot->Write();
    c1->Write();
    c2->Write();
    c1->Print("particles.pdf");
    c2->Print("invariant_mass.pdf");
    analysis->Close();
}