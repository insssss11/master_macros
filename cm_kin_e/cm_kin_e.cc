double cmke(double *x, double *p)
{
	double m1 = p[0], m2 = p[1];
	TLorentzVector v1(0., 0., sqrt(2*x[0]*(x[0] + m1)), m1 + x[0]),
				   v2(0., 0., 0., m2);
	// double M = TMath::Sqrt(m1*m1 + m2*m2 + 2*m1*m2 + 2*m1*x[0]);
	double M = (v1 + v2).Mag();
	double mom = TMath::Sqrt((M*M - (m1 + m2)*(m1 + m2))*(M*M - (m1 - m2)*(m1 - m2)))/(2*M);
	
 	return TMath::Sqrt(m1*m1 + mom*mom) - m1 + TMath::Sqrt(m2*m2 + mom*mom) - m2;
}

// T_cm to T_lab
double GetBeamT(double *x, double *p)
{
    double T = x[0], m1 = p[0], m2 = p[1];
    return (T*T + 2*(m1 + m2)*T)/(2*m2);
}

int cm_kin_e()
{
	gStyle->SetTitleFont(132, "XYZ");
	gStyle->SetTitleFont(132, "T");
	gStyle->SetLineWidth(2);
	constexpr double amu = 931.5016, me = 0.509982;
	constexpr double a1 = 12, a2 = 4.002602;
	constexpr double m1 = amu*a1 - 6*me, m2 = amu*a2 - 2*me;
	TCanvas *c1 = new TCanvas("c1", "c1", 1000, 900);
	c1->SetGrid();
	
	TF1 *f1 = new TF1("C12", GetBeamT, 0, 3., 2);
	f1->SetParameters(m1, m2);
	f1->SetTitle("#it{E}_{cm} vs #it{E}_{beam};#it{E}_{cm} [MeV];#it{E}_{beam} [MeV]");
	f1->Draw();
	TF1 *f2 = new TF1("alpha", GetBeamT, 0, 3., 2);
	f2->SetParameters(m2, m1);
	f2->SetLineColor(kBlue);
	f2->Draw("SAME");
	
	TLatex *lt = new TLatex();
	lt->SetTextFont(132);
	lt->SetTextColor(kRed);
	lt->DrawLatexNDC(0.40, 0.7, "^{12}C beam");
	lt->SetTextColor(kBlue);
	lt->DrawLatexNDC(0.65, 0.25, "#alpha beam");
	return 0;
}
