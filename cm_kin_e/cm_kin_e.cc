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

int cm_kin_e()
{
	gStyle->SetLineWidth(2);
	constexpr double amu = 931.5016;
	constexpr double a1 = 12, a2 = 4.002602;
	TF1 *f1 = new TF1("cmke", cmke, 0, 1e1, 2);
	TCanvas *c1 = new TCanvas("c1", "c1", 1000, 900);
	c1->SetGrid();
	
	f1->SetParameters(amu*a1, amu*a2);
	f1->SetTitle("Total kinetic energy in the lab and C.M frame of ^{12}C + ^{4}H;#it{E}_{beam} [MeV];#it{E}_{C.M} [MeV]");
	f1->Draw();
	return 0;
}
