using namespace TMath;

double getbeta(double Tcm, double m1, double m2)
{ 

	double E1 = ((Tcm + 2*m1)*(Tcm + 2*m2) - 2*m1*m2)/(2*m2);
	cout << E1 + m2 << endl;
	double mom = Sqrt(E1*E1 - m1*m1);
	cout << "Tbeam = " << Sqrt(m1*m1 + mom*mom) - m1 << endl;
	return mom/(E1 + m2);

}

TGraph *draw_scat_angle_e(double Tcm = 0.9)
{

	gStyle->SetLineWidth(2);
	constexpr int np = 100;
	constexpr double amu = 931.5016;
	constexpr double me = 0.5109;
	// Carbon, Helium, Oxygen, Gamma
	constexpr double a1 = 12, a2 = 4.002602, a3 = 15.999, a4 = 0;
	constexpr double m1 = a1*amu - 2*me, m2 = a2*amu - 2*me, m3 = a3*amu - 8*me, m4 = a4*amu;

	double M = m1 + m2 + Tcm;
	double mom = TMath::Sqrt((M*M - (m3 + m4)*(m3 + m4))*(M*M - (m3 - m4)*(m3 - m4)))/(2*M);
	double beta = getbeta(Tcm, m1, m2);
	// E3cm + E4cm = m1 + m2 + Tcm 
	cout << Sqrt(m3*m3 + mom*mom) + Sqrt(m4*m4 + mom*mom) << " " << beta << endl;
	double ang[np], T[np];
	for(int i = 0;i < np;++i)
	{		
		double ang_cm = i*3.141592/(np);
		TLorentzVector v3_cm(Sin(ang_cm)*mom, 0., Cos(ang_cm)*mom, Sqrt(m3*m3 + mom*mom));
		TLorentzVector v4_cm(-Sin(ang_cm)*mom, 0., -Cos(ang_cm)*mom, Sqrt(m4*m4 + mom*mom));
		v3_cm.Boost(0., 0., beta);
		v4_cm.Boost(0., 0., beta);
		ang[i] = v3_cm.Theta()*RadToDeg();
		T[i] = m3*(v3_cm.Gamma() - 1);
		// cout << T[i] + m3 + v4_cm.E() << endl;
		// sqrt((p3 + p4)*(p3 + p4)) = Ecm = m1 + m2 + Tcm
		// cout << Sqrt(m3*m3 + m4*m4 + 2*v3_cm*v4_cm) << " " << m1 + m2 + Tcm << endl;
	}
	double m = (T[0] + T[np - 1])/2;
	for(auto &t : T)
	{
		t = 100*(t - m)/m;

	}
	auto gr = new TGraph(np, ang, T);
	return gr;
}
TGraph *draw_scat_angle_b(double Tb = 3.6)
{
	constexpr int np = 100;
	constexpr double amu = 931.5016;
	constexpr double me = 0.510998;
	// Carbon, Helium, Oxygen, Gamma
	constexpr double a1 = 12, a2 = 4.002602, a3 = 15.999, a4 = 0;
	constexpr double m1 = a1*amu, m2 = a2*amu - 2*me, m3 = a3*amu - 8*me, m4 = a4*amu;

	double M = Sqrt(m1*m1 + m2*m2 + 2*(Tb + m1)*m2);
	double mom = TMath::Sqrt((M*M - (m3 + m4)*(m3 + m4))*(M*M - (m3 - m4)*(m3 - m4)))/(2*M);
	double beta = Sqrt(Tb*(2*m1 + Tb))/(Tb + m1 + m2);
	cout << Sqrt(m3*m3 + mom*mom) + Sqrt(m4*m4 + mom*mom) << " " << beta << endl;
	double ang[np + 1], T[np + 1];
	TVector3 vec(0., 0., 1.);
	for(int i = 0;i < np + 1;++i)
	{		
		double ang_cm = i*3.141592/(np);
		TLorentzVector v3_cm(Sin(ang_cm)*mom, 0., Cos(ang_cm)*mom, Sqrt(m3*m3 + mom*mom));
		TLorentzVector v4_cm(-Sin(ang_cm)*mom, 0., -Cos(ang_cm)*mom, Sqrt(m4*m4 + mom*mom));
		v3_cm.Boost(0., 0., beta);
		v4_cm.Boost(0., 0., beta);
		ang[i] = v3_cm.Angle(vec)*RadToDeg();
		T[i] = m3*(v3_cm.Gamma() - 1);
	}
	double m = (T[0] + T[np - 1])/2;
	cout << "------------------------------------------------------" << endl;
	cout << "Medium Energy : " << m << " at " << Tb << endl;
	cout << "------------------------------------------------------" << endl;
	for(auto &t : T)
	{
		t = 100*(t - m)/m;
	}
	auto gr = new TGraph(np +1, ang, T);
	return gr;
}

int scat_angle_e()
{
	gStyle->SetTitleFont(132, "XYZ");
    gStyle->SetTitleFont(132, "T");
    gStyle->SetTextFont(132);
	gStyle->SetLineWidth(2);
	// mg->GetHistogram()->GetXaxis()->SetRangeUser(0., 2.);
	// mg->Add(gr1);
	// mg->Add(gr2);
	// mg->Add(gr3);
	// mg->Draw("AC");
	// gr2->Draw();
	auto gr2 = draw_scat_angle_b(9.6);
	auto c1 = new TCanvas("c2", "c2", 1000, 900);
	auto mg = new TMultiGraph();
	double Ek[] = {4.8, 6.0, 7.2, 8.4, 9.6};
	int lineStyle[] = {1, 10, 9, 2, 4, 3};
	for(int i = 0;i < 6;++i)
	{
		auto g = draw_scat_angle_b(Ek[i]);
		g->SetLineStyle(lineStyle[i]);
		g->SetLineWidth(2);
		mg->Add(g);
	}
	TLatex *tl = new TLatex();
	c1->SetLeftMargin(0.12);
	mg->SetTitle(";Lab scattering angle of ^{16}O recoil [deg];#it{E}_{lab} spread [%]");
	mg->GetXaxis()->SetLimits(0, 1.6);
	// mg->GetYaxis()->SetTitleOffset(1.5);
	// mg->GetHistogram()->GetXaxis()->SetRangeUser(0., 2.);
	// mg->Add(gr1);
	// mg->Add(gr2);
	// mg->Add(gr3);
	mg->Draw("AC");
	tl->SetTextSize(30);
	tl->SetTextFont(13);
	tl->DrawLatex(0.2, 2.0, "#it{E}(^{12}C)=9.6 MeV");
	tl->DrawLatex(0.2, 1.5, "#it{E}(^{12}C)=8.4 MeV");
	tl->DrawLatex(0.2, 1.0, "#it{E}(^{12}C)=7.2 MeV");
	tl->DrawLatex(0.2, 0.5, "#it{E}(^{12}C)=6.0 MeV");
	tl->DrawLatex(0.2, 0.0, "#it{E}(^{12}C)=4.8 MeV");
	tl->SetTextSize(22);
	// gr2->Draw();

	TLine *t = new TLine();
	t->SetLineStyle(2);
	t->SetLineWidth(2);

	t->DrawLine(1.091, -3.85, 1.091, 0.);
	tl->DrawLatex(0.92, -4.02, "7.2 MeV");
	
	t->DrawLine(1.127, -3.6, 1.127, 0.);
	tl->DrawLatex(1.15, -3.74, "5.4 MeV");
	
	t->DrawLine(1.178, -3.2, 1.178, 0.);
	tl->DrawLatex(1.20, -3.26, "4.5 MeV");
	
	t->DrawLine(1.238, -2.7, 1.238, 0.);
	tl->DrawLatex(1.25, -2.78, "3.6 MeV");
	
	t->DrawLine(1.332, -2.2, 1.332, 0.);
	tl->DrawLatex(1.34, -2.29, "1.8 MeV");

	return 0;
}
