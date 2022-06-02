
// energy in keV

constexpr double fsc = 1./137.036;
constexpr double pi = 3.141592;
constexpr double k = 8.617333262e-8;

constexpr double amu = 931.5016e3;
constexpr double mO = 15.999*amu;
constexpr int zO = 8;
constexpr double mC = 12*amu;
constexpr int zC = 6;
constexpr double mHe = 4.0015*amu;
constexpr int zHe = 2;
constexpr double u = (mC*mHe)/(mC + mHe);
constexpr double Eg = 2*u*(pi*fsc*zHe*zC)*(pi*fsc*zHe*zC);

using namespace TMath;

double maxwell(double *x, double *p)
{
	return x[0]*Exp(-x[0]/(p[0]*k));
}

double tunnel(double *x, double *p)
{
	return Exp(-Sqrt(Eg/x[0]));
}

double gamow(double *x, double *p)
{
	return maxwell(x, p)*tunnel(x, p);
}

int draw_gamow_peak()
{
	auto f1 = new TF1("Maxwell-Boltzmann", maxwell, 1., 1000., 1);
	auto f2 = new TF1("QM tunneling", tunnel, 1., 1000., 0);
	auto f3 = new TF1("Gamov Peak", gamow, 1., 1000., 1);

	f1->SetParameter(0, 2e8);
	f3->SetParameter(0, 2e8);

	f1->GetYaxis()->SetRangeUser(1e-28, 10);
	f1->SetTitle("Gamov peak at #it{T}=2#times10^{8} K;#it{E}_{C.M} [keV];relative probability");
	f1->GetXaxis()->SetTitleSize(0.055);
	f1->GetXaxis()->SetTitleOffset(0.8);
	f1->GetYaxis()->SetTitleSize(0.055);
	f1->GetYaxis()->SetTitleOffset(1.00);


	f1->SetLineStyle(9);
	f2->SetLineStyle(9);

	f1->SetLineColor(9);
	f2->SetLineColor(6);
	f3->SetLineColor(2);

	f1->SetLineWidth(2);
	f2->SetLineWidth(2);
	f3->SetLineWidth(2);

	TCanvas *c1 = new TCanvas("c1", "c1", 1200, 1000);
	c1->SetLogy();
	c1->SetLeftMargin(1.35);
	c1->SetBottomMargin(1.35);
	f1->Draw();

	f2->Draw("SAME");
	f3->Draw("SAME");
	return 1;
}
