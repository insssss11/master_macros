using namespace TMath;

constexpr double fsc = 1/137.;
constexpr double amu = 931.5016;
constexpr double me = 0.510998;

// inputs
// file : s-factor data file
// x    : Tcm
// y    : s-factor
void ReadData(const string &file, vector<double> &x, vector<double> &y)
{
	x.clear();y.clear();
	char buffer1[256], buffer2[256];
	ifstream fileIn(file, ios_base::in);
	while(!fileIn.eof())
	{
		fileIn.getline(buffer1, 256, ',');
		fileIn.getline(buffer2, 256);
		if(atof(buffer1) == 0.)
			continue;
		x.push_back(atof(buffer1));
		y.push_back(atof(buffer2));
		// cout << x.back() << " " << y.back() << endl;
	}
}

// T_cm to T_lab
double GetBeamT(double *x, double *p)
{
	double T = x[0], m1 = p[0], m2 = p[1];
	return (T*T + 2*(m1 + m2)*T)/(2*m2);
}

// outputs
// e  : T_lab
// xs : crosssection
// inputs
// x  : T_cm
// y  : s-factor
void GetCrossSection(vector<double> &xs, const vector<double> &x, const vector<double> &y, int z1, int z2, double a1, double a2)
{
	xs.clear();
	double m1 = amu*a1 - me*z1;
	double m2 = amu*a2 - me*z2;
	double mu = m1*m2/(m1 + m2);
	double Eg = 2*mu*fsc*fsc*3.141592*3.141592*z1*z1*z2*z2;
	for(size_t n = 0;n < x.size();++n)
	{
		xs.push_back(
				1e-3*Exp(-Sqrt(Eg/(x.at(n))))
				*y.at(n)/x.at(n)
				);
	}
}

void Extrapolate(vector<double> &xs, vector<double> &x, size_t np, double xmin, double x_extra, double y_extra, int z1, int z2, double a1, double a2)
{
	xs.clear();
	x.clear();
	double m1 = amu*a1 - me*z1;
	double m2 = amu*a2 - me*z2;
	double mu = m1*m2/(m1 + m2);
	double Eg = 2*mu*fsc*fsc*3.141592*3.141592*z1*z1*z2*z2;
	cout << y_extra << endl;
	for(size_t n = 0;n < np;++n)
	{
		x.push_back(xmin + n*(x_extra - xmin)/np);
		xs.push_back(
				1e-3*Exp(-Sqrt(Eg/(x.at(n))))
				*y_extra/x.at(n)
				);
	}
}

int crosssection()
{
	gStyle->SetTitleFont(132, "XYZ");
	gStyle->SetTitleFont(132, "T");
	gStyle->SetTextFont(132);
	
	TCanvas *c1 = new TCanvas("c1", "c1", 900, 900);
	c1->SetTopMargin(0.14);
	
	constexpr int z1 = 6, z2 = 2;
	constexpr double a1 = 12, a2 = 4.002602;

	vector<double> x, y, xs;
	ReadData("s_factors.csv", x, y);
	GetCrossSection(xs, x, y, z1, z2, a1, a2);

	double xmin = 0.3, xmax = 5.0;
	double ymin = 1e-15, ymax = 1e-4;

	auto gr0 = new TGraph(xs.size(), x.data(), xs.data());
	gr0->SetLineWidth(3);
	gr0->SetMarkerStyle(23);
	gr0->GetXaxis()->SetLimits(xmin, xmax);
	gr0->GetHistogram()->SetMinimum(ymin);
	gr0->GetHistogram()->SetMaximum(ymax);
	gr0->GetYaxis()->SetTitleOffset(1.4);;
	gr0->SetTitle("Cross section of ^{12}C(#alpha,#gamma)^{16}O reaction;#it{E}_{cm} [MeV];#sigma(#it{E}) [b]");
	gr0->Draw("AP");
	
	vector<double> x_ext, xs_ext;
	Extrapolate(xs_ext, x_ext, 10, xmin, x.front(), y.front(), z1, z2, a1, a2);
	auto gr1 = new TGraph(xs_ext.size(), x_ext.data(), xs_ext.data());
	gr1->SetLineWidth(2);
	gr1->SetLineStyle(9);
	gr1->Draw("CSAME");
	c1->Update();

	double masses[] = {amu*a1 - me*z1, amu*a2 - me*z2};
	
	TF1 *axisFunc = new TF1("axisFunc", GetBeamT, GetBeamT(&xmin, masses), GetBeamT(&xmax, masses), 2);
	axisFunc->SetParameters(masses);
	
	// axis of beam energy in the lab frame on the top side.
	TGaxis *axis = new TGaxis(xmin, gPad->GetUymax(), xmax, gPad->GetUymax(), "axisFunc", 505, "-");
	axis->SetTitle("#it{E}_{lab} [MeV]");
	// axis->SetTitleOffset(0.7);
	axis->SetTextFont(132);
	axis->SetLabelFont(42);
	axis->SetLabelSize(0.035);
	axis->SetTitleOffset(0.5);
	axis->Draw();
	c1->SetLogy();
	c1->SetGrid();
	return 0;
}
	
