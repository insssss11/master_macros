using namespace TMath;

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
double GetBeamT(double T, double m1, double m2)
{
	return (T*T + 2*(m1 + m2)*T)/(2*m2);
}

// outputs
// e  : T_lab
// xs : crosssection
// inputs
// x  : T_cm
// y  : s-factor
void GetCrossSection(vector<double> &e, vector<double> &xs, const vector<double> &x, const vector<double> &y, int z1, int z2, double a1, double a2)
{
	e.clear();
	xs.clear();
	constexpr double fsc = 1/137.;
	constexpr double amu = 931.5016;
	double mu = amu*a1*a2/(a1 + a2);
	double Eg = 2*mu*fsc*fsc*3.141592*3.141592*z1*z1*z2*z2;
	for(size_t n = 0;n < x.size();++n)
	{
		xs.push_back(
				1e-3*Exp(-Sqrt(Eg/(x.at(n))))
				*y.at(n)/x.at(n)
				);
		e.push_back(GetBeamT(x.at(n), amu*a1, amu*a2));
	}
}


int crosssection()
{
	gStyle->SetTitleFont(132, "XYZ");
	gStyle->SetTitleFont(132, "T");
	gStyle->SetTextFont(132);
	TCanvas *c1 = new TCanvas("c1", "c1", 900, 900);
	c1->SetLogy();
	constexpr int z1 = 6, z2 = 2;
	constexpr double a1 = 12, a2 = 4.002602;

	vector<double> e, x, y, xs;
	ReadData("s_factors.csv", x, y);

	GetCrossSection(e, xs, x, y, z1, z2, a1, a2);
	auto gr = new TGraph(xs.size(), e.data(), xs.data());
	gr->SetLineWidth(3);
	gr->GetYaxis()->SetTitleOffset(1.2);
	gr->SetTitle("#sigma(#it{E}_{lab}) of ^{12}C(#alpha,#gamma)^{16}O reaction;#it{E}_{Lab} [MeV];#sigma(E) [b]");
	gr->SetMarkerStyle(23);
	gr->Draw("AP");
	return 0;
}
	
