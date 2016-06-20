#include <iostream>

using namespace std;

//----------------------------------
// Variables
//----------------------------------

TTree *ntp_svxseg;

TH1F *h_prec_x;
TH1F *h_prec_y;
TH1F *h_prec_z;

TH1F *h_ew_x;
TH1F *h_ew_y;
TH1F *h_ew_z;

//----------------------------------
// Functions
//----------------------------------

void ComputeBeamSpot()
{
	TFile *fin = new TFile("Data/423844_data_narrowvtx_ew.root");
	ntp_svxseg = (TTree*) fin->Get("ntp_svxseg");

	//Get precise vertex distribution
	ntp_svxseg->Draw("vtx_prec[0]>>h_prec_x(400,-0.5,0.5)", "", "goff");
	h_prec_x = (TH1F*) gDirectory->FindObject("h_prec_x");

	ntp_svxseg->Draw("vtx_prec[1]>>h_prec_y(400,-0.5,0.5)", "", "goff");
	h_prec_y = (TH1F*) gDirectory->FindObject("h_prec_y");

	ntp_svxseg->Draw("vtx_prec[2]>>h_prec_z(400,-20,20)", "", "goff");
	h_prec_z = (TH1F*) gDirectory->FindObject("h_prec_z");

	ntp_svxseg->Draw("vtx_prec_E[0] - vtx_prec_W[0]>>h_ew_x(400,-0.5,0.5)", "", "goff");
	h_ew_x = (TH1F*) gDirectory->FindObject("h_ew_x");

	ntp_svxseg->Draw("vtx_prec_E[1] - vtx_prec_W[1]>>h_ew_y(400,-0.5,0.5)", "", "goff");
	h_ew_y = (TH1F*) gDirectory->FindObject("h_ew_y");

	ntp_svxseg->Draw("vtx_prec_E[2] - vtx_prec_W[2]>>h_ew_z(400,-0.5,0.5)", "", "goff");
	h_ew_z = (TH1F*) gDirectory->FindObject("h_ew_z");

	//Fit histograms with Gaussians
	double r1, r2;
	double p0, p1, p2;

	// --> Precise X
	p0 = h_prec_x->GetMaximum();
	p1 = h_prec_x->GetBinCenter(h_prec_x->GetMaximumBin());
	p2 = 0.35 * h_prec_x->GetRMS();

	r1 = p1 - p2;
	r2 = p1 + p2;

	TF1 *fSignal = new TF1("f_prec_x", "gaus", r1, r2);
	fSignal->SetParameters(p0, p1, p2);

	h_prec_x->Fit("f_prec_x", "Q0R");

	// --> Precise Y
	p0 = h_prec_y->GetMaximum();
	p1 = h_prec_y->GetBinCenter(h_prec_y->GetMaximumBin());
	p2 = 0.35 * h_prec_y->GetRMS();

	r1 = p1 - p2;
	r2 = p1 + p2;

	fSignal = new TF1("f_prec_y", "gaus", r1, r2);
	fSignal->SetParameters(p0, p1, p2);

	h_prec_y->Fit("f_prec_y", "Q0R");

	// --> Precise E-W X
	p0 = h_ew_x->GetMaximum();
	p1 = h_ew_x->GetBinCenter(h_ew_x->GetMaximumBin());
	p2 = 0.35 * h_ew_x->GetRMS();

	r1 = p1 - p2;
	r2 = p1 + p2;

	fSignal = new TF1("f_ew_x", "gaus", r1, r2);
	fSignal->SetParameters(p0, p1, p2);

	h_ew_x->Fit("f_ew_x", "Q0R");

	// --> Precise E-W Y
	p0 = h_ew_y->GetMaximum();
	p1 = h_ew_y->GetBinCenter(h_ew_y->GetMaximumBin());
	p2 = 0.35 * h_ew_y->GetRMS();

	r1 = p1 - p2;
	r2 = p1 + p2;

	fSignal = new TF1("f_ew_y", "gaus", r1, r2);
	fSignal->SetParameters(p0, p1, p2);

	h_ew_y->Fit("f_ew_y", "Q0R");

	// --> Precise E-W Z
	p0 = h_ew_z->GetMaximum();
	p1 = h_ew_z->GetBinCenter(h_ew_z->GetMaximumBin());
	p2 = 0.35 * h_ew_z->GetRMS();

	r1 = p1 - p2;
	r2 = p1 + p2;

	fSignal = new TF1("f_ew_z", "gaus", r1, r2);
	fSignal->SetParameters(p0, p1, p2);

	h_ew_z->Fit("f_ew_z", "Q0R");

	//Plot things
	gStyle->SetOptStat(0);

	TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
	h_prec_x->SetTitle("");
	h_prec_x->GetXaxis()->SetTitle("PRECISE_{x}");
	h_prec_x->GetXaxis()->SetLabelFont(62);
	h_prec_x->GetXaxis()->SetTitleFont(62);
	h_prec_x->GetYaxis()->SetLabelFont(62);
	h_prec_x->GetYaxis()->SetTitleFont(62);
	h_prec_x->Draw();

	TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
	h_prec_y->SetTitle("");
	h_prec_y->GetXaxis()->SetTitle("PRECISE_{y}");
	h_prec_y->GetXaxis()->SetLabelFont(62);
	h_prec_y->GetXaxis()->SetTitleFont(62);
	h_prec_y->GetYaxis()->SetLabelFont(62);
	h_prec_y->GetYaxis()->SetTitleFont(62);
	h_prec_y->Draw();

	TCanvas *c3 = new TCanvas("c3", "c3", 600, 600);
	h_prec_z->SetTitle("");
	h_prec_z->GetXaxis()->SetTitle("PRECISE_{z}");
	h_prec_z->GetXaxis()->SetLabelFont(62);
	h_prec_z->GetXaxis()->SetTitleFont(62);
	h_prec_z->GetYaxis()->SetLabelFont(62);
	h_prec_z->GetYaxis()->SetTitleFont(62);
	h_prec_z->Draw();

	TCanvas *c4 = new TCanvas("c4", "c4", 600, 600);
	h_ew_x->SetTitle("");
	h_ew_x->GetXaxis()->SetTitle("(PRECISE_E - PRECISE_W)_{x}");
	h_ew_x->GetXaxis()->SetLabelFont(62);
	h_ew_x->GetXaxis()->SetTitleFont(62);
	h_ew_x->GetYaxis()->SetLabelFont(62);
	h_ew_x->GetYaxis()->SetTitleFont(62);
	h_ew_x->Draw();

	TCanvas *c5 = new TCanvas("c5", "c5", 600, 600);
	h_ew_y->SetTitle("");
	h_ew_y->GetXaxis()->SetTitle("(PRECISE_E - PRECISE_W)_{y}");
	h_ew_y->GetXaxis()->SetLabelFont(62);
	h_ew_y->GetXaxis()->SetTitleFont(62);
	h_ew_y->GetYaxis()->SetLabelFont(62);
	h_ew_y->GetYaxis()->SetTitleFont(62);
	h_ew_y->Draw();

	TCanvas *c6 = new TCanvas("c6", "c6", 600, 600);
	h_ew_z->SetTitle("");
	h_ew_z->GetXaxis()->SetTitle("(PRECISE_E - PRECISE_W)_{z}");
	h_ew_z->GetXaxis()->SetLabelFont(62);
	h_ew_z->GetXaxis()->SetTitleFont(62);
	h_ew_z->GetYaxis()->SetLabelFont(62);
	h_ew_z->GetYaxis()->SetTitleFont(62);
	h_ew_z->Draw();
}