#include <iostream>

using namespace std;

//----------------------------------
// Variables
//----------------------------------

TTree *ntp_svxseg;

TH1F *h_prec_x;
TH1F *h_prec_y;
TH1F *h_prec_z;

TH1F *h_lof_x;
TH1F *h_lof_y;
TH1F *h_lof_z;

//----------------------------------
// Functions
//----------------------------------

void PlotVertexDistributions()
{
	TFile *fin = new TFile("Data/423844_1_1_1_data_lof.root");
	ntp_svxseg = (TTree*) fin->Get("ntp_svxseg");

	//Get precise vertex distribution
	ntp_svxseg->Draw("vtx_prec[0]>>h_prec_x(400,-0.5,0.5)", "", "goff");
	h_prec_x = (TH1F*) gDirectory->FindObject("h_prec_x");

	ntp_svxseg->Draw("vtx_prec[1]>>h_prec_y(400,-0.5,0.5)", "", "goff");
	h_prec_y = (TH1F*) gDirectory->FindObject("h_prec_y");

	ntp_svxseg->Draw("vtx_prec[2]>>h_prec_z(400,-20,20)", "", "goff");
	h_prec_z = (TH1F*) gDirectory->FindObject("h_prec_z");

	//Get the LOF vertex distributions
	ntp_svxseg->Draw("vtx_lof[0]>>h_lof_x(400,-0.5,0.5)", "", "goff");
	h_lof_x = (TH1F*) gDirectory->FindObject("h_lof_x");

	ntp_svxseg->Draw("vtx_lof[1]>>h_lof_y(400,-0.5,0.5)", "", "goff");
	h_lof_y = (TH1F*) gDirectory->FindObject("h_lof_y");

	ntp_svxseg->Draw("vtx_lof[2]>>h_lof_z(400,-20,20)", "", "goff");
	h_lof_z = (TH1F*) gDirectory->FindObject("h_lof_z");

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

	// --> LOF X
	p0 = h_lof_x->GetMaximum();
	p1 = h_lof_x->GetBinCenter(h_lof_x->GetMaximumBin());
	p2 = 0.35 * h_lof_x->GetRMS();

	r1 = p1 - p2;
	r2 = p1 + p2;

	fSignal = new TF1("f_lof_x", "gaus", r1, r2);
	fSignal->SetParameters(p0, p1, p2);

	h_lof_x->Fit("f_lof_x", "Q0R");

	// --> LOF Y
	p0 = h_lof_y->GetMaximum();
	p1 = h_lof_y->GetBinCenter(h_lof_y->GetMaximumBin());
	p2 = 0.35 * h_lof_y->GetRMS();

	r1 = p1 - p2;
	r2 = p1 + p2;

	fSignal = new TF1("f_lof_y", "gaus", r1, r2);
	fSignal->SetParameters(p0, p1, p2);

	h_lof_y->Fit("f_lof_y", "Q0R");

	//Plot things
	gStyle->SetOptStat(0);

	TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
	h_prec_x->SetTitle("");
	h_prec_x->GetXaxis()->SetTitle("PRECISE_{x} [cm]");
	h_prec_x->GetXaxis()->SetLabelFont(62);
	h_prec_x->GetXaxis()->SetTitleFont(62);
	h_prec_x->GetYaxis()->SetLabelFont(62);
	h_prec_x->GetYaxis()->SetTitleFont(62);
	h_prec_x->GetXaxis()->SetRangeUser(-0.1, 0.35);
	h_prec_x->Draw();
	TF1 *fit_prec_x = h_prec_x->GetFunction("f_prec_x");
	fit_prec_x->Draw("same");
	TLatex *tl_prec_x_mean = new TLatex(0.15, 0.8, Form("#mu = %g", fit_prec_x->GetParameter(1)));
	tl_prec_x_mean->SetNDC(kTRUE);
	tl_prec_x_mean->Draw("same");
	TLatex *tl_prec_x_width = new TLatex(0.15, 0.75, Form("#sigma = %g", fit_prec_x->GetParameter(2)));
	tl_prec_x_width->SetNDC(kTRUE);
	tl_prec_x_width->Draw("same");

	TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
	h_prec_y->SetTitle("");
	h_prec_y->GetXaxis()->SetTitle("PRECISE_{y} [cm]");
	h_prec_y->GetXaxis()->SetLabelFont(62);
	h_prec_y->GetXaxis()->SetTitleFont(62);
	h_prec_y->GetYaxis()->SetLabelFont(62);
	h_prec_y->GetYaxis()->SetTitleFont(62);
	h_prec_y->GetXaxis()->SetRangeUser(-0.1, 0.25);
	h_prec_y->Draw();
	TF1 *fit_prec_y = h_prec_y->GetFunction("f_prec_y");
	fit_prec_y->Draw("same");
	TLatex *tl_prec_y_mean = new TLatex(0.15, 0.8, Form("#mu = %g", fit_prec_y->GetParameter(1)));
	tl_prec_y_mean->SetNDC(kTRUE);
	tl_prec_y_mean->Draw("same");
	TLatex *tl_prec_y_width = new TLatex(0.15, 0.75, Form("#sigma = %g", fit_prec_y->GetParameter(2)));
	tl_prec_y_width->SetNDC(kTRUE);
	tl_prec_y_width->Draw("same");

	TCanvas *c3 = new TCanvas("c3", "c3", 600, 600);
	h_prec_z->SetTitle("");
	h_prec_z->GetXaxis()->SetTitle("PRECISE_{z} [cm]");
	h_prec_z->GetXaxis()->SetLabelFont(62);
	h_prec_z->GetXaxis()->SetTitleFont(62);
	h_prec_z->GetYaxis()->SetLabelFont(62);
	h_prec_z->GetYaxis()->SetTitleFont(62);
	h_prec_z->Draw();

	TCanvas *c4 = new TCanvas("c4", "c4", 600, 600);
	h_lof_x->SetTitle("");
	h_lof_x->GetXaxis()->SetTitle("LOF_{x} [cm]");
	h_lof_x->GetXaxis()->SetLabelFont(62);
	h_lof_x->GetXaxis()->SetTitleFont(62);
	h_lof_x->GetYaxis()->SetLabelFont(62);
	h_lof_x->GetYaxis()->SetTitleFont(62);
	h_lof_x->GetXaxis()->SetRangeUser(-0.1, 0.35);
	h_lof_x->Draw();
	TF1 *fit_lof_x = h_lof_x->GetFunction("f_lof_x");
	fit_lof_x->Draw("same");
	TLatex *tl_lof_x_mean = new TLatex(0.15, 0.8, Form("#mu = %g", fit_lof_x->GetParameter(1)));
	tl_lof_x_mean->SetNDC(kTRUE);
	tl_lof_x_mean->Draw("same");
	TLatex *tl_lof_x_width = new TLatex(0.15, 0.75, Form("#sigma = %g", fit_lof_x->GetParameter(2)));
	tl_lof_x_width->SetNDC(kTRUE);
	tl_lof_x_width->Draw("same");

	TCanvas *c5 = new TCanvas("c5", "c5", 600, 600);
	h_lof_y->SetTitle("");
	h_lof_y->GetXaxis()->SetTitle("LOF_{y} [cm]");
	h_lof_y->GetXaxis()->SetLabelFont(62);
	h_lof_y->GetXaxis()->SetTitleFont(62);
	h_lof_y->GetYaxis()->SetLabelFont(62);
	h_lof_y->GetYaxis()->SetTitleFont(62);
	h_lof_y->GetXaxis()->SetRangeUser(-0.1, 0.25);
	h_lof_y->Draw();
	TF1 *fit_lof_y = h_lof_y->GetFunction("f_lof_y");
	fit_lof_y->Draw("same");
	TLatex *tl_lof_y_mean = new TLatex(0.15, 0.8, Form("#mu = %g", fit_lof_y->GetParameter(1)));
	tl_lof_y_mean->SetNDC(kTRUE);
	tl_lof_y_mean->Draw("same");
	TLatex *tl_lof_y_width = new TLatex(0.15, 0.75, Form("#sigma = %g", fit_lof_y->GetParameter(2)));
	tl_lof_y_width->SetNDC(kTRUE);
	tl_lof_y_width->Draw("same");

	TCanvas *c6 = new TCanvas("c6", "c6", 600, 600);
	h_lof_z->SetTitle("");
	h_lof_z->GetXaxis()->SetTitle("LOF_{z} [cm]");
	h_lof_z->GetXaxis()->SetLabelFont(62);
	h_lof_z->GetXaxis()->SetTitleFont(62);
	h_lof_z->GetYaxis()->SetLabelFont(62);
	h_lof_z->GetYaxis()->SetTitleFont(62);
	h_lof_z->Draw();
}