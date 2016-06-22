#include <iostream>

using namespace std;

//----------------------------------
// Variables
//----------------------------------

TTree *ntp_svxseg;

TH1F *h_prec_x;
TH1F *h_prec_x_synth;

TH1F *h_prec_y_synth;

TH1F *h_prec_x_e;
TH1F *h_prec_y_e;

TH1F *h_prec_x_w;
TH1F *h_prec_y_w;

TH1F *h_prec_diff_x;
TH1F *h_prec_diff_y;

string trackCut = "nvtxtrk_prec_w == 3 && nvtxtrk_prec_e == 3";

//----------------------------------
// Functions
//----------------------------------

void PlotVertexDistributions()
{
	TFile *fin = new TFile("Data/423844_data_narrowvtx_ew.root");
	ntp_svxseg = (TTree*) fin->Get("ntp_svxseg");

	//Get precise vertex distribution in x
	ntp_svxseg->Draw("vtx_prec[0]>>h_prec_x(200,-0.4,0.4)", trackCut.c_str(), "goff");
	h_prec_x = (TH1F*) gDirectory->FindObject("h_prec_x");

	ntp_svxseg->Draw("(vtx_prec_E[0] + vtx_prec_W[0])/2.0>>h_prec_x_synth(200,-0.4,0.4)", trackCut.c_str(), "goff");
	h_prec_x_synth = (TH1F*) gDirectory->FindObject("h_prec_x_synth");

	ntp_svxseg->Draw("(vtx_prec_E[1] + vtx_prec_W[1])/2.0>>h_prec_y_synth(200,-0.4,0.4)", trackCut.c_str(), "goff");
	h_prec_y_synth = (TH1F*) gDirectory->FindObject("h_prec_y_synth");

	//Get precise vertex distribution in the east only
	ntp_svxseg->Draw("vtx_prec_E[0]>>h_prec_x_e(200,-0.4,0.4)", trackCut.c_str(), "goff");
	h_prec_x_e = (TH1F*) gDirectory->FindObject("h_prec_x_e");

	ntp_svxseg->Draw("vtx_prec_E[1]>>h_prec_y_e(200,-0.4,0.4)", trackCut.c_str(), "goff");
	h_prec_y_e = (TH1F*) gDirectory->FindObject("h_prec_y_e");

	//Get precise vertex distribution in the west only
	ntp_svxseg->Draw("vtx_prec_W[0]>>h_prec_x_w(200,-0.4,0.4)", trackCut.c_str(), "goff");
	h_prec_x_w = (TH1F*) gDirectory->FindObject("h_prec_x_w");

	ntp_svxseg->Draw("vtx_prec_W[1]>>h_prec_y_w(200,-0.4,0.4)", trackCut.c_str(), "goff");
	h_prec_y_w = (TH1F*) gDirectory->FindObject("h_prec_y_w");

	//Get the vertex difference distributions
	ntp_svxseg->Draw("vtx_prec_E[0]-vtx_prec_W[0]>>h_prec_diff_x(400,-0.5,0.5)", trackCut.c_str(), "goff");
	h_prec_diff_x = (TH1F*) gDirectory->FindObject("h_prec_diff_x");

	ntp_svxseg->Draw("vtx_prec_E[1]-vtx_prec_W[1]>>h_prec_diff_y(400,-0.5,0.5)", trackCut.c_str(), "goff");
	h_prec_diff_y = (TH1F*) gDirectory->FindObject("h_prec_diff_y");

	//Fit histograms with Gaussians
	double r1, r2;
	double p0, p1, p2;

	// --> Precise X
	p0 = h_prec_x->GetMaximum();
	p1 = h_prec_x->GetMean();
	p2 = 0.7 * h_prec_x->GetRMS();

	r1 = p1 - p2;
	r2 = p1 + p2;

	TF1 *fSignal = new TF1("f_prec_x", "gaus", r1, r2);
	fSignal->SetParameters(p0, p1, p2);

	h_prec_x->Fit("f_prec_x", "Q0R");

	// --> Precise X Synthetic
	p0 = h_prec_x_synth->GetMaximum();
	p1 = h_prec_x_synth->GetMean();
	p2 = 0.7 * h_prec_x_synth->GetRMS();

	r1 = p1 - p2;
	r2 = p1 + p2;

	fSignal = new TF1("f_prec_x_synth", "gaus", r1, r2);
	fSignal->SetParameters(p0, p1, p2);

	h_prec_x_synth->Fit("f_prec_x_synth", "Q0R");

	// --> Precise X E
	p0 = h_prec_x_e->GetMaximum();
	p1 = h_prec_x_e->GetMean();
	p2 = 0.7 * h_prec_x_e->GetRMS();

	r1 = p1 - p2;
	r2 = p1 + p2;

	fSignal = new TF1("f_prec_x_e", "gaus", r1, r2);
	fSignal->SetParameters(p0, p1, p2);

	h_prec_x_e->Fit("f_prec_x_e", "Q0R");

	// --> Precise Y Synthetic
	p0 = h_prec_y_synth->GetMaximum();
	p1 = h_prec_y_synth->GetMean();
	p2 = 0.7 * h_prec_y_synth->GetRMS();

	r1 = p1 - p2;
	r2 = p1 + p2;

	fSignal = new TF1("f_prec_y_synth", "gaus", r1, r2);
	fSignal->SetParameters(p0, p1, p2);

	h_prec_y_synth->Fit("f_prec_y_synth", "Q0R");

	// --> Precise Y E
	p0 = h_prec_y_e->GetMaximum();
	p1 = h_prec_y_e->GetMean();
	p2 = 0.7 * h_prec_y_e->GetRMS();

	r1 = p1 - p2;
	r2 = p1 + p2;

	fSignal = new TF1("f_prec_y_e", "gaus", r1, r2);
	fSignal->SetParameters(p0, p1, p2);

	h_prec_y_e->Fit("f_prec_y_e", "Q0R");

	// --> Precise X W
	p0 = h_prec_x_w->GetMaximum();
	p1 = h_prec_x_w->GetMean();
	p2 = 0.7 * h_prec_x_w->GetRMS();

	r1 = p1 - p2;
	r2 = p1 + p2;

	fSignal = new TF1("f_prec_x_w", "gaus", r1, r2);
	fSignal->SetParameters(p0, p1, p2);

	h_prec_x_w->Fit("f_prec_x_w", "Q0R");

	// --> Precise Y W
	p0 = h_prec_y_w->GetMaximum();
	p1 = h_prec_y_w->GetMean();
	p2 = 0.7 * h_prec_y_w->GetRMS();

	r1 = p1 - p2;
	r2 = p1 + p2;

	fSignal = new TF1("f_prec_y_w", "gaus", r1, r2);
	fSignal->SetParameters(p0, p1, p2);

	h_prec_y_w->Fit("f_prec_y_w", "Q0R");

	// --> Precise EW X Diff
	p0 = h_prec_diff_x->GetMaximum();
	p1 = h_prec_diff_x->GetMean();
	p2 = 0.7 * h_prec_diff_x->GetRMS();

	r1 = p1 - p2;
	r2 = p1 + p2;

	fSignal = new TF1("f_prec_diff_x", "gaus", r1, r2);
	fSignal->SetParameters(p0, p1, p2);

	h_prec_diff_x->Fit("f_prec_diff_x", "Q0R");

	// --> Precise EW Y Diff
	p0 = h_prec_diff_y->GetMaximum();
	p1 = h_prec_diff_y->GetMean();
	p2 = 0.7 * h_prec_diff_y->GetRMS();

	r1 = p1 - p2;
	r2 = p1 + p2;

	fSignal = new TF1("f_prec_diff_y", "gaus", r1, r2);
	fSignal->SetParameters(p0, p1, p2);

	h_prec_diff_y->Fit("f_prec_diff_y", "Q0R");

	//Plot things
	gStyle->SetOptStat(0);
	/*

		TCanvas *c2 = new TCanvas("c2", "c2", 600, 600);
		h_prec_y->SetTitle("");
		h_prec_y->GetXaxis()->SetTitle("ITER_{y} [cm]");
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
		*/

	TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
	h_prec_x->SetTitle("");
	h_prec_x->GetXaxis()->SetTitle("ITER_{x} [cm]");
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

	TCanvas *c1_1 = new TCanvas("c1_1", "c1_1", 600, 600);
	h_prec_x_synth->SetTitle("");
	h_prec_x_synth->GetXaxis()->SetTitle("ITER_SYNTH_{x} [cm]");
	h_prec_x_synth->GetXaxis()->SetLabelFont(62);
	h_prec_x_synth->GetXaxis()->SetTitleFont(62);
	h_prec_x_synth->GetYaxis()->SetLabelFont(62);
	h_prec_x_synth->GetYaxis()->SetTitleFont(62);
	h_prec_x_synth->GetXaxis()->SetRangeUser(-0.1, 0.35);
	h_prec_x_synth->Draw();
	TF1 *fit_prec_x_synth = h_prec_x_synth->GetFunction("f_prec_x_synth");
	fit_prec_x_synth->Draw("same");
	TLatex *tl_prec_x_synth_mean = new TLatex(0.15, 0.8, Form("#mu = %g", fit_prec_x_synth->GetParameter(1)));
	tl_prec_x_synth_mean->SetNDC(kTRUE);
	tl_prec_x_synth_mean->Draw("same");
	TLatex *tl_prec_x_synth_width = new TLatex(0.15, 0.75, Form("#sigma = %g", fit_prec_x_synth->GetParameter(2)));
	tl_prec_x_synth_width->SetNDC(kTRUE);
	tl_prec_x_synth_width->Draw("same");

	TCanvas *c1_2 = new TCanvas("c1_2", "c1_2", 600, 600);
	h_prec_y_synth->SetTitle("");
	h_prec_y_synth->GetXaxis()->SetTitle("ITER_SYNTH_{y} [cm]");
	h_prec_y_synth->GetXaxis()->SetLabelFont(62);
	h_prec_y_synth->GetXaxis()->SetTitleFont(62);
	h_prec_y_synth->GetYaxis()->SetLabelFont(62);
	h_prec_y_synth->GetYaxis()->SetTitleFont(62);
	h_prec_y_synth->GetXaxis()->SetRangeUser(-0.1, 0.35);
	h_prec_y_synth->Draw();
	TF1 *fit_prec_y_synth = h_prec_y_synth->GetFunction("f_prec_y_synth");
	fit_prec_y_synth->Draw("same");
	TLatex *tl_prec_y_synth_mean = new TLatex(0.15, 0.8, Form("#mu = %g", fit_prec_y_synth->GetParameter(1)));
	tl_prec_y_synth_mean->SetNDC(kTRUE);
	tl_prec_y_synth_mean->Draw("same");
	TLatex *tl_prec_y_synth_width = new TLatex(0.15, 0.75, Form("#sigma = %g", fit_prec_y_synth->GetParameter(2)));
	tl_prec_y_synth_width->SetNDC(kTRUE);
	tl_prec_y_synth_width->Draw("same");

	TCanvas *c7 = new TCanvas("c7", "c7", 600, 600);
	h_prec_diff_x->SetTitle("");
	h_prec_diff_x->GetXaxis()->SetTitle("ITER-EW_{x} [cm]");
	h_prec_diff_x->GetXaxis()->SetLabelFont(62);
	h_prec_diff_x->GetXaxis()->SetTitleFont(62);
	h_prec_diff_x->GetYaxis()->SetLabelFont(62);
	h_prec_diff_x->GetYaxis()->SetTitleFont(62);
	//h_prec_diff_x->GetXaxis()->SetRangeUser(-0.1, 0.25);
	h_prec_diff_x->Draw();
	TF1 *fit_prec_diff_x = h_prec_diff_x->GetFunction("f_prec_diff_x");
	fit_prec_diff_x->Draw("same");
	TLatex *tl_prec_diff_x_mean = new TLatex(0.15, 0.8, Form("#mu = %g", fit_prec_diff_x->GetParameter(1)));
	tl_prec_diff_x_mean->SetNDC(kTRUE);
	tl_prec_diff_x_mean->Draw("same");
	TLatex *tl_prec_diff_x_width = new TLatex(0.15, 0.75, Form("#sigma = %g", fit_prec_diff_x->GetParameter(2)));
	tl_prec_diff_x_width->SetNDC(kTRUE);
	tl_prec_diff_x_width->Draw("same");

	TCanvas *c8 = new TCanvas("c8", "c8", 600, 600);
	h_prec_diff_y->SetTitle("");
	h_prec_diff_y->GetXaxis()->SetTitle("ITER-EW_{y} [cm]");
	h_prec_diff_y->GetXaxis()->SetLabelFont(62);
	h_prec_diff_y->GetXaxis()->SetTitleFont(62);
	h_prec_diff_y->GetYaxis()->SetLabelFont(62);
	h_prec_diff_y->GetYaxis()->SetTitleFont(62);
	//h_prec_diff_y->GetXaxis()->SetRangeUser(-0.1, 0.25);
	h_prec_diff_y->Draw();
	TF1 *fit_prec_diff_y = h_prec_diff_y->GetFunction("f_prec_diff_y");
	fit_prec_diff_y->Draw("same");
	TLatex *tl_prec_diff_y_mean = new TLatex(0.15, 0.8, Form("#mu = %g", fit_prec_diff_y->GetParameter(1)));
	tl_prec_diff_y_mean->SetNDC(kTRUE);
	tl_prec_diff_y_mean->Draw("same");
	TLatex *tl_prec_diff_y_width = new TLatex(0.15, 0.75, Form("#sigma = %g", fit_prec_diff_y->GetParameter(2)));
	tl_prec_diff_y_width->SetNDC(kTRUE);
	tl_prec_diff_y_width->Draw("same");

	TCanvas *c9 = new TCanvas("c9", "c9", 600, 600);
	h_prec_x_e->SetTitle("");
	h_prec_x_e->GetXaxis()->SetTitle("ITER-E_{x} [cm]");
	h_prec_x_e->GetXaxis()->SetLabelFont(62);
	h_prec_x_e->GetXaxis()->SetTitleFont(62);
	h_prec_x_e->GetYaxis()->SetLabelFont(62);
	h_prec_x_e->GetYaxis()->SetTitleFont(62);
	//h_prec_x_e->GetXaxis()->SetRangeUser(-0.1, 0.25);
	h_prec_x_e->Draw();
	TF1 *fit_prec_x_e = h_prec_x_e->GetFunction("f_prec_x_e");
	fit_prec_x_e->Draw("same");
	TLatex *tl_prec_x_e_mean = new TLatex(0.15, 0.8, Form("#mu = %g", fit_prec_x_e->GetParameter(1)));
	tl_prec_x_e_mean->SetNDC(kTRUE);
	tl_prec_x_e_mean->Draw("same");
	TLatex *tl_prec_x_e_width = new TLatex(0.15, 0.75, Form("#sigma = %g", fit_prec_x_e->GetParameter(2)));
	tl_prec_x_e_width->SetNDC(kTRUE);
	tl_prec_x_e_width->Draw("same");

	TCanvas *c10 = new TCanvas("c10", "c10", 600, 600);
	h_prec_y_e->SetTitle("");
	h_prec_y_e->GetXaxis()->SetTitle("ITER-E_{y} [cm]");
	h_prec_y_e->GetXaxis()->SetLabelFont(62);
	h_prec_y_e->GetXaxis()->SetTitleFont(62);
	h_prec_y_e->GetYaxis()->SetLabelFont(62);
	h_prec_y_e->GetYaxis()->SetTitleFont(62);
	//h_prec_y_e->GetXaxis()->SetRangeUser(-0.1, 0.25);
	h_prec_y_e->Draw();
	TF1 *fit_prec_y_e = h_prec_y_e->GetFunction("f_prec_y_e");
	fit_prec_y_e->Draw("same");
	TLatex *tl_prec_y_e_mean = new TLatex(0.15, 0.8, Form("#mu = %g", fit_prec_y_e->GetParameter(1)));
	tl_prec_y_e_mean->SetNDC(kTRUE);
	tl_prec_y_e_mean->Draw("same");
	TLatex *tl_prec_y_e_width = new TLatex(0.15, 0.75, Form("#sigma = %g", fit_prec_y_e->GetParameter(2)));
	tl_prec_y_e_width->SetNDC(kTRUE);
	tl_prec_y_e_width->Draw("same");

	TCanvas *c11 = new TCanvas("c11", "c11", 600, 600);
	h_prec_x_w->SetTitle("");
	h_prec_x_w->GetXaxis()->SetTitle("ITER-W_{x} [cm]");
	h_prec_x_w->GetXaxis()->SetLabelFont(62);
	h_prec_x_w->GetXaxis()->SetTitleFont(62);
	h_prec_x_w->GetYaxis()->SetLabelFont(62);
	h_prec_x_w->GetYaxis()->SetTitleFont(62);
	//h_prec_x_w->GetXaxis()->SetRangeUser(-0.1, 0.25);
	h_prec_x_w->Draw();
	TF1 *fit_prec_x_w = h_prec_x_w->GetFunction("f_prec_x_w");
	fit_prec_x_w->Draw("same");
	TLatex *tl_prec_x_w_mean = new TLatex(0.15, 0.8, Form("#mu = %g", fit_prec_x_w->GetParameter(1)));
	tl_prec_x_w_mean->SetNDC(kTRUE);
	tl_prec_x_w_mean->Draw("same");
	TLatex *tl_prec_x_w_width = new TLatex(0.15, 0.75, Form("#sigma = %g", fit_prec_x_w->GetParameter(2)));
	tl_prec_x_w_width->SetNDC(kTRUE);
	tl_prec_x_w_width->Draw("same");

	TCanvas *c12 = new TCanvas("c12", "c12", 600, 600);
	h_prec_y_w->SetTitle("");
	h_prec_y_w->GetXaxis()->SetTitle("ITER-W_{y} [cm]");
	h_prec_y_w->GetXaxis()->SetLabelFont(62);
	h_prec_y_w->GetXaxis()->SetTitleFont(62);
	h_prec_y_w->GetYaxis()->SetLabelFont(62);
	h_prec_y_w->GetYaxis()->SetTitleFont(62);
	//h_prec_y_w->GetXaxis()->SetRangeUser(-0.1, 0.25);
	h_prec_y_w->Draw();
	TF1 *fit_prec_y_w = h_prec_y_w->GetFunction("f_prec_y_w");
	fit_prec_y_w->Draw("same");
	TLatex *tl_prec_y_w_mean = new TLatex(0.15, 0.8, Form("#mu = %g", fit_prec_y_w->GetParameter(1)));
	tl_prec_y_w_mean->SetNDC(kTRUE);
	tl_prec_y_w_mean->Draw("same");
	TLatex *tl_prec_y_w_width = new TLatex(0.15, 0.75, Form("#sigma = %g", fit_prec_y_w->GetParameter(2)));
	tl_prec_y_w_width->SetNDC(kTRUE);
	tl_prec_y_w_width->Draw("same");

	//Find beam spot size
	float s_x    = fit_prec_x->GetParameter(2);
	float s_x_synth    = fit_prec_x_synth->GetParameter(2);

	float s_x_e  = fit_prec_x_e->GetParameter(2);
	float s_y_e  = fit_prec_y_e->GetParameter(2);

	float s_x_w  = fit_prec_x_w->GetParameter(2);
	float s_y_w  = fit_prec_y_w->GetParameter(2);

	float s_x_ew = fit_prec_diff_x->GetParameter(2);
	float s_y_ew = fit_prec_diff_y->GetParameter(2);

	float s_b_x  = TMath::Sqrt(s_x_synth * s_x_synth - 0.25 * s_x_ew * s_x_ew);

	cout << "sbx = " << s_b_x*10000 << endl;

}