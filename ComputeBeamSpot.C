//------------------------------------------------------
// Macro to compute the beam center and beam spot size
// using the precise vertex distribution measured in
// two independent arms of the VTX detector.
//
// The measurement is carried out as a function of the
// number of tracks in each arm
//
// JDOK
//------------------------------------------------------

#include <iostream>

using namespace std;

//----------------------------------
// Variables
//----------------------------------

//Number of tracks
const int NUMSEG = 4;
string nseg_cuts_lof[4] = {"nvtxtrk_lof_e==nvtxtrk_lof_e && nvtxtrk_lof_w==nvtxtrk_lof_w", "nvtxtrk_lof_e == 2 && nvtxtrk_lof_w == 2", "nvtxtrk_lof_e == 3 && nvtxtrk_lof_w == 3", "nvtxtrk_lof_e == 4 && nvtxtrk_lof_w == 4"};
string nseg_cuts_prec[4] = {"nvtxtrk_prec_e==nvtxtrk_prec_e && nvtxtrk_prec_w==nvtxtrk_prec_w", "nvtxtrk_prec_e == 2 && nvtxtrk_prec_w == 2", "nvtxtrk_prec_e == 3 && nvtxtrk_prec_w == 3", "nvtxtrk_prec_e == 4 && nvtxtrk_prec_w == 4"};

//Tree read in from file
TTree *ntp_svxseg;

//Distribution of vertices from the iterative algorithm
TH1F *h_prec_x[NUMSEG];
TH1F *h_prec_y[NUMSEG];
TH1F *h_prec_z[NUMSEG];

//Distribution of vertices from the LOF algorithm
TH1F *h_lof_x[NUMSEG];
TH1F *h_lof_y[NUMSEG];
TH1F *h_lof_z[NUMSEG];

//East-west vertex difference distributions from iterative algorithm
TH1F *h_ew_prec_x[NUMSEG];
TH1F *h_ew_prec_y[NUMSEG];
TH1F *h_ew_prec_z[NUMSEG];

//East-west vertex difference distributions from LOF algorithm
TH1F *h_ew_lof_x[NUMSEG];
TH1F *h_ew_lof_y[NUMSEG];
TH1F *h_ew_lof_z[NUMSEG];

//Fit functions for each vertex and vertex difference as a function of the number of tracks
TF1 *f_gauss_fits_vtx_x[NUMSEG];
TF1 *f_gauss_fits_vtx_y[NUMSEG];

TF1 *f_gauss_fits_diff_x[NUMSEG];
TF1 *f_gauss_fits_diff_y[NUMSEG];
TF1 *f_gauss_fits_diff_z[NUMSEG];

//Resolution
float resol[NUMSEG] = {0};

//Beam spot size
float bspt[NUMSEG] = {0};

//----------------------------------
// Functions
//----------------------------------

void plotHistograms()
{
	//Plot things
	gStyle->SetOptStat(0);

	TCanvas *cVertex[NUMSEG];

	for (int i = 0; i < NUMSEG; i++)
	{
		cVertex[i] = new TCanvas(Form("cVertex_%i", i), Form("cVertex_%i", i), 1000, 400);
		cVertex[i]->Divide(3, 1);

		cVertex[i]->cd(1);
		h_prec_x[i]->SetTitle("");
		h_prec_x[i]->GetXaxis()->SetTitle("ITERATIVE_{x} [cm]");
		h_prec_x[i]->GetXaxis()->SetLabelFont(62);
		h_prec_x[i]->GetXaxis()->SetTitleFont(62);
		h_prec_x[i]->GetYaxis()->SetLabelFont(62);
		h_prec_x[i]->GetYaxis()->SetTitleFont(62);
		h_prec_x[i]->Draw();
		f_gauss_fits_vtx_x[i]->Draw("same");

		cVertex[i]->cd(2);
		h_prec_y[i]->SetTitle("");
		h_prec_y[i]->GetXaxis()->SetTitle("ITERATIVE_{y} [cm]");
		h_prec_y[i]->GetXaxis()->SetLabelFont(62);
		h_prec_y[i]->GetXaxis()->SetTitleFont(62);
		h_prec_y[i]->GetYaxis()->SetLabelFont(62);
		h_prec_y[i]->GetYaxis()->SetTitleFont(62);
		h_prec_y[i]->Draw();
		f_gauss_fits_vtx_y[i]->Draw("same");

		cVertex[i]->cd(3);
		h_prec_z[i]->SetTitle("");
		h_prec_z[i]->GetXaxis()->SetTitle("ITERATIVE_{z} [cm]");
		h_prec_z[i]->GetXaxis()->SetLabelFont(62);
		h_prec_z[i]->GetXaxis()->SetTitleFont(62);
		h_prec_z[i]->GetYaxis()->SetLabelFont(62);
		h_prec_z[i]->GetYaxis()->SetTitleFont(62);
		h_prec_z[i]->Draw();
	}

	TCanvas *cVertexDiff[NUMSEG];

	for (int i = 0; i < NUMSEG; i++)
	{
		cVertexDiff[i] = new TCanvas(Form("cVertexDiff_%i", i), Form("cVertexDiff_%i", i), 1000, 400);
		cVertexDiff[i]->Divide(3, 1);

		cVertexDiff[i]->cd(1);
		h_ew_prec_x[i]->SetTitle("");
		h_ew_prec_x[i]->GetXaxis()->SetTitle("(ITER_E - ITER_W)_{x} [cm]");
		h_ew_prec_x[i]->GetXaxis()->SetLabelFont(62);
		h_ew_prec_x[i]->GetXaxis()->SetTitleFont(62);
		h_ew_prec_x[i]->GetYaxis()->SetLabelFont(62);
		h_ew_prec_x[i]->GetYaxis()->SetTitleFont(62);
		h_ew_prec_x[i]->Draw();
		f_gauss_fits_diff_x[i]->Draw("same");

		cVertexDiff[i]->cd(2);
		h_ew_prec_y[i]->SetTitle("");
		h_ew_prec_y[i]->GetXaxis()->SetTitle("(ITER_E - ITER_W)_{y} [cm]");
		h_ew_prec_y[i]->GetXaxis()->SetLabelFont(62);
		h_ew_prec_y[i]->GetXaxis()->SetTitleFont(62);
		h_ew_prec_y[i]->GetYaxis()->SetLabelFont(62);
		h_ew_prec_y[i]->GetYaxis()->SetTitleFont(62);
		h_ew_prec_y[i]->Draw();
		f_gauss_fits_diff_y[i]->Draw("same");

		cVertexDiff[i]->cd(3);
		h_ew_prec_z[i]->SetTitle("");
		h_ew_prec_z[i]->GetXaxis()->SetTitle("(ITER_E - ITER_W)_{z} [cm]");
		h_ew_prec_z[i]->GetXaxis()->SetLabelFont(62);
		h_ew_prec_z[i]->GetXaxis()->SetTitleFont(62);
		h_ew_prec_z[i]->GetYaxis()->SetLabelFont(62);
		h_ew_prec_z[i]->GetYaxis()->SetTitleFont(62);
		h_ew_prec_z[i]->Draw();
		f_gauss_fits_diff_z[i]->Draw("same");
	}
}

void fitHistograms()
{
	//Fit histograms with Gaussians
	double r1, r2;
	double p0, p1, p2;

	for (int i = 0; i < NUMSEG; i++)
	{
		// --> Precise X
		p0 = h_prec_x[i]->GetMaximum();
		p1 = h_prec_x[i]->GetBinCenter(h_prec_x[i]->GetMaximumBin());
		p2 = 0.5 * h_prec_x[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_x[i] = new TF1(Form("f_prec_x_%i", i), "gaus", r1, r2);
		f_gauss_fits_vtx_x[i]->SetParameters(p0, p1, p2);

		h_prec_x[i]->Fit(Form("f_prec_x_%i", i), "Q0R");

		// --> Precise Y
		p0 = h_prec_y[i]->GetMaximum();
		p1 = h_prec_y[i]->GetBinCenter(h_prec_y[i]->GetMaximumBin());
		p2 = 0.5 * h_prec_y[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_y[i] = new TF1(Form("f_prec_y_%i", i), "gaus", r1, r2);
		f_gauss_fits_vtx_y[i]->SetParameters(p0, p1, p2);

		h_prec_y[i]->Fit(Form("f_prec_y_%i", i), "Q0R");

		// --> Precise E-W X
		p0 = h_ew_prec_x[i]->GetMaximum();
		p1 = h_ew_prec_x[i]->GetBinCenter(h_ew_prec_x[i]->GetMaximumBin());
		p2 = 0.6 * h_ew_prec_x[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_diff_x[i] = new TF1(Form("f_ew_prec_x_%i", i), "gaus", r1, r2);
		f_gauss_fits_diff_x[i]->SetParameters(p0, p1, p2);

		h_ew_prec_x[i]->Fit(Form("f_ew_prec_x_%i", i), "Q0R");

		// --> Precise E-W Y
		p0 = h_ew_prec_y[i]->GetMaximum();
		p1 = h_ew_prec_y[i]->GetBinCenter(h_ew_prec_y[i]->GetMaximumBin());
		p2 = 0.6 * h_ew_prec_y[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_diff_y[i] = new TF1(Form("f_ew_prec_y_%i", i), "gaus", r1, r2);
		f_gauss_fits_diff_y[i]->SetParameters(p0, p1, p2);

		h_ew_prec_y[i]->Fit(Form("f_ew_prec_y_%i", i), "Q0R");

		// --> Precise E-W Z
		p0 = h_ew_prec_z[i]->GetMaximum();
		p1 = h_ew_prec_z[i]->GetBinCenter(h_ew_prec_z[i]->GetMaximumBin());
		p2 = 0.6 * h_ew_prec_z[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_diff_z[i] = new TF1(Form("f_ew_prec_z_%i", i), "gaus", r1, r2);
		f_gauss_fits_diff_z[i]->SetParameters(p0, p1, p2);

		h_ew_prec_z[i]->Fit(Form("f_ew_prec_z_%i", i), "Q0R");
	}
}

void readHistograms()
{
	TFile *fin = new TFile("Data/423844_data_narrowvtx_ew.root");
	ntp_svxseg = (TTree*) fin->Get("ntp_svxseg");

	for (int i = 0; i < NUMSEG; i++)
	{
		//Get vertex distribution from iterative algorithm
		ntp_svxseg->Draw(Form("vtx_prec[0]>>h_prec_x_%i(400,-0.5,0.5)", i), nseg_cuts_prec[i].c_str(), "goff");
		h_prec_x[i] = (TH1F*) gDirectory->FindObject(Form("h_prec_x_%i", i));

		ntp_svxseg->Draw(Form("vtx_prec[1]>>h_prec_y_%i(400,-0.5,0.5)", i), nseg_cuts_prec[i].c_str(), "goff");
		h_prec_y[i] = (TH1F*) gDirectory->FindObject(Form("h_prec_y_%i", i));

		ntp_svxseg->Draw(Form("vtx_prec[2]>>h_prec_z_%i(400,-20,20)", i), nseg_cuts_prec[i].c_str(), "goff");
		h_prec_z[i] = (TH1F*) gDirectory->FindObject(Form("h_prec_z_%i", i));

		ntp_svxseg->Draw(Form("vtx_prec_E[0] - vtx_prec_W[0]>>h_ew_prec_x_%i(400,-0.5,0.5)", i), nseg_cuts_prec[i].c_str(), "goff");
		h_ew_prec_x[i] = (TH1F*) gDirectory->FindObject(Form("h_ew_prec_x_%i", i));

		ntp_svxseg->Draw(Form("vtx_prec_E[1] - vtx_prec_W[1]>>h_ew_prec_y_%i(400,-0.5,0.5)", i), nseg_cuts_prec[i].c_str(), "goff");
		h_ew_prec_y[i] = (TH1F*) gDirectory->FindObject(Form("h_ew_prec_y_%i", i));

		ntp_svxseg->Draw(Form("vtx_prec_E[2] - vtx_prec_W[2]>>h_ew_prec_z_%i(400,-0.5,0.5)", i), nseg_cuts_prec[i].c_str(), "goff");
		h_ew_prec_z[i] = (TH1F*) gDirectory->FindObject(Form("h_ew_prec_z_%i", i));

		//Get vertex distributions from the LOF algorithm
		ntp_svxseg->Draw(Form("vtx_lof[0]>>h_lof_x_%i(400,-0.5,0.5)", i), nseg_cuts_lof[i].c_str(), "goff");
		h_lof_x[i] = (TH1F*) gDirectory->FindObject(Form("h_lof_x_%i", i));

		ntp_svxseg->Draw(Form("vtx_lof[1]>>h_lof_y_%i(400,-0.5,0.5)", i), nseg_cuts_lof[i].c_str(), "goff");
		h_lof_y[i] = (TH1F*) gDirectory->FindObject(Form("h_lof_y_%i", i));

		ntp_svxseg->Draw(Form("vtx_lof[2]>>h_lof_z_%i(400,-20,20)", i), nseg_cuts_lof[i].c_str(), "goff");
		h_lof_z[i] = (TH1F*) gDirectory->FindObject(Form("h_lof_z_%i", i));

		ntp_svxseg->Draw(Form("vtx_lof_E[0] - vtx_lof_W[0]>>h_ew_lof_x_%i(400,-0.5,0.5)", i), nseg_cuts_lof[i].c_str(), "goff");
		h_ew_lof_x[i] = (TH1F*) gDirectory->FindObject(Form("h_ew_lof_x_%i", i));

		ntp_svxseg->Draw(Form("vtx_lof_E[1] - vtx_lof_W[1]>>h_ew_lof_y_%i(400,-0.5,0.5)", i), nseg_cuts_lof[i].c_str(), "goff");
		h_ew_lof_y[i] = (TH1F*) gDirectory->FindObject(Form("h_ew_lof_y_%i", i));

		ntp_svxseg->Draw(Form("vtx_lof_E[2] - vtx_lof_W[2]>>h_ew_lof_z_%i(400,-0.5,0.5)", i), nseg_cuts_lof[i].c_str(), "goff");
		h_ew_lof_z[i] = (TH1F*) gDirectory->FindObject(Form("h_ew_lof_z_%i", i));
	}
}

void ComputeBeamSpot()
{
	readHistograms();
	fitHistograms();
	plotHistograms();
}