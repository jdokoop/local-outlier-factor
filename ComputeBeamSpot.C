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

//Number of E/W track combinations
const int NUMSEG = 8;
//Number of total vertex tracks
const int NUMTRACKS = 7;

string nseg_cuts_lof[NUMSEG] = {"nvtxtrk_lof_e>1 && nvtxtrk_lof_w>1", "nvtxtrk_lof_e == 2 && nvtxtrk_lof_w == 2", "nvtxtrk_lof_e == 2 && nvtxtrk_lof_w == 3", "nvtxtrk_lof_e == 3 && nvtxtrk_lof_w == 2", "nvtxtrk_lof_e == 3 && nvtxtrk_lof_w == 3", "nvtxtrk_lof_e == 3 && nvtxtrk_lof_w == 4", "nvtxtrk_lof_e == 4 && nvtxtrk_lof_w == 3", "nvtxtrk_lof_e == 4 && nvtxtrk_lof_w == 4"};
string nseg_cuts_prec[NUMSEG] = {"nvtxtrk_prec_e>1 && nvtxtrk_prec_w>1", "nvtxtrk_prec_e == 2 && nvtxtrk_prec_w == 2", "nvtxtrk_prec_e == 2 && nvtxtrk_prec_w == 3", "nvtxtrk_prec_e == 3 && nvtxtrk_prec_w == 2", "nvtxtrk_prec_e == 3 && nvtxtrk_prec_w == 3", "nvtxtrk_prec_e == 3 && nvtxtrk_prec_w == 4", "nvtxtrk_prec_e == 4 && nvtxtrk_prec_w == 3", "nvtxtrk_prec_e == 4 && nvtxtrk_prec_w == 4"};
string bin_label[NUMSEG] = {"ANY", "2x2", "2x3", "3x2", "3x3", "3x4", "4x3", "4x4"};

//Tree read in from file
TTree *ntp_svxseg;

//Relative frequency of events with different number of total tracks
float eventTrackFrequency[NUMTRACKS];
TGraphErrors *g_track_frequency;

//Distribution of vertices from the iterative algorithm
TH1F *h_prec_x[NUMSEG];
TH1F *h_prec_y[NUMSEG];
TH1F *h_prec_z[NUMSEG];

TH1F *h_prec_x_E[NUMSEG];
TH1F *h_prec_x_W[NUMSEG];
TH1F *h_prec_y_E[NUMSEG];
TH1F *h_prec_y_W[NUMSEG];

//Synthetic distribution of vertices from the interative algorith
TH1F *h_prec_synth_x[NUMSEG];
TH1F *h_prec_synth_y[NUMSEG];

//Distribution of vertices from the LOF algorithm
TH1F *h_lof_x[NUMSEG];
TH1F *h_lof_y[NUMSEG];
TH1F *h_lof_z[NUMSEG];

TH1F *h_lof_x_E[NUMSEG];
TH1F *h_lof_x_W[NUMSEG];
TH1F *h_lof_y_E[NUMSEG];
TH1F *h_lof_y_W[NUMSEG];

//Synthetic distribution of vertices from the LOF algorith
TH1F *h_lof_synth_x[NUMSEG];
TH1F *h_lof_synth_y[NUMSEG];

//East-west vertex difference distributions from iterative algorithm
TH1F *h_ew_prec_x[NUMSEG];
TH1F *h_ew_prec_y[NUMSEG];
TH1F *h_ew_prec_z[NUMSEG];

//East-west vertex difference distributions from LOF algorithm
TH1F *h_ew_lof_x[NUMSEG];
TH1F *h_ew_lof_y[NUMSEG];
TH1F *h_ew_lof_z[NUMSEG];

//Fit functions for each vertex and vertex difference as a function of the number of tracks
TF1 *f_gauss_fits_vtx_prec_x[NUMSEG];
TF1 *f_gauss_fits_vtx_prec_y[NUMSEG];

TF1 *f_gauss_fits_vtx_prec_total_x[NUMSEG];
TF1 *f_gauss_fits_vtx_prec_total_y[NUMSEG];

TF1 *f_gauss_fits_vtx_prec_x_E[NUMSEG];
TF1 *f_gauss_fits_vtx_prec_x_W[NUMSEG];
TF1 *f_gauss_fits_vtx_prec_y_E[NUMSEG];
TF1 *f_gauss_fits_vtx_prec_y_W[NUMSEG];

TF1 *f_gauss_fits_vtx_lof_x[NUMSEG];
TF1 *f_gauss_fits_vtx_lof_y[NUMSEG];

TF1 *f_gauss_fits_vtx_lof_total_x[NUMSEG];
TF1 *f_gauss_fits_vtx_lof_total_y[NUMSEG];

TF1 *f_gauss_fits_vtx_lof_x_E[NUMSEG];
TF1 *f_gauss_fits_vtx_lof_x_W[NUMSEG];
TF1 *f_gauss_fits_vtx_lof_y_E[NUMSEG];
TF1 *f_gauss_fits_vtx_lof_y_W[NUMSEG];

TF1 *f_gauss_fits_diff_prec_x[NUMSEG];
TF1 *f_gauss_fits_diff_prec_y[NUMSEG];
TF1 *f_gauss_fits_diff_prec_z[NUMSEG];

TF1 *f_gauss_fits_diff_lof_x[NUMSEG];
TF1 *f_gauss_fits_diff_lof_y[NUMSEG];
TF1 *f_gauss_fits_diff_lof_z[NUMSEG];

TF1 *f_gauss_fits_vtx_prec_synth_x[NUMSEG];
TF1 *f_gauss_fits_vtx_prec_synth_y[NUMSEG];

TF1 *f_gauss_fits_vtx_lof_synth_x[NUMSEG];
TF1 *f_gauss_fits_vtx_lof_synth_y[NUMSEG];

//Resolution
float resol_lof[2][NUMSEG] = {0};
float resol_prec[2][NUMSEG] = {0};
float resol_lof_err[2][NUMSEG] = {0};
float resol_prec_err[2][NUMSEG] = {0};
TGraphErrors *g_resol_prec_x;
TGraphErrors *g_resol_prec_y;
TGraphErrors *g_bspt_prec_x;
TGraphErrors *g_bspt_prec_y;
TGraph *g_resol_prec_x_aux;
TGraph *g_resol_prec_y_aux;
TGraph *g_bspt_prec_x_aux;
TGraph *g_bspt_prec_y_aux;
TGraphErrors *g_bspt_prec_x_2;
TGraphErrors *g_bspt_prec_y_2;

//Beam spot size
float bspt_lof[2][NUMSEG] = {0};
float bspt_prec[2][NUMSEG] = {0};
float bspt_lof_err[2][NUMSEG] = {0};
float bspt_prec_err[2][NUMSEG] = {0};
TGraphErrors *g_resol_lof_x;
TGraphErrors *g_resol_lof_y;
TGraphErrors *g_bspt_lof_x;
TGraphErrors *g_bspt_lof_y;
TGraph *g_resol_lof_x_aux;
TGraph *g_resol_lof_y_aux;
TGraph *g_bspt_lof_x_aux;
TGraph *g_bspt_lof_y_aux;
TGraphErrors *g_bspt_lof_x_2;
TGraphErrors *g_bspt_lof_y_2;

//Beam spot size via another independent equation
float bspt_lof_2[2][NUMSEG] = {0};
float bspt_prec_2[2][NUMSEG] = {0};

//Precise vertex distributions as a function of the total number of tracks
TH1F *h_prec_total_x[NUMTRACKS];
TH1F *h_prec_total_y[NUMTRACKS];

TH1F *h_lof_total_x[NUMTRACKS];
TH1F *h_lof_total_y[NUMTRACKS];

//Resolution as a function of the total number of vertex tracks per event
float resol_prec_total_x[NUMTRACKS];
float resol_prec_total_y[NUMTRACKS];

float resol_lof_total_x[NUMTRACKS];
float resol_lof_total_y[NUMTRACKS];

TGraphErrors *g_resol_prec_total_x;
TGraphErrors *g_resol_prec_total_y;

TGraphErrors *g_resol_lof_total_x;
TGraphErrors *g_resol_lof_total_y;

//----------------------------------
// Functions
//----------------------------------

void calculateResolution()
{
	for (int i = 0; i < NUMSEG; i++)
	{
		//Errors on quantities
		float err_synth_prec_x = f_gauss_fits_vtx_prec_synth_x[i]->GetParError(2);
		float err_synth_prec_y = f_gauss_fits_vtx_prec_synth_y[i]->GetParError(2);
		float err_diff_prec_x  = f_gauss_fits_diff_prec_x[i]->GetParError(2);
		float err_diff_prec_y  = f_gauss_fits_diff_prec_y[i]->GetParError(2);

		float err_synth_lof_x = f_gauss_fits_vtx_lof_synth_x[i]->GetParError(2);
		float err_synth_lof_y = f_gauss_fits_vtx_lof_synth_y[i]->GetParError(2);
		float err_diff_lof_x  = f_gauss_fits_diff_lof_x[i]->GetParError(2);
		float err_diff_lof_y  = f_gauss_fits_diff_lof_y[i]->GetParError(2);

		float err_prec_x      = f_gauss_fits_vtx_prec_x[i]->GetParError(2);
		float err_prec_y      = f_gauss_fits_vtx_prec_y[i]->GetParError(2);

		float err_lof_x      = f_gauss_fits_vtx_lof_x[i]->GetParError(2);
		float err_lof_y      = f_gauss_fits_vtx_lof_y[i]->GetParError(2);

		//Calculation
		bspt_prec[0][i] = TMath::Sqrt(f_gauss_fits_vtx_prec_synth_x[i]->GetParameter(2) * f_gauss_fits_vtx_prec_synth_x[i]->GetParameter(2) - 0.25 * f_gauss_fits_diff_prec_x[i]->GetParameter(2) * f_gauss_fits_diff_prec_x[i]->GetParameter(2));
		bspt_prec[1][i] = TMath::Sqrt(f_gauss_fits_vtx_prec_synth_y[i]->GetParameter(2) * f_gauss_fits_vtx_prec_synth_y[i]->GetParameter(2) - 0.25 * f_gauss_fits_diff_prec_y[i]->GetParameter(2) * f_gauss_fits_diff_prec_y[i]->GetParameter(2));

		bspt_prec_2[0][i] = TMath::Sqrt(0.5 * (f_gauss_fits_vtx_prec_x_E[i]->GetParameter(2) * f_gauss_fits_vtx_prec_x_E[i]->GetParameter(2) + f_gauss_fits_vtx_prec_x_W[i]->GetParameter(2) * f_gauss_fits_vtx_prec_x_W[i]->GetParameter(2) - f_gauss_fits_diff_prec_x[i]->GetParameter(2) * f_gauss_fits_diff_prec_x[i]->GetParameter(2)));
		bspt_prec_2[1][i] = TMath::Sqrt(0.5 * (f_gauss_fits_vtx_prec_y_E[i]->GetParameter(2) * f_gauss_fits_vtx_prec_y_E[i]->GetParameter(2) + f_gauss_fits_vtx_prec_y_W[i]->GetParameter(2) * f_gauss_fits_vtx_prec_y_W[i]->GetParameter(2) - f_gauss_fits_diff_prec_y[i]->GetParameter(2) * f_gauss_fits_diff_prec_y[i]->GetParameter(2)));

		resol_prec[0][i] = TMath::Sqrt(f_gauss_fits_vtx_prec_x[i]->GetParameter(2) * f_gauss_fits_vtx_prec_x[i]->GetParameter(2) - bspt_prec[0][i] * bspt_prec[0][i]);
		resol_prec[1][i] = TMath::Sqrt(f_gauss_fits_vtx_prec_y[i]->GetParameter(2) * f_gauss_fits_vtx_prec_y[i]->GetParameter(2) - bspt_prec[1][i] * bspt_prec[1][i]);

		bspt_prec_err[0][i] = TMath::Sqrt(pow(f_gauss_fits_vtx_prec_synth_x[i]->GetParameter(2) * err_synth_prec_x, 2) + pow(2 * f_gauss_fits_diff_prec_x[i]->GetParameter(2) * err_diff_prec_x, 2)) / TMath::Sqrt(2 * (pow(f_gauss_fits_vtx_prec_synth_x[i]->GetParameter(2), 2) + pow(f_gauss_fits_diff_prec_x[i]->GetParameter(2), 2)));
		bspt_prec_err[1][i] = TMath::Sqrt(pow(f_gauss_fits_vtx_prec_synth_y[i]->GetParameter(2) * err_synth_prec_y, 2) + pow(2 * f_gauss_fits_diff_prec_y[i]->GetParameter(2) * err_diff_prec_y, 2)) / TMath::Sqrt(2 * (pow(f_gauss_fits_vtx_prec_synth_y[i]->GetParameter(2), 2) + pow(f_gauss_fits_diff_prec_y[i]->GetParameter(2), 2)));

		resol_prec_err[0][i] = TMath::Sqrt(pow(f_gauss_fits_vtx_prec_x[i]->GetParameter(2) * err_prec_x, 2) + pow(2 * bspt_prec[0][i] * bspt_prec_err[0][i], 2)) / TMath::Sqrt(2 * (pow(f_gauss_fits_vtx_prec_x[i]->GetParameter(2), 2) + pow(bspt_prec[0][i], 2)));
		resol_prec_err[1][i] = TMath::Sqrt(pow(f_gauss_fits_vtx_prec_y[i]->GetParameter(2) * err_prec_y, 2) + pow(2 * bspt_prec[1][i] * bspt_prec_err[1][i], 2)) / TMath::Sqrt(2 * (pow(f_gauss_fits_vtx_prec_y[i]->GetParameter(2), 2) + pow(bspt_prec[1][i], 2)));

		bspt_lof[0][i] = TMath::Sqrt(f_gauss_fits_vtx_lof_synth_x[i]->GetParameter(2) * f_gauss_fits_vtx_lof_synth_x[i]->GetParameter(2) - 0.25 * f_gauss_fits_diff_lof_x[i]->GetParameter(2) * f_gauss_fits_diff_lof_x[i]->GetParameter(2));
		bspt_lof[1][i] = TMath::Sqrt(f_gauss_fits_vtx_lof_synth_y[i]->GetParameter(2) * f_gauss_fits_vtx_lof_synth_y[i]->GetParameter(2) - 0.25 * f_gauss_fits_diff_lof_y[i]->GetParameter(2) * f_gauss_fits_diff_lof_y[i]->GetParameter(2));

		bspt_lof_2[0][i] = TMath::Sqrt(0.5 * (f_gauss_fits_vtx_lof_x_E[i]->GetParameter(2) * f_gauss_fits_vtx_lof_x_E[i]->GetParameter(2) + f_gauss_fits_vtx_lof_x_W[i]->GetParameter(2) * f_gauss_fits_vtx_lof_x_W[i]->GetParameter(2) - f_gauss_fits_diff_lof_x[i]->GetParameter(2) * f_gauss_fits_diff_lof_x[i]->GetParameter(2)));
		bspt_lof_2[1][i] = TMath::Sqrt(0.5 * (f_gauss_fits_vtx_lof_y_E[i]->GetParameter(2) * f_gauss_fits_vtx_lof_y_E[i]->GetParameter(2) + f_gauss_fits_vtx_lof_y_W[i]->GetParameter(2) * f_gauss_fits_vtx_lof_y_W[i]->GetParameter(2) - f_gauss_fits_diff_lof_y[i]->GetParameter(2) * f_gauss_fits_diff_lof_y[i]->GetParameter(2)));

		resol_lof[0][i] = TMath::Sqrt(f_gauss_fits_vtx_lof_x[i]->GetParameter(2) * f_gauss_fits_vtx_lof_x[i]->GetParameter(2) - bspt_lof[0][i] * bspt_lof[0][i]);
		resol_lof[1][i] = TMath::Sqrt(f_gauss_fits_vtx_lof_y[i]->GetParameter(2) * f_gauss_fits_vtx_lof_y[i]->GetParameter(2) - bspt_lof[1][i] * bspt_lof[1][i]);

		bspt_lof_err[0][i] = TMath::Sqrt(pow(f_gauss_fits_vtx_lof_synth_x[i]->GetParameter(2) * err_synth_lof_x, 2) + pow(2 * f_gauss_fits_diff_lof_x[i]->GetParameter(2) * err_diff_lof_x, 2)) / TMath::Sqrt(2 * (pow(f_gauss_fits_vtx_lof_synth_x[i]->GetParameter(2), 2) + pow(f_gauss_fits_diff_lof_x[i]->GetParameter(2), 2)));
		bspt_lof_err[1][i] = TMath::Sqrt(pow(f_gauss_fits_vtx_lof_synth_y[i]->GetParameter(2) * err_synth_lof_y, 2) + pow(2 * f_gauss_fits_diff_lof_y[i]->GetParameter(2) * err_diff_lof_y, 2)) / TMath::Sqrt(2 * (pow(f_gauss_fits_vtx_lof_synth_y[i]->GetParameter(2), 2) + pow(f_gauss_fits_diff_lof_y[i]->GetParameter(2), 2)));

		resol_lof_err[0][i] = TMath::Sqrt(pow(f_gauss_fits_vtx_lof_x[i]->GetParameter(2) * err_lof_x, 2) + pow(2 * bspt_lof[0][i] * bspt_lof_err[0][i], 2)) / TMath::Sqrt(2 * (pow(f_gauss_fits_vtx_lof_x[i]->GetParameter(2), 2) + pow(bspt_lof[0][i], 2)));
		resol_lof_err[1][i] = TMath::Sqrt(pow(f_gauss_fits_vtx_lof_y[i]->GetParameter(2) * err_lof_y, 2) + pow(2 * bspt_lof[1][i] * bspt_lof_err[1][i], 2)) / TMath::Sqrt(2 * (pow(f_gauss_fits_vtx_lof_y[i]->GetParameter(2), 2) + pow(bspt_lof[1][i], 2)));
	}
}

void plotHistograms()
{
	//Plot things
	gStyle->SetOptStat(0);

	TCanvas *cVertex[NUMSEG];
	TCanvas *cVertexLOF[NUMSEG];

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
		f_gauss_fits_vtx_prec_x[i]->Draw("same");

		cVertex[i]->cd(2);
		h_prec_y[i]->SetTitle("");
		h_prec_y[i]->GetXaxis()->SetTitle("ITERATIVE_{y} [cm]");
		h_prec_y[i]->GetXaxis()->SetLabelFont(62);
		h_prec_y[i]->GetXaxis()->SetTitleFont(62);
		h_prec_y[i]->GetYaxis()->SetLabelFont(62);
		h_prec_y[i]->GetYaxis()->SetTitleFont(62);
		h_prec_y[i]->Draw();
		f_gauss_fits_vtx_prec_y[i]->Draw("same");

		cVertex[i]->cd(3);
		h_prec_z[i]->SetTitle("");
		h_prec_z[i]->GetXaxis()->SetTitle("ITERATIVE_{z} [cm]");
		h_prec_z[i]->GetXaxis()->SetLabelFont(62);
		h_prec_z[i]->GetXaxis()->SetTitleFont(62);
		h_prec_z[i]->GetYaxis()->SetLabelFont(62);
		h_prec_z[i]->GetYaxis()->SetTitleFont(62);
		h_prec_z[i]->Draw();

		cVertexLOF[i] = new TCanvas(Form("cVertexLOF_%i", i), Form("cVertexLOF_%i", i), 1000, 400);
		cVertexLOF[i]->Divide(3, 1);

		cVertexLOF[i]->cd(1);
		h_lof_x[i]->SetTitle("");
		h_lof_x[i]->GetXaxis()->SetTitle("LOF_{x} [cm]");
		h_lof_x[i]->GetXaxis()->SetLabelFont(62);
		h_lof_x[i]->GetXaxis()->SetTitleFont(62);
		h_lof_x[i]->GetYaxis()->SetLabelFont(62);
		h_lof_x[i]->GetYaxis()->SetTitleFont(62);
		h_lof_x[i]->Draw();
		f_gauss_fits_vtx_lof_x[i]->Draw("same");

		cVertexLOF[i]->cd(2);
		h_lof_y[i]->SetTitle("");
		h_lof_y[i]->GetXaxis()->SetTitle("LOF_{y} [cm]");
		h_lof_y[i]->GetXaxis()->SetLabelFont(62);
		h_lof_y[i]->GetXaxis()->SetTitleFont(62);
		h_lof_y[i]->GetYaxis()->SetLabelFont(62);
		h_lof_y[i]->GetYaxis()->SetTitleFont(62);
		h_lof_y[i]->Draw();
		f_gauss_fits_vtx_lof_y[i]->Draw("same");

		cVertexLOF[i]->cd(3);
		h_lof_z[i]->SetTitle("");
		h_lof_z[i]->GetXaxis()->SetTitle("LOF_{z} [cm]");
		h_lof_z[i]->GetXaxis()->SetLabelFont(62);
		h_lof_z[i]->GetXaxis()->SetTitleFont(62);
		h_lof_z[i]->GetYaxis()->SetLabelFont(62);
		h_lof_z[i]->GetYaxis()->SetTitleFont(62);
		h_lof_z[i]->Draw();
	}

	TCanvas *cVertexDiff[NUMSEG];
	TCanvas *cVertexDiffLOF[NUMSEG];

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
		f_gauss_fits_diff_prec_x[i]->Draw("same");

		cVertexDiff[i]->cd(2);
		h_ew_prec_y[i]->SetTitle("");
		h_ew_prec_y[i]->GetXaxis()->SetTitle("(ITER_E - ITER_W)_{y} [cm]");
		h_ew_prec_y[i]->GetXaxis()->SetLabelFont(62);
		h_ew_prec_y[i]->GetXaxis()->SetTitleFont(62);
		h_ew_prec_y[i]->GetYaxis()->SetLabelFont(62);
		h_ew_prec_y[i]->GetYaxis()->SetTitleFont(62);
		h_ew_prec_y[i]->Draw();
		f_gauss_fits_diff_prec_y[i]->Draw("same");

		cVertexDiff[i]->cd(3);
		h_ew_prec_z[i]->SetTitle("");
		h_ew_prec_z[i]->GetXaxis()->SetTitle("(ITER_E - ITER_W)_{z} [cm]");
		h_ew_prec_z[i]->GetXaxis()->SetLabelFont(62);
		h_ew_prec_z[i]->GetXaxis()->SetTitleFont(62);
		h_ew_prec_z[i]->GetYaxis()->SetLabelFont(62);
		h_ew_prec_z[i]->GetYaxis()->SetTitleFont(62);
		h_ew_prec_z[i]->Draw();
		f_gauss_fits_diff_prec_z[i]->Draw("same");

		cVertexDiffLOF[i] = new TCanvas(Form("cVertexDiffLOF_%i", i), Form("cVertexDiffLOF_%i", i), 1000, 400);
		cVertexDiffLOF[i]->Divide(3, 1);

		cVertexDiffLOF[i]->cd(1);
		h_ew_lof_x[i]->SetTitle("");
		h_ew_lof_x[i]->GetXaxis()->SetTitle("(LOF_E - LOF_W)_{x} [cm]");
		h_ew_lof_x[i]->GetXaxis()->SetLabelFont(62);
		h_ew_lof_x[i]->GetXaxis()->SetTitleFont(62);
		h_ew_lof_x[i]->GetYaxis()->SetLabelFont(62);
		h_ew_lof_x[i]->GetYaxis()->SetTitleFont(62);
		h_ew_lof_x[i]->Draw();
		f_gauss_fits_diff_lof_x[i]->Draw("same");

		cVertexDiffLOF[i]->cd(2);
		h_ew_lof_y[i]->SetTitle("");
		h_ew_lof_y[i]->GetXaxis()->SetTitle("(LOF_E - LOF_W)_{y} [cm]");
		h_ew_lof_y[i]->GetXaxis()->SetLabelFont(62);
		h_ew_lof_y[i]->GetXaxis()->SetTitleFont(62);
		h_ew_lof_y[i]->GetYaxis()->SetLabelFont(62);
		h_ew_lof_y[i]->GetYaxis()->SetTitleFont(62);
		h_ew_lof_y[i]->Draw();
		f_gauss_fits_diff_lof_y[i]->Draw("same");

		cVertexDiffLOF[i]->cd(3);
		h_ew_lof_z[i]->SetTitle("");
		h_ew_lof_z[i]->GetXaxis()->SetTitle("(LOF_E - LOF_W)_{z} [cm]");
		h_ew_lof_z[i]->GetXaxis()->SetLabelFont(62);
		h_ew_lof_z[i]->GetXaxis()->SetTitleFont(62);
		h_ew_lof_z[i]->GetYaxis()->SetLabelFont(62);
		h_ew_lof_z[i]->GetYaxis()->SetTitleFont(62);
		h_ew_lof_z[i]->Draw();
		f_gauss_fits_diff_lof_z[i]->Draw("same");
	}
}

void fitHistograms()
{
	//Fit histograms with Gaussians
	//Do the fit in two iterations: First, seed fit with curve-derived parameters
	//Then, take parameters from first fit to fit again over a narrower range
	double r1, r2;
	double p0, p1, p2;

	double vertexRMS = 1.5;
	double vertexDiffRMS = 1.5;

	for (int i = 0; i < NUMTRACKS; i++)
	{
		// --> Precise Total X
		p0 = h_prec_total_x[i]->GetMaximum();
		p1 = h_prec_total_x[i]->GetMean();
		p2 = vertexRMS * h_prec_total_x[i]->GetRMS();
		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_prec_total_x[i] = new TF1(Form("f_prec_total_x_%i", i), "gaus", r1, r2);
		f_gauss_fits_vtx_prec_total_x[i]->SetParameters(p0, p1, p2);

		h_prec_total_x[i]->Fit(Form("f_prec_total_x_%i", i), "Q0R");

		p0 = f_gauss_fits_vtx_prec_total_x[i]->GetParameter(0);
		p1 = f_gauss_fits_vtx_prec_total_x[i]->GetParameter(1);
		p2 = vertexRMS * f_gauss_fits_vtx_prec_total_x[i]->GetParameter(2);

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_prec_total_x[i]->SetParameters(p0, p1, p2);
		f_gauss_fits_vtx_prec_total_x[i]->SetRange(r1, r2);
		h_prec_total_x[i]->Fit(Form("f_prec_total_x_%i", i), "Q0R");

		// --> Precise Total Y
		p0 = h_prec_total_y[i]->GetMaximum();
		p1 = h_prec_total_y[i]->GetMean();
		p2 = vertexRMS * h_prec_total_y[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_prec_total_y[i] = new TF1(Form("f_prec_total_y_%i", i), "gaus", r1, r2);
		f_gauss_fits_vtx_prec_total_y[i]->SetParameters(p0, p1, p2);

		h_prec_total_y[i]->Fit(Form("f_prec_total_y_%i", i), "Q0R");

		p0 = f_gauss_fits_vtx_prec_total_y[i]->GetParameter(0);
		p1 = f_gauss_fits_vtx_prec_total_y[i]->GetParameter(1);
		p2 = vertexRMS * f_gauss_fits_vtx_prec_total_y[i]->GetParameter(2);

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_prec_total_y[i]->SetParameters(p0, p1, p2);
		f_gauss_fits_vtx_prec_total_y[i]->SetRange(r1, r2);
		h_prec_total_y[i]->Fit(Form("f_prec_total_y_%i", i), "Q0R");

		// --> LOF Total X
		p0 = h_lof_total_x[i]->GetMaximum();
		p1 = h_lof_total_x[i]->GetMean();
		p2 = vertexRMS * h_lof_total_x[i]->GetRMS();
		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_lof_total_x[i] = new TF1(Form("f_lof_total_x_%i", i), "gaus", r1, r2);
		f_gauss_fits_vtx_lof_total_x[i]->SetParameters(p0, p1, p2);

		h_lof_total_x[i]->Fit(Form("f_lof_total_x_%i", i), "Q0R");

		p0 = f_gauss_fits_vtx_lof_total_x[i]->GetParameter(0);
		p1 = f_gauss_fits_vtx_lof_total_x[i]->GetParameter(1);
		p2 = vertexRMS * f_gauss_fits_vtx_lof_total_x[i]->GetParameter(2);

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_lof_total_x[i]->SetParameters(p0, p1, p2);
		f_gauss_fits_vtx_lof_total_x[i]->SetRange(r1, r2);
		h_lof_total_x[i]->Fit(Form("f_lof_total_x_%i", i), "Q0R");

		// --> LOF Total Y
		p0 = h_lof_total_y[i]->GetMaximum();
		p1 = h_lof_total_y[i]->GetMean();
		p2 = vertexRMS * h_lof_total_y[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_lof_total_y[i] = new TF1(Form("f_lof_total_y_%i", i), "gaus", r1, r2);
		f_gauss_fits_vtx_lof_total_y[i]->SetParameters(p0, p1, p2);

		h_lof_total_y[i]->Fit(Form("f_lof_total_y_%i", i), "Q0R");

		p0 = f_gauss_fits_vtx_lof_total_y[i]->GetParameter(0);
		p1 = f_gauss_fits_vtx_lof_total_y[i]->GetParameter(1);
		p2 = vertexRMS * f_gauss_fits_vtx_lof_total_y[i]->GetParameter(2);

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_lof_total_y[i]->SetParameters(p0, p1, p2);
		f_gauss_fits_vtx_lof_total_y[i]->SetRange(r1, r2);
		h_lof_total_y[i]->Fit(Form("f_lof_total_y_%i", i), "Q0R");
	}

	for (int i = 0; i < NUMSEG; i++)
	{
		if (NUMSEG == 2 || NUMSEG == 3)
		{
			vertexRMS = 1.5;
			vertexDiffRMS = 1.5;
		}

		// --> Precise X
		p0 = h_prec_x[i]->GetMaximum();
		p1 = h_prec_x[i]->GetMean();
		p2 = vertexRMS * h_prec_x[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_prec_x[i] = new TF1(Form("f_prec_x_%i", i), "gaus", r1, r2);
		f_gauss_fits_vtx_prec_x[i]->SetParameters(p0, p1, p2);

		h_prec_x[i]->Fit(Form("f_prec_x_%i", i), "Q0R");

		p0 = f_gauss_fits_vtx_prec_x[i]->GetParameter(0);
		p1 = f_gauss_fits_vtx_prec_x[i]->GetParameter(1);
		p2 = vertexRMS * f_gauss_fits_vtx_prec_x[i]->GetParameter(2);

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_prec_x[i]->SetParameters(p0, p1, p2);
		f_gauss_fits_vtx_prec_x[i]->SetRange(r1, r2);
		h_prec_x[i]->Fit(Form("f_prec_x_%i", i), "Q0R");

		// --> Precise X E
		p0 = h_prec_x_E[i]->GetMaximum();
		p1 = h_prec_x_E[i]->GetMean();
		p2 = vertexRMS * h_prec_x_E[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_prec_x_E[i] = new TF1(Form("f_prec_x_E_%i", i), "gaus", r1, r2);
		f_gauss_fits_vtx_prec_x_E[i]->SetParameters(p0, p1, p2);

		h_prec_x_E[i]->Fit(Form("f_prec_x_E_%i", i), "Q0R");

		p0 = f_gauss_fits_vtx_prec_x_E[i]->GetParameter(0);
		p1 = f_gauss_fits_vtx_prec_x_E[i]->GetParameter(1);
		p2 = vertexRMS * f_gauss_fits_vtx_prec_x_E[i]->GetParameter(2);

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_prec_x_E[i]->SetParameters(p0, p1, p2);
		f_gauss_fits_vtx_prec_x_E[i]->SetRange(r1, r2);
		h_prec_x_E[i]->Fit(Form("f_prec_x_E_%i", i), "Q0R");

		// --> Precise X W
		p0 = h_prec_x_W[i]->GetMaximum();
		p1 = h_prec_x_W[i]->GetMean();
		p2 = vertexRMS * h_prec_x_W[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_prec_x_W[i] = new TF1(Form("f_prec_x_W_%i", i), "gaus", r1, r2);
		f_gauss_fits_vtx_prec_x_W[i]->SetParameters(p0, p1, p2);

		h_prec_x_W[i]->Fit(Form("f_prec_x_W_%i", i), "Q0R");

		p0 = f_gauss_fits_vtx_prec_x_W[i]->GetParameter(0);
		p1 = f_gauss_fits_vtx_prec_x_W[i]->GetParameter(1);
		p2 = vertexRMS * f_gauss_fits_vtx_prec_x_W[i]->GetParameter(2);

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_prec_x_W[i]->SetParameters(p0, p1, p2);
		f_gauss_fits_vtx_prec_x_W[i]->SetRange(r1, r2);
		h_prec_x_W[i]->Fit(Form("f_prec_x_W_%i", i), "Q0R");

		// --> Synthetic Precise X
		p0 = h_prec_synth_x[i]->GetMaximum();
		p1 = h_prec_synth_x[i]->GetMean();
		p2 = vertexRMS * h_prec_synth_x[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_prec_synth_x[i] = new TF1(Form("f_prec_synth_x_%i", i), "gaus", r1, r2);
		f_gauss_fits_vtx_prec_synth_x[i]->SetParameters(p0, p1, p2);

		h_prec_synth_x[i]->Fit(Form("f_prec_synth_x_%i", i), "Q0R");

		p0 = f_gauss_fits_vtx_prec_synth_x[i]->GetParameter(0);
		p1 = f_gauss_fits_vtx_prec_synth_x[i]->GetParameter(1);
		p2 = vertexRMS * f_gauss_fits_vtx_prec_synth_x[i]->GetParameter(2);

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_prec_synth_x[i]->SetParameters(p0, p1, p2);
		f_gauss_fits_vtx_prec_synth_x[i]->SetRange(r1, r2);
		h_prec_synth_x[i]->Fit(Form("f_prec_synth_x_%i", i), "Q0R");

		// --> Precise Y
		p0 = h_prec_y[i]->GetMaximum();
		p1 = h_prec_y[i]->GetMean();
		p2 = vertexRMS * h_prec_y[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_prec_y[i] = new TF1(Form("f_prec_y_%i", i), "gaus", r1, r2);
		f_gauss_fits_vtx_prec_y[i]->SetParameters(p0, p1, p2);

		h_prec_y[i]->Fit(Form("f_prec_y_%i", i), "Q0R");

		p0 = f_gauss_fits_vtx_prec_y[i]->GetParameter(0);
		p1 = f_gauss_fits_vtx_prec_y[i]->GetParameter(1);
		p2 = vertexRMS * f_gauss_fits_vtx_prec_y[i]->GetParameter(2);

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_prec_y[i]->SetParameters(p0, p1, p2);
		f_gauss_fits_vtx_prec_y[i]->SetRange(r1, r2);
		h_prec_y[i]->Fit(Form("f_prec_y_%i", i), "Q0R");

		// --> Precise Y E
		p0 = h_prec_y_E[i]->GetMaximum();
		p1 = h_prec_y_E[i]->GetMean();
		p2 = vertexRMS * h_prec_y_E[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_prec_y_E[i] = new TF1(Form("f_prec_y_E_%i", i), "gaus", r1, r2);
		f_gauss_fits_vtx_prec_y_E[i]->SetParameters(p0, p1, p2);

		h_prec_y_E[i]->Fit(Form("f_prec_y_E_%i", i), "Q0R");

		p0 = f_gauss_fits_vtx_prec_y_E[i]->GetParameter(0);
		p1 = f_gauss_fits_vtx_prec_y_E[i]->GetParameter(1);
		p2 = vertexRMS * f_gauss_fits_vtx_prec_y_E[i]->GetParameter(2);

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_prec_y_E[i]->SetParameters(p0, p1, p2);
		f_gauss_fits_vtx_prec_y_E[i]->SetRange(r1, r2);
		h_prec_y_E[i]->Fit(Form("f_prec_y_E_%i", i), "Q0R");

		// --> Precise Y W
		p0 = h_prec_y_W[i]->GetMaximum();
		p1 = h_prec_y_W[i]->GetMean();
		p2 = vertexRMS * h_prec_y_W[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_prec_y_W[i] = new TF1(Form("f_prec_y_W_%i", i), "gaus", r1, r2);
		f_gauss_fits_vtx_prec_y_W[i]->SetParameters(p0, p1, p2);

		h_prec_y_W[i]->Fit(Form("f_prec_y_W_%i", i), "Q0R");

		p0 = f_gauss_fits_vtx_prec_y_W[i]->GetParameter(0);
		p1 = f_gauss_fits_vtx_prec_y_W[i]->GetParameter(1);
		p2 = vertexRMS * f_gauss_fits_vtx_prec_y_W[i]->GetParameter(2);

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_prec_y_W[i]->SetParameters(p0, p1, p2);
		f_gauss_fits_vtx_prec_y_W[i]->SetRange(r1, r2);
		h_prec_y_W[i]->Fit(Form("f_prec_y_W_%i", i), "Q0R");

		// --> Synthetic Precise Y
		p0 = h_prec_synth_y[i]->GetMaximum();
		p1 = h_prec_synth_y[i]->GetMean();
		p2 = vertexRMS * h_prec_synth_y[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_prec_synth_y[i] = new TF1(Form("f_prec_synth_y_%i", i), "gaus", r1, r2);
		f_gauss_fits_vtx_prec_synth_y[i]->SetParameters(p0, p1, p2);

		h_prec_synth_y[i]->Fit(Form("f_prec_synth_y_%i", i), "Q0R");

		p0 = f_gauss_fits_vtx_prec_synth_y[i]->GetParameter(0);
		p1 = f_gauss_fits_vtx_prec_synth_y[i]->GetParameter(1);
		p2 = vertexRMS * f_gauss_fits_vtx_prec_synth_y[i]->GetParameter(2);

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_prec_synth_y[i]->SetParameters(p0, p1, p2);
		f_gauss_fits_vtx_prec_synth_y[i]->SetRange(r1, r2);
		h_prec_synth_y[i]->Fit(Form("f_prec_synth_y_%i", i), "Q0R");

		// --> Precise E-W X
		p0 = h_ew_prec_x[i]->GetMaximum();
		p1 = h_ew_prec_x[i]->GetMean();
		p2 = vertexDiffRMS * h_ew_prec_x[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_diff_prec_x[i] = new TF1(Form("f_ew_prec_x_%i", i), "gaus", r1, r2);
		f_gauss_fits_diff_prec_x[i]->SetParameters(p0, p1, p2);

		h_ew_prec_x[i]->Fit(Form("f_ew_prec_x_%i", i), "Q0R");

		p0 = f_gauss_fits_diff_prec_x[i]->GetParameter(0);
		p1 = f_gauss_fits_diff_prec_x[i]->GetParameter(1);
		p2 = vertexRMS * f_gauss_fits_diff_prec_x[i]->GetParameter(2);

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_diff_prec_x[i]->SetParameters(p0, p1, p2);
		f_gauss_fits_diff_prec_x[i]->SetRange(r1, r2);
		h_ew_prec_x[i]->Fit(Form("f_ew_prec_x_%i", i), "Q0R");

		// --> Precise E-W Y
		p0 = h_ew_prec_y[i]->GetMaximum();
		p1 = h_ew_prec_y[i]->GetMean();
		p2 = vertexDiffRMS * h_ew_prec_y[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_diff_prec_y[i] = new TF1(Form("f_ew_prec_y_%i", i), "gaus", r1, r2);
		f_gauss_fits_diff_prec_y[i]->SetParameters(p0, p1, p2);

		h_ew_prec_y[i]->Fit(Form("f_ew_prec_y_%i", i), "Q0R");

		p0 = f_gauss_fits_diff_prec_y[i]->GetParameter(0);
		p1 = f_gauss_fits_diff_prec_y[i]->GetParameter(1);
		p2 = vertexRMS * f_gauss_fits_diff_prec_y[i]->GetParameter(2);

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_diff_prec_y[i]->SetParameters(p0, p1, p2);
		f_gauss_fits_diff_prec_y[i]->SetRange(r1, r2);
		h_ew_prec_y[i]->Fit(Form("f_ew_prec_y_%i", i), "Q0R");

		// --> Precise E-W Z
		p0 = h_ew_prec_z[i]->GetMaximum();
		p1 = h_ew_prec_z[i]->GetMean();
		p2 = vertexDiffRMS * h_ew_prec_z[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_diff_prec_z[i] = new TF1(Form("f_ew_prec_z_%i", i), "gaus", r1, r2);
		f_gauss_fits_diff_prec_z[i]->SetParameters(p0, p1, p2);

		h_ew_prec_z[i]->Fit(Form("f_ew_prec_z_%i", i), "Q0R");

		p0 = f_gauss_fits_diff_prec_z[i]->GetParameter(0);
		p1 = f_gauss_fits_diff_prec_z[i]->GetParameter(1);
		p2 = vertexRMS * f_gauss_fits_diff_prec_z[i]->GetParameter(2);

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_diff_prec_z[i]->SetParameters(p0, p1, p2);
		f_gauss_fits_diff_prec_z[i]->SetRange(r1, r2);
		h_ew_prec_z[i]->Fit(Form("f_ew_prec_z_%i", i), "Q0R");

		// --> LOF X
		p0 = h_lof_x[i]->GetMaximum();
		p1 = h_lof_x[i]->GetMean();
		p2 = vertexRMS * h_lof_x[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_lof_x[i] = new TF1(Form("f_lof_x_%i", i), "gaus", r1, r2);
		f_gauss_fits_vtx_lof_x[i]->SetParameters(p0, p1, p2);

		h_lof_x[i]->Fit(Form("f_lof_x_%i", i), "Q0R");

		p0 = f_gauss_fits_vtx_lof_x[i]->GetParameter(0);
		p1 = f_gauss_fits_vtx_lof_x[i]->GetParameter(1);
		p2 = vertexRMS * f_gauss_fits_vtx_lof_x[i]->GetParameter(2);

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_lof_x[i]->SetParameters(p0, p1, p2);
		f_gauss_fits_vtx_lof_x[i]->SetRange(r1, r2);
		h_lof_x[i]->Fit(Form("f_lof_x_%i", i), "Q0R");

		// --> Synthetic LOF X
		p0 = h_lof_synth_x[i]->GetMaximum();
		p1 = h_lof_synth_x[i]->GetMean();
		p2 = vertexRMS * h_lof_synth_x[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_lof_synth_x[i] = new TF1(Form("f_lof_synth_x_%i", i), "gaus", r1, r2);
		f_gauss_fits_vtx_lof_synth_x[i]->SetParameters(p0, p1, p2);

		h_lof_synth_x[i]->Fit(Form("f_lof_synth_x_%i", i), "Q0R");

		p0 = f_gauss_fits_vtx_lof_synth_x[i]->GetParameter(0);
		p1 = f_gauss_fits_vtx_lof_synth_x[i]->GetParameter(1);
		p2 = vertexRMS * f_gauss_fits_vtx_lof_synth_x[i]->GetParameter(2);

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_lof_synth_x[i]->SetParameters(p0, p1, p2);
		f_gauss_fits_vtx_lof_synth_x[i]->SetRange(r1, r2);
		h_lof_synth_x[i]->Fit(Form("f_lof_synth_x_%i", i), "Q0R");

		// --> LOF X E
		p0 = h_lof_x_E[i]->GetMaximum();
		p1 = h_lof_x_E[i]->GetMean();
		p2 = vertexRMS * h_lof_x_E[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_lof_x_E[i] = new TF1(Form("f_lof_x_E_%i", i), "gaus", r1, r2);
		f_gauss_fits_vtx_lof_x_E[i]->SetParameters(p0, p1, p2);

		h_lof_x_E[i]->Fit(Form("f_lof_x_E_%i", i), "Q0R");

		p0 = f_gauss_fits_vtx_lof_x_E[i]->GetParameter(0);
		p1 = f_gauss_fits_vtx_lof_x_E[i]->GetParameter(1);
		p2 = vertexRMS * f_gauss_fits_vtx_lof_x_E[i]->GetParameter(2);

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_lof_x_E[i]->SetParameters(p0, p1, p2);
		f_gauss_fits_vtx_lof_x_E[i]->SetRange(r1, r2);
		h_lof_x_E[i]->Fit(Form("f_lof_x_E_%i", i), "Q0R");

		// --> LOF X W
		p0 = h_lof_x_W[i]->GetMaximum();
		p1 = h_lof_x_W[i]->GetMean();
		p2 = vertexRMS * h_lof_x_W[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_lof_x_W[i] = new TF1(Form("f_lof_x_W_%i", i), "gaus", r1, r2);
		f_gauss_fits_vtx_lof_x_W[i]->SetParameters(p0, p1, p2);

		h_lof_x_W[i]->Fit(Form("f_lof_x_W_%i", i), "Q0R");

		p0 = f_gauss_fits_vtx_lof_x_W[i]->GetParameter(0);
		p1 = f_gauss_fits_vtx_lof_x_W[i]->GetParameter(1);
		p2 = vertexRMS * f_gauss_fits_vtx_lof_x_W[i]->GetParameter(2);

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_lof_x_W[i]->SetParameters(p0, p1, p2);
		f_gauss_fits_vtx_lof_x_W[i]->SetRange(r1, r2);
		h_lof_x_W[i]->Fit(Form("f_lof_x_W_%i", i), "Q0R");

		// --> LOF Y
		p0 = h_lof_y[i]->GetMaximum();
		p1 = h_lof_y[i]->GetMean();
		p2 = vertexRMS * h_lof_y[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_lof_y[i] = new TF1(Form("f_lof_y_%i", i), "gaus", r1, r2);
		f_gauss_fits_vtx_lof_y[i]->SetParameters(p0, p1, p2);

		h_lof_y[i]->Fit(Form("f_lof_y_%i", i), "Q0R");

		p0 = f_gauss_fits_vtx_lof_y[i]->GetParameter(0);
		p1 = f_gauss_fits_vtx_lof_y[i]->GetParameter(1);
		p2 = vertexRMS * f_gauss_fits_vtx_lof_y[i]->GetParameter(2);

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_lof_y[i]->SetParameters(p0, p1, p2);
		f_gauss_fits_vtx_lof_y[i]->SetRange(r1, r2);
		h_lof_y[i]->Fit(Form("f_lof_y_%i", i), "Q0R");

		// --> Synthetic LOF Y
		p0 = h_lof_synth_y[i]->GetMaximum();
		p1 = h_lof_synth_y[i]->GetMean();
		p2 = vertexRMS * h_lof_synth_y[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_lof_synth_y[i] = new TF1(Form("f_lof_synth_y_%i", i), "gaus", r1, r2);
		f_gauss_fits_vtx_lof_synth_y[i]->SetParameters(p0, p1, p2);

		h_lof_synth_y[i]->Fit(Form("f_lof_synth_y_%i", i), "Q0R");

		p0 = f_gauss_fits_vtx_lof_synth_y[i]->GetParameter(0);
		p1 = f_gauss_fits_vtx_lof_synth_y[i]->GetParameter(1);
		p2 = vertexRMS * f_gauss_fits_vtx_lof_synth_y[i]->GetParameter(2);

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_lof_synth_y[i]->SetParameters(p0, p1, p2);
		f_gauss_fits_vtx_lof_synth_y[i]->SetRange(r1, r2);
		h_lof_synth_y[i]->Fit(Form("f_lof_synth_y_%i", i), "Q0R");

		// --> LOF Y E
		p0 = h_lof_y_E[i]->GetMaximum();
		p1 = h_lof_y_E[i]->GetMean();
		p2 = vertexRMS * h_lof_y_E[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_lof_y_E[i] = new TF1(Form("f_lof_y_E_%i", i), "gaus", r1, r2);
		f_gauss_fits_vtx_lof_y_E[i]->SetParameters(p0, p1, p2);

		h_lof_y_E[i]->Fit(Form("f_lof_y_E_%i", i), "Q0R");

		p0 = f_gauss_fits_vtx_lof_y_E[i]->GetParameter(0);
		p1 = f_gauss_fits_vtx_lof_y_E[i]->GetParameter(1);
		p2 = vertexRMS * f_gauss_fits_vtx_lof_y_E[i]->GetParameter(2);

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_lof_y_E[i]->SetParameters(p0, p1, p2);
		f_gauss_fits_vtx_lof_y_E[i]->SetRange(r1, r2);
		h_lof_y_E[i]->Fit(Form("f_lof_y_E_%i", i), "Q0R");

		// --> LOF Y W
		p0 = h_lof_y_W[i]->GetMaximum();
		p1 = h_lof_y_W[i]->GetMean();
		p2 = vertexRMS * h_lof_y_W[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_lof_y_W[i] = new TF1(Form("f_lof_y_W_%i", i), "gaus", r1, r2);
		f_gauss_fits_vtx_lof_y_W[i]->SetParameters(p0, p1, p2);

		h_lof_y_W[i]->Fit(Form("f_lof_y_W_%i", i), "Q0R");

		p0 = f_gauss_fits_vtx_lof_y_W[i]->GetParameter(0);
		p1 = f_gauss_fits_vtx_lof_y_W[i]->GetParameter(1);
		p2 = vertexRMS * f_gauss_fits_vtx_lof_y_W[i]->GetParameter(2);

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_vtx_lof_y_W[i]->SetParameters(p0, p1, p2);
		f_gauss_fits_vtx_lof_y_W[i]->SetRange(r1, r2);
		h_lof_y_W[i]->Fit(Form("f_lof_y_W_%i", i), "Q0R");

		// --> LOF E-W X
		p0 = h_ew_lof_x[i]->GetMaximum();
		p1 = h_ew_lof_x[i]->GetMean();
		p2 = vertexDiffRMS * h_ew_lof_x[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_diff_lof_x[i] = new TF1(Form("f_ew_lof_x_%i", i), "gaus", r1, r2);
		f_gauss_fits_diff_lof_x[i]->SetParameters(p0, p1, p2);

		h_ew_lof_x[i]->Fit(Form("f_ew_lof_x_%i", i), "Q0R");

		p0 = f_gauss_fits_diff_lof_x[i]->GetParameter(0);
		p1 = f_gauss_fits_diff_lof_x[i]->GetParameter(1);
		p2 = vertexRMS * f_gauss_fits_diff_lof_x[i]->GetParameter(2);

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_diff_lof_x[i]->SetParameters(p0, p1, p2);
		f_gauss_fits_diff_lof_x[i]->SetRange(r1, r2);
		h_ew_lof_x[i]->Fit(Form("f_ew_lof_x_%i", i), "Q0R");

		// --> LOF E-W Y
		p0 = h_ew_lof_y[i]->GetMaximum();
		p1 = h_ew_lof_y[i]->GetMean();
		p2 = vertexDiffRMS * h_ew_lof_y[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_diff_lof_y[i] = new TF1(Form("f_ew_lof_y_%i", i), "gaus", r1, r2);
		f_gauss_fits_diff_lof_y[i]->SetParameters(p0, p1, p2);

		h_ew_lof_y[i]->Fit(Form("f_ew_lof_y_%i", i), "Q0R");

		p0 = f_gauss_fits_diff_lof_y[i]->GetParameter(0);
		p1 = f_gauss_fits_diff_lof_y[i]->GetParameter(1);
		p2 = vertexRMS * f_gauss_fits_diff_lof_y[i]->GetParameter(2);

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_diff_lof_y[i]->SetParameters(p0, p1, p2);
		f_gauss_fits_diff_lof_y[i]->SetRange(r1, r2);
		h_ew_lof_y[i]->Fit(Form("f_ew_lof_y_%i", i), "Q0R");

		// --> LOF E-W Z
		p0 = h_ew_lof_z[i]->GetMaximum();
		p1 = h_ew_lof_z[i]->GetMean();
		p2 = vertexDiffRMS * h_ew_lof_z[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_diff_lof_z[i] = new TF1(Form("f_ew_lof_z_%i", i), "gaus", r1, r2);
		f_gauss_fits_diff_lof_z[i]->SetParameters(p0, p1, p2);

		h_ew_lof_z[i]->Fit(Form("f_ew_lof_z_%i", i), "Q0R");

		p0 = f_gauss_fits_diff_lof_z[i]->GetParameter(0);
		p1 = f_gauss_fits_diff_lof_z[i]->GetParameter(1);
		p2 = vertexRMS * f_gauss_fits_diff_lof_z[i]->GetParameter(2);

		r1 = p1 - p2;
		r2 = p1 + p2;

		f_gauss_fits_diff_lof_z[i]->SetParameters(p0, p1, p2);
		f_gauss_fits_diff_lof_z[i]->SetRange(r1, r2);
		h_ew_lof_z[i]->Fit(Form("f_ew_lof_z_%i", i), "Q0R");
	}
}

void readHistograms()
{
	TFile *fin = new TFile("Data/423844_data_narrowvtx_ew.root");
	ntp_svxseg = (TTree*) fin->Get("ntp_event");

	for (int i = 0; i < NUMTRACKS; i++)
	{
		ntp_svxseg->Draw(Form("vtx_prec[0]>>h_prec_total_x_%i(150,-0.05,0.3)", i), Form("nvtxtrk_prec == %i", i + 2), "goff");
		h_prec_total_x[i] = (TH1F*) gDirectory->FindObject(Form("h_prec_total_x_%i", i));

		ntp_svxseg->Draw(Form("vtx_prec[1]>>h_prec_total_y_%i(150,-0.05,0.3)", i), Form("nvtxtrk_prec == %i", i + 2), "goff");
		h_prec_total_y[i] = (TH1F*) gDirectory->FindObject(Form("h_prec_total_y_%i", i));

		ntp_svxseg->Draw(Form("vtx_lof[0]>>h_lof_total_x_%i(150,-0.05,0.3)", i), Form("nvtxtrk_lof == %i", i + 2), "goff");
		h_lof_total_x[i] = (TH1F*) gDirectory->FindObject(Form("h_lof_total_x_%i", i));

		ntp_svxseg->Draw(Form("vtx_lof[1]>>h_lof_total_y_%i(150,-0.05,0.3)", i), Form("nvtxtrk_lof == %i", i + 2), "goff");
		h_lof_total_y[i] = (TH1F*) gDirectory->FindObject(Form("h_lof_total_y_%i", i));
	}

	for (int i = 0; i < NUMSEG; i++)
	{
		//Get vertex distribution from iterative algorithm
		ntp_svxseg->Draw(Form("vtx_prec[0]>>h_prec_x_%i(150,-0.05,0.3)", i), nseg_cuts_prec[i].c_str(), "goff");
		h_prec_x[i] = (TH1F*) gDirectory->FindObject(Form("h_prec_x_%i", i));

		ntp_svxseg->Draw(Form("vtx_prec[1]>>h_prec_y_%i(150,-0.05,0.3)", i), nseg_cuts_prec[i].c_str(), "goff");
		h_prec_y[i] = (TH1F*) gDirectory->FindObject(Form("h_prec_y_%i", i));

		ntp_svxseg->Draw(Form("vtx_prec[2]>>h_prec_z_%i(150,-20,20)", i), nseg_cuts_prec[i].c_str(), "goff");
		h_prec_z[i] = (TH1F*) gDirectory->FindObject(Form("h_prec_z_%i", i));

		ntp_svxseg->Draw(Form("vtx_prec_E[0]>>h_prec_x_E_%i(150,-0.05,0.3)", i), nseg_cuts_prec[i].c_str(), "goff");
		h_prec_x_E[i] = (TH1F*) gDirectory->FindObject(Form("h_prec_x_E_%i", i));

		ntp_svxseg->Draw(Form("vtx_prec_E[1]>>h_prec_y_E_%i(150,-0.05,0.3)", i), nseg_cuts_prec[i].c_str(), "goff");
		h_prec_y_E[i] = (TH1F*) gDirectory->FindObject(Form("h_prec_y_E_%i", i));

		ntp_svxseg->Draw(Form("vtx_prec_W[0]>>h_prec_x_W_%i(150,-0.05,0.3)", i), nseg_cuts_prec[i].c_str(), "goff");
		h_prec_x_W[i] = (TH1F*) gDirectory->FindObject(Form("h_prec_x_W_%i", i));

		ntp_svxseg->Draw(Form("vtx_prec_W[1]>>h_prec_y_W_%i(150,-0.05,0.3)", i), nseg_cuts_prec[i].c_str(), "goff");
		h_prec_y_W[i] = (TH1F*) gDirectory->FindObject(Form("h_prec_y_W_%i", i));

		ntp_svxseg->Draw(Form("vtx_prec_E[0] - vtx_prec_W[0]>>h_ew_prec_x_%i(150,-0.3,0.3)", i), nseg_cuts_prec[i].c_str(), "goff");
		h_ew_prec_x[i] = (TH1F*) gDirectory->FindObject(Form("h_ew_prec_x_%i", i));

		ntp_svxseg->Draw(Form("vtx_prec_E[1] - vtx_prec_W[1]>>h_ew_prec_y_%i(150,-0.3,0.3)", i), nseg_cuts_prec[i].c_str(), "goff");
		h_ew_prec_y[i] = (TH1F*) gDirectory->FindObject(Form("h_ew_prec_y_%i", i));

		ntp_svxseg->Draw(Form("vtx_prec_E[2] - vtx_prec_W[2]>>h_ew_prec_z_%i(150,-0.3,0.3)", i), nseg_cuts_prec[i].c_str(), "goff");
		h_ew_prec_z[i] = (TH1F*) gDirectory->FindObject(Form("h_ew_prec_z_%i", i));

		ntp_svxseg->Draw(Form("(vtx_prec_E[0] + vtx_prec_W[0])/2.0>>h_prec_synth_x_%i(150,-0.05,0.3)", i), nseg_cuts_prec[i].c_str(), "goff");
		h_prec_synth_x[i] = (TH1F*) gDirectory->FindObject(Form("h_prec_synth_x_%i", i));

		ntp_svxseg->Draw(Form("(vtx_prec_E[1] + vtx_prec_W[1])/2.0>>h_prec_synth_y_%i(150,-0.05,0.3)", i), nseg_cuts_prec[i].c_str(), "goff");
		h_prec_synth_y[i] = (TH1F*) gDirectory->FindObject(Form("h_prec_synth_y_%i", i));

		//Get vertex distributions from the LOF algorithm
		ntp_svxseg->Draw(Form("vtx_lof[0]>>h_lof_x_%i(150,-0.05,0.3)", i), nseg_cuts_lof[i].c_str(), "goff");
		h_lof_x[i] = (TH1F*) gDirectory->FindObject(Form("h_lof_x_%i", i));

		ntp_svxseg->Draw(Form("vtx_lof[1]>>h_lof_y_%i(150,-0.05,0.3)", i), nseg_cuts_lof[i].c_str(), "goff");
		h_lof_y[i] = (TH1F*) gDirectory->FindObject(Form("h_lof_y_%i", i));

		ntp_svxseg->Draw(Form("vtx_lof[2]>>h_lof_z_%i(150,-20,20)", i), nseg_cuts_lof[i].c_str(), "goff");
		h_lof_z[i] = (TH1F*) gDirectory->FindObject(Form("h_lof_z_%i", i));

		ntp_svxseg->Draw(Form("vtx_lof_E[0]>>h_lof_x_E_%i(150,-0.05,0.3)", i), nseg_cuts_lof[i].c_str(), "goff");
		h_lof_x_E[i] = (TH1F*) gDirectory->FindObject(Form("h_lof_x_E_%i", i));

		ntp_svxseg->Draw(Form("vtx_lof_E[1]>>h_lof_y_E_%i(150,-0.05,0.3)", i), nseg_cuts_lof[i].c_str(), "goff");
		h_lof_y_E[i] = (TH1F*) gDirectory->FindObject(Form("h_lof_y_E_%i", i));

		ntp_svxseg->Draw(Form("vtx_lof_W[0]>>h_lof_x_W_%i(150,-0.05,0.3)", i), nseg_cuts_lof[i].c_str(), "goff");
		h_lof_x_W[i] = (TH1F*) gDirectory->FindObject(Form("h_lof_x_W_%i", i));

		ntp_svxseg->Draw(Form("vtx_lof_W[1]>>h_lof_y_W_%i(150,-0.05,0.3)", i), nseg_cuts_lof[i].c_str(), "goff");
		h_lof_y_W[i] = (TH1F*) gDirectory->FindObject(Form("h_lof_y_W_%i", i));

		ntp_svxseg->Draw(Form("vtx_lof_E[0] - vtx_lof_W[0]>>h_ew_lof_x_%i(150,-0.3,0.3)", i), nseg_cuts_lof[i].c_str(), "goff");
		h_ew_lof_x[i] = (TH1F*) gDirectory->FindObject(Form("h_ew_lof_x_%i", i));

		ntp_svxseg->Draw(Form("vtx_lof_E[1] - vtx_lof_W[1]>>h_ew_lof_y_%i(150,-0.3,0.3)", i), nseg_cuts_lof[i].c_str(), "goff");
		h_ew_lof_y[i] = (TH1F*) gDirectory->FindObject(Form("h_ew_lof_y_%i", i));

		ntp_svxseg->Draw(Form("vtx_lof_E[2] - vtx_lof_W[2]>>h_ew_lof_z_%i(150,-0.3,0.3)", i), nseg_cuts_lof[i].c_str(), "goff");
		h_ew_lof_z[i] = (TH1F*) gDirectory->FindObject(Form("h_ew_lof_z_%i", i));

		ntp_svxseg->Draw(Form("(vtx_lof_E[0] + vtx_lof_W[0])/2.0>>h_lof_synth_x_%i(150,-0.05,0.3)", i), nseg_cuts_lof[i].c_str(), "goff");
		h_lof_synth_x[i] = (TH1F*) gDirectory->FindObject(Form("h_lof_synth_x_%i", i));

		ntp_svxseg->Draw(Form("(vtx_lof_E[1] + vtx_lof_W[1])/2.0>>h_lof_synth_y_%i(150,-0.05,0.3)", i), nseg_cuts_lof[i].c_str(), "goff");
		h_lof_synth_y[i] = (TH1F*) gDirectory->FindObject(Form("h_lof_synth_y_%i", i));
	}
}

void plotResolution()
{
	float x[NUMSEG];
	float x_aux[NUMSEG - 1];
	float err_x[NUMSEG] = {0};

	for (int i = 0; i < NUMSEG; i++)
	{
		x[i] = i + 1;

		if (i > 0)
		{
			x_aux[i - 1] = i + 1;
		}
	}

	float resol_prec_x[NUMSEG];
	float resol_prec_y[NUMSEG];
	float resol_prec_x_aux[NUMSEG - 1];
	float resol_prec_y_aux[NUMSEG - 1];
	float resol_prec_x_err[NUMSEG];
	float resol_prec_y_err[NUMSEG];

	float resol_lof_x[NUMSEG];
	float resol_lof_y[NUMSEG];
	float resol_lof_x_aux[NUMSEG - 1];
	float resol_lof_y_aux[NUMSEG - 1];
	float resol_lof_x_err[NUMSEG];
	float resol_lof_y_err[NUMSEG];

	for (int i = 0; i < NUMSEG; i++)
	{
		resol_prec_x[i] = resol_prec[0][i];
		resol_prec_y[i] = resol_prec[1][i];

		resol_lof_x[i] = resol_lof[0][i];
		resol_lof_y[i] = resol_lof[1][i];

		if (i > 0)
		{
			resol_prec_x_aux[i - 1] = resol_prec[0][i];
			resol_prec_y_aux[i - 1] = resol_prec[1][i];

			resol_lof_x_aux[i - 1] = resol_lof[0][i];
			resol_lof_y_aux[i - 1] = resol_lof[1][i];
		}

		resol_prec_x_err[i] = resol_prec_err[0][i];
		resol_prec_y_err[i] = resol_prec_err[1][i];

		resol_lof_x_err[i] = resol_lof_err[0][i];
		resol_lof_y_err[i] = resol_lof_err[1][i];
	}

	g_resol_prec_x = new TGraphErrors(NUMSEG, x, resol_prec_x, err_x, resol_prec_x_err);
	g_resol_prec_y = new TGraphErrors(NUMSEG, x, resol_prec_y, err_x, resol_prec_y_err);

	g_resol_lof_x = new TGraphErrors(NUMSEG, x, resol_lof_x, err_x, resol_lof_x_err);
	g_resol_lof_y = new TGraphErrors(NUMSEG, x, resol_lof_y, err_x, resol_lof_y_err);

	g_resol_prec_x_aux = new TGraphErrors(NUMSEG - 1, x_aux, resol_prec_x_aux);
	g_resol_prec_y_aux = new TGraphErrors(NUMSEG - 1, x_aux, resol_prec_y_aux);

	g_resol_lof_x_aux = new TGraphErrors(NUMSEG - 1, x_aux, resol_lof_x_aux);
	g_resol_lof_y_aux = new TGraphErrors(NUMSEG - 1, x_aux, resol_lof_y_aux);

	TCanvas *cResol = new TCanvas("cResol", "cResol", 900, 500);
	cResol->Divide(2, 1);
	cResol->cd(1);
	gPad->SetGridy();
	g_resol_prec_x->SetTitle("ITERATIVE");
	g_resol_prec_x->SetMarkerStyle(20);
	g_resol_prec_x->GetXaxis()->SetTitle("Number of Vertex Tracks in Each Arm");
	g_resol_prec_x->GetXaxis()->SetTitleFont(62);
	g_resol_prec_x->GetXaxis()->SetLabelFont(62);
	g_resol_prec_x->GetYaxis()->SetTitle("#sigma_{res} [cm]");
	g_resol_prec_x->GetYaxis()->SetTitleOffset(1.85);
	g_resol_prec_x->GetYaxis()->SetRangeUser(0, 0.025);
	g_resol_prec_x->GetYaxis()->SetTitleFont(62);
	g_resol_prec_x->GetYaxis()->SetLabelFont(62);
	g_resol_prec_x->Draw("AP");
	g_resol_prec_y->SetMarkerStyle(20);
	g_resol_prec_y->SetLineColor(kBlue);
	g_resol_prec_y->SetMarkerColor(kBlue);
	g_resol_prec_y->Draw("P,same");
	g_resol_prec_y_aux->SetMarkerStyle(20);
	g_resol_prec_y_aux->SetLineStyle(2);
	g_resol_prec_y_aux->SetLineColor(kBlue);
	g_resol_prec_y_aux->SetMarkerColor(kBlue);
	g_resol_prec_y_aux->Draw("LP,same");
	g_resol_prec_x_aux->SetMarkerStyle(20);
	g_resol_prec_x_aux->SetLineStyle(2);
	g_resol_prec_x_aux->SetLineColor(kBlack);
	g_resol_prec_x_aux->SetMarkerColor(kBlack);
	g_resol_prec_x_aux->Draw("LP,same");

	for (int i = 0; i < NUMSEG; i++)
	{
		g_resol_prec_x->GetXaxis()->SetBinLabel(g_resol_prec_x->GetXaxis()->FindBin(i + 1), bin_label[i].c_str());
		g_resol_prec_y->GetXaxis()->SetBinLabel(g_resol_prec_y->GetXaxis()->FindBin(i + 1), bin_label[i].c_str());

	}

	TLegend *tlg_resol_prec = new TLegend(0.6, 0.2, 0.85, 0.3);
	tlg_resol_prec->AddEntry(g_resol_prec_x, "X-Vertex", "LP");
	tlg_resol_prec->AddEntry(g_resol_prec_y, "Y-Vertex", "LP");
	tlg_resol_prec->SetLineColor(kWhite);
	tlg_resol_prec->SetTextSize(0.05);
	tlg_resol_prec->Draw("same");

	TLine *tlineDivAny = new TLine(1.5, 0, 1.5, 0.025);
	tlineDivAny->SetLineStyle(2);
	tlineDivAny->Draw("same");

	cResol->cd(2);
	gPad->SetGridy();
	g_resol_lof_x->SetTitle("LOF");
	g_resol_lof_x->SetMarkerStyle(20);
	g_resol_lof_x->GetXaxis()->SetTitle("Number of Vertex Tracks in Each Arm");
	g_resol_lof_x->GetXaxis()->SetTitleFont(62);
	g_resol_lof_x->GetXaxis()->SetLabelFont(62);
	g_resol_lof_x->GetYaxis()->SetTitle("#sigma_{res} [cm]");
	g_resol_lof_x->GetYaxis()->SetTitleOffset(1.85);
	g_resol_lof_x->GetYaxis()->SetRangeUser(0, 0.025);
	g_resol_lof_x->GetYaxis()->SetTitleFont(62);
	g_resol_lof_x->GetYaxis()->SetLabelFont(62);
	g_resol_lof_x->Draw("AP");
	g_resol_lof_y->SetMarkerStyle(20);
	g_resol_lof_y->SetLineColor(kBlue);
	g_resol_lof_y->SetMarkerColor(kBlue);
	g_resol_lof_y->Draw("P,same");
	g_resol_lof_y_aux->SetMarkerStyle(20);
	g_resol_lof_y_aux->SetLineStyle(2);
	g_resol_lof_y_aux->SetLineColor(kBlue);
	g_resol_lof_y_aux->SetMarkerColor(kBlue);
	g_resol_lof_y_aux->Draw("LP,same");
	g_resol_lof_x_aux->SetMarkerStyle(20);
	g_resol_lof_x_aux->SetLineStyle(2);
	g_resol_lof_x_aux->SetLineColor(kBlack);
	g_resol_lof_x_aux->SetMarkerColor(kBlack);
	g_resol_lof_x_aux->Draw("LP,same");

	for (int i = 0; i < NUMSEG; i++)
	{
		g_resol_lof_x->GetXaxis()->SetBinLabel(g_resol_lof_x->GetXaxis()->FindBin(i + 1), bin_label[i].c_str());
		g_resol_lof_y->GetXaxis()->SetBinLabel(g_resol_lof_y->GetXaxis()->FindBin(i + 1), bin_label[i].c_str());

	}

	tlineDivAny->Draw("same");

	//Plot resolution as a function of the total number of tracks with a fixed beam spread
	TCanvas *cTotalResolution = new TCanvas("cTotalResolution", "cTotalResolution", 900, 500);
	cTotalResolution->Divide(2, 1);

	cTotalResolution->cd(1);
	g_resol_prec_total_x->SetTitle("");
	g_resol_prec_total_x->SetMarkerStyle(21);
	g_resol_prec_total_x->SetMarkerColor(kRed);
	g_resol_prec_total_x->SetLineColor(kRed);
	g_resol_prec_total_x->GetXaxis()->SetTitle("Total Vertex Tracks");
	g_resol_prec_total_x->GetXaxis()->SetTitleFont(62);
	g_resol_prec_total_x->GetXaxis()->SetLabelFont(62);
	g_resol_prec_total_x->GetYaxis()->SetTitle("#sigma_{res,x} [cm]");
	g_resol_prec_total_x->GetYaxis()->SetTitleOffset(1.85);
	g_resol_prec_total_x->GetYaxis()->SetRangeUser(0, 0.04);
	g_resol_prec_total_x->GetYaxis()->SetTitleFont(62);
	g_resol_prec_total_x->GetYaxis()->SetLabelFont(62);
	g_resol_prec_total_x->Draw("ALP");

	g_resol_lof_total_x->SetMarkerStyle(21);
	g_resol_lof_total_x->SetMarkerColor(kBlue);
	g_resol_lof_total_x->SetLineColor(kBlue);
	g_resol_lof_total_x->GetXaxis()->SetTitleFont(62);
	g_resol_lof_total_x->GetXaxis()->SetLabelFont(62);
	g_resol_lof_total_x->GetYaxis()->SetRangeUser(0, 0.04);
	g_resol_lof_total_x->GetYaxis()->SetTitleFont(62);
	g_resol_lof_total_x->GetYaxis()->SetLabelFont(62);
	g_resol_lof_total_x->Draw("LP,same");

	cTotalResolution->cd(2);
	g_resol_prec_total_y->SetTitle("");
	g_resol_prec_total_y->SetMarkerStyle(21);
	g_resol_prec_total_y->SetMarkerColor(kRed);
	g_resol_prec_total_y->SetLineColor(kRed);
	g_resol_prec_total_y->GetXaxis()->SetTitle("Total Vertex Tracks");
	g_resol_prec_total_y->GetXaxis()->SetTitleFont(62);
	g_resol_prec_total_y->GetXaxis()->SetLabelFont(62);
	g_resol_prec_total_y->GetYaxis()->SetTitle("#sigma_{res,y} [cm]");
	g_resol_prec_total_y->GetYaxis()->SetTitleOffset(1.85);
	g_resol_prec_total_y->GetYaxis()->SetRangeUser(0, 0.04);
	g_resol_prec_total_y->GetYaxis()->SetTitleFont(62);
	g_resol_prec_total_y->GetYaxis()->SetLabelFont(62);
	g_resol_prec_total_y->Draw("ALP");

	g_resol_lof_total_y->SetMarkerStyle(21);
	g_resol_lof_total_y->SetMarkerColor(kBlue);
	g_resol_lof_total_y->SetLineColor(kBlue);
	g_resol_lof_total_y->GetXaxis()->SetTitleFont(62);
	g_resol_lof_total_y->GetXaxis()->SetLabelFont(62);
	g_resol_lof_total_y->GetYaxis()->SetRangeUser(0, 0.04);
	g_resol_lof_total_y->GetYaxis()->SetTitleFont(62);
	g_resol_lof_total_y->GetYaxis()->SetLabelFont(62);
	g_resol_lof_total_y->Draw("LP,same");
}

void plotBeamSpot()
{
	float x[NUMSEG];
	float x_aux[NUMSEG - 1];

	float err_x[NUMSEG] = {0};

	for (int i = 0; i < NUMSEG; i++)
	{
		x[i] = i + 1;

		if (i > 0)
		{
			x_aux[i - 1] = i + 1;
		}
	}

	float bspt_prec_x[NUMSEG];
	float bspt_prec_y[NUMSEG];
	float bspt_prec_x_2[NUMSEG];
	float bspt_prec_y_2[NUMSEG];
	float bspt_prec_x_aux[NUMSEG - 1];
	float bspt_prec_y_aux[NUMSEG - 1];

	float bspt_prec_er_x[NUMSEG];
	float bspt_prec_er_y[NUMSEG];

	float bspt_lof_x[NUMSEG];
	float bspt_lof_y[NUMSEG];
	float bspt_lof_x_2[NUMSEG];
	float bspt_lof_y_2[NUMSEG];
	float bspt_lof_x_aux[NUMSEG - 1];
	float bspt_lof_y_aux[NUMSEG - 1];

	float bspt_lof_er_x[NUMSEG];
	float bspt_lof_er_y[NUMSEG];

	for (int i = 0; i < NUMSEG; i++)
	{
		bspt_prec_x[i] = bspt_prec[0][i];
		bspt_prec_y[i] = bspt_prec[1][i];

		bspt_prec_x_2[i] = bspt_prec_2[0][i];
		bspt_prec_y_2[i] = bspt_prec_2[1][i];

		bspt_prec_er_x[i] = bspt_prec_err[0][i];
		bspt_prec_er_y[i] = bspt_prec_err[1][i];

		bspt_lof_x[i] = bspt_lof[0][i];
		bspt_lof_y[i] = bspt_lof[1][i];

		bspt_lof_x_2[i] = bspt_lof_2[0][i];
		bspt_lof_y_2[i] = bspt_lof_2[1][i];

		bspt_lof_er_x[i] = bspt_lof_err[0][i];
		bspt_lof_er_y[i] = bspt_lof_err[1][i];

		cout << "**** " <<  bspt_prec_er_x[i] << ", " << bspt_prec_er_y[i] << endl;

		if (i > 0)
		{
			bspt_prec_x_aux[i - 1] = bspt_prec[0][i];
			bspt_prec_y_aux[i - 1] = bspt_prec[1][i];

			bspt_lof_x_aux[i - 1] = bspt_lof[0][i];
			bspt_lof_y_aux[i - 1] = bspt_lof[1][i];
		}
	}

	g_bspt_prec_x = new TGraphErrors(NUMSEG, x, bspt_prec_x, err_x, bspt_prec_er_x);
	g_bspt_prec_y = new TGraphErrors(NUMSEG, x, bspt_prec_y, err_x, bspt_prec_er_y);

	g_bspt_prec_x_2 = new TGraphErrors(NUMSEG, x, bspt_prec_x_2);
	g_bspt_prec_y_2 = new TGraphErrors(NUMSEG, x, bspt_prec_y_2);

	g_bspt_lof_x = new TGraphErrors(NUMSEG, x, bspt_lof_x, err_x, bspt_lof_er_x);
	g_bspt_lof_y = new TGraphErrors(NUMSEG, x, bspt_lof_y, err_x, bspt_lof_er_y);

	g_bspt_lof_x_2 = new TGraphErrors(NUMSEG, x, bspt_lof_x_2);
	g_bspt_lof_y_2 = new TGraphErrors(NUMSEG, x, bspt_lof_y_2);

	g_bspt_prec_x_aux = new TGraphErrors(NUMSEG - 1, x_aux, bspt_prec_x_aux);
	g_bspt_prec_y_aux = new TGraphErrors(NUMSEG - 1, x_aux, bspt_prec_y_aux);

	g_bspt_lof_x_aux = new TGraphErrors(NUMSEG - 1, x_aux, bspt_lof_x_aux);
	g_bspt_lof_y_aux = new TGraphErrors(NUMSEG - 1, x_aux, bspt_lof_y_aux);


	TCanvas *cBspt = new TCanvas("cBspt", "cBspt", 900, 500);
	cBspt->Divide(2, 1);
	cBspt->cd(1);
	gPad->SetGridy();
	g_bspt_prec_x->SetTitle("ITERATIVE");
	g_bspt_prec_x->SetMarkerStyle(20);
	g_bspt_prec_x->GetXaxis()->SetTitle("Number of Vertex Tracks in Each Arm");
	g_bspt_prec_x->GetXaxis()->SetTitleFont(62);
	g_bspt_prec_x->GetXaxis()->SetLabelFont(62);
	g_bspt_prec_x->GetYaxis()->SetTitle("#sigma_{b} [cm]");
	g_bspt_prec_x->GetYaxis()->SetTitleOffset(1.85);
	g_bspt_prec_x->GetYaxis()->SetRangeUser(0, 0.025);
	g_bspt_prec_x->GetYaxis()->SetTitleFont(62);
	g_bspt_prec_x->GetYaxis()->SetLabelFont(62);
	g_bspt_prec_x->Draw("AP");

	g_bspt_prec_y->SetMarkerStyle(20);
	g_bspt_prec_y->SetLineColor(kBlue);
	g_bspt_prec_y->SetMarkerColor(kBlue);
	g_bspt_prec_y->Draw("P,same");

	g_bspt_prec_y_aux->SetMarkerStyle(20);
	g_bspt_prec_y_aux->SetLineStyle(2);
	g_bspt_prec_y_aux->SetLineColor(kBlue);
	g_bspt_prec_y_aux->SetMarkerColor(kBlue);
	g_bspt_prec_y_aux->Draw("LP,same");

	g_bspt_prec_x_aux->SetMarkerStyle(20);
	g_bspt_prec_x_aux->SetLineStyle(2);
	g_bspt_prec_x_aux->SetLineColor(kBlack);
	g_bspt_prec_x_aux->SetMarkerColor(kBlack);
	g_bspt_prec_x_aux->Draw("LP,same");

	g_bspt_prec_x_2->SetTitle("");
	g_bspt_prec_x_2->SetMarkerStyle(4);
	g_bspt_prec_x_2->Draw("P,same");


	for (int i = 0; i < NUMSEG; i++)
	{
		g_bspt_prec_x->GetXaxis()->SetBinLabel(g_bspt_prec_x->GetXaxis()->FindBin(i + 1), bin_label[i].c_str());
		g_bspt_prec_y->GetXaxis()->SetBinLabel(g_bspt_prec_y->GetXaxis()->FindBin(i + 1), bin_label[i].c_str());

	}

	TLegend *tlg_bspt_prec = new TLegend(0.6, 0.2, 0.85, 0.3);
	tlg_bspt_prec->AddEntry(g_bspt_prec_x, "X-Vertex", "LP");
	tlg_bspt_prec->AddEntry(g_bspt_prec_y, "Y-Vertex", "LP");
	tlg_bspt_prec->SetLineColor(kWhite);
	tlg_bspt_prec->SetTextSize(0.05);
	tlg_bspt_prec->Draw("same");

	TLine *tlineDivAny = new TLine(1.5, 0, 1.5, 0.025);
	tlineDivAny->SetLineStyle(2);
	tlineDivAny->Draw("same");

	cBspt->cd(2);
	gPad->SetGridy();
	g_bspt_lof_x->SetTitle("LOF");
	g_bspt_lof_x->SetMarkerStyle(20);
	g_bspt_lof_x->GetXaxis()->SetTitle("Number of Vertex Tracks in Each Arm");
	g_bspt_lof_x->GetXaxis()->SetTitleFont(62);
	g_bspt_lof_x->GetXaxis()->SetLabelFont(62);
	g_bspt_lof_x->GetYaxis()->SetTitle("#sigma_{b} [cm]");
	g_bspt_lof_x->GetYaxis()->SetTitleOffset(1.85);
	g_bspt_lof_x->GetYaxis()->SetRangeUser(0, 0.025);
	g_bspt_lof_x->GetYaxis()->SetTitleFont(62);
	g_bspt_lof_x->GetYaxis()->SetLabelFont(62);
	g_bspt_lof_x->Draw("AP");
	g_bspt_lof_y->SetMarkerStyle(20);
	g_bspt_lof_y->SetLineColor(kBlue);
	g_bspt_lof_y->SetMarkerColor(kBlue);
	g_bspt_lof_y->Draw("P,same");
	g_bspt_lof_y_aux->SetMarkerStyle(20);
	g_bspt_lof_y_aux->SetLineStyle(2);
	g_bspt_lof_y_aux->SetLineColor(kBlue);
	g_bspt_lof_y_aux->SetMarkerColor(kBlue);
	g_bspt_lof_y_aux->Draw("LP,same");
	g_bspt_lof_x_aux->SetMarkerStyle(20);
	g_bspt_lof_x_aux->SetLineStyle(2);
	g_bspt_lof_x_aux->SetLineColor(kBlack);
	g_bspt_lof_x_aux->SetMarkerColor(kBlack);
	g_bspt_lof_x_aux->Draw("LP,same");

	for (int i = 0; i < NUMSEG; i++)
	{
		g_bspt_lof_x->GetXaxis()->SetBinLabel(g_bspt_lof_x->GetXaxis()->FindBin(i + 1), bin_label[i].c_str());
		g_bspt_lof_y->GetXaxis()->SetBinLabel(g_bspt_lof_y->GetXaxis()->FindBin(i + 1), bin_label[i].c_str());
	}

	tlineDivAny->Draw("same");
}

void calculateTotalResolution()
{
	//Pick a fixed beam spread corresponding to the no-track-cut selection
	float sigmaBeamPrecX = bspt_prec[0][0];
	float sigmaBeamPrecY = bspt_prec[1][0];
	float sigmaBeamLOFX = bspt_lof[0][0];
	float sigmaBeamLOFY = bspt_lof[1][0];
	float trackVals[NUMTRACKS];

	//Calculate the resolution as a function of the total number of tracks
	for (int i = 0; i < NUMTRACKS; i++)
	{
		float s_x = f_gauss_fits_vtx_prec_total_x[i]->GetParameter(2);
		float s_y = f_gauss_fits_vtx_prec_total_y[i]->GetParameter(2);
		resol_prec_total_x[i] = TMath::Sqrt(s_x * s_x - sigmaBeamPrecX * sigmaBeamPrecX);
		resol_prec_total_y[i] = TMath::Sqrt(s_y * s_y - sigmaBeamPrecY * sigmaBeamPrecY);

		s_x = f_gauss_fits_vtx_lof_total_x[i]->GetParameter(2);
		s_y = f_gauss_fits_vtx_lof_total_y[i]->GetParameter(2);
		resol_lof_total_x[i] = TMath::Sqrt(s_x * s_x - sigmaBeamLOFX * sigmaBeamLOFX);
		resol_lof_total_y[i] = TMath::Sqrt(s_y * s_y - sigmaBeamLOFY * sigmaBeamLOFY);

		trackVals[i] = i + 2;
	}

	g_resol_prec_total_x = new TGraphErrors(NUMTRACKS, trackVals, resol_prec_total_x);
	g_resol_prec_total_y = new TGraphErrors(NUMTRACKS, trackVals, resol_prec_total_y);

	g_resol_lof_total_x = new TGraphErrors(NUMTRACKS, trackVals, resol_lof_total_x);
	g_resol_lof_total_y = new TGraphErrors(NUMTRACKS, trackVals, resol_lof_total_y);
}

void plotTrackFrequency()
{
	float trackVals[NUMTRACKS];

	for(int i=0; i<NUMTRACKS; i++)
	{
		trackVals[i] = i+2;
		eventTrackFrequency[i] = (float) ntp_svxseg->GetEntries(Form("nvtxtrk_prec == %i",i+2)) / ntp_svxseg->GetEntries();
	}

	g_track_frequency = new TGraphErrors(NUMTRACKS, trackVals, eventTrackFrequency);

	TCanvas *cTrackFrequency = new TCanvas("cTrackFrequency", "cTrackFrequency", 600, 600);
	g_track_frequency->SetTitle("");
	g_track_frequency->SetMarkerStyle(20);
	g_track_frequency->GetXaxis()->SetTitle("Total Vertex Tracks");
	g_track_frequency->GetXaxis()->SetTitleFont(62);
	g_track_frequency->GetXaxis()->SetLabelFont(62);
	g_track_frequency->GetYaxis()->SetTitle("Fraction of Events");
	g_track_frequency->GetYaxis()->SetTitleOffset(1.85);
	g_track_frequency->GetYaxis()->SetRangeUser(0, 0.025);
	g_track_frequency->GetYaxis()->SetTitleFont(62);
	g_track_frequency->GetYaxis()->SetLabelFont(62);
	g_track_frequency->Draw("ALP");
}

void printParameters()
{
	//With no track requirement
	cout << "-----------------------------------------" << endl;
	cout << " NO TRACK REQUIREMENT" << endl;
	cout << "-----------------------------------------" << endl << endl;

	cout << "**BC_prec_x  = " << 10000 * f_gauss_fits_vtx_prec_x[0]->GetParameter(1) << " +/- " << 10000 * f_gauss_fits_vtx_prec_x[0]->GetParError(1) << endl;
	cout << "**BC_prec_y  = " << 10000 * f_gauss_fits_vtx_prec_y[0]->GetParameter(1) << " +/- " << 10000 * f_gauss_fits_vtx_prec_y[0]->GetParError(1) << endl << endl;

	cout << "**BC_lof_x  = " << 10000 * f_gauss_fits_vtx_lof_x[0]->GetParameter(1) << " +/- " << 10000 * f_gauss_fits_vtx_lof_x[0]->GetParError(1) << endl;
	cout << "**BC_lof_y  = " << 10000 * f_gauss_fits_vtx_lof_y[0]->GetParameter(1) << " +/- " << 10000 * 10000 * f_gauss_fits_vtx_lof_y[0]->GetParError(1) << endl << endl;

	cout << "s_x_prec     = " << 10000 * f_gauss_fits_vtx_prec_x[0]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_prec_x[0]->GetParError(2) << endl;
	cout << "s_x_prec_E   = " << 10000 * f_gauss_fits_vtx_prec_x_E[0]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_prec_x_E[0]->GetParError(2) << endl;
	cout << "s_x_prec_W   = " << 10000 * f_gauss_fits_vtx_prec_x_W[0]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_prec_x_W[0]->GetParError(2) << endl;
	cout << "s_x_prec_E-W = " << 10000 * f_gauss_fits_diff_prec_x[0]->GetParameter(2) << " +/- " << 10000 * 10000 * f_gauss_fits_vtx_prec_synth_x[0]->GetParError(2) << endl;
	cout << "s_x_prec_Syn = " << 10000 * f_gauss_fits_vtx_prec_synth_x[0]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_prec_synth_x[0]->GetParError(2) << endl << endl;

	cout << "s_x_lof     = " << 10000 * f_gauss_fits_vtx_lof_x[0]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_lof_x[0]->GetParError(2) << endl;
	cout << "s_x_lof_E   = " << 10000 * f_gauss_fits_vtx_lof_x_E[0]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_lof_x_E[0]->GetParError(2) << endl;
	cout << "s_x_lof_W   = " << 10000 * f_gauss_fits_vtx_lof_x_W[0]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_lof_x_W[0]->GetParError(2) << endl;
	cout << "s_x_lof_E-W = " << 10000 * f_gauss_fits_diff_lof_x[0]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_diff_lof_x[0]->GetParError(2) << endl;
	cout << "s_x_lof_Syn = " << 10000 * f_gauss_fits_vtx_lof_synth_x[0]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_lof_synth_x[0]->GetParError(2) << endl << endl;

	cout << "s_y_prec     = " << 10000 * f_gauss_fits_vtx_prec_y[0]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_prec_y[0]->GetParError(2) << endl;
	cout << "s_y_prec_E   = " << 10000 * f_gauss_fits_vtx_prec_y_E[0]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_prec_y_E[0]->GetParError(2) << endl;
	cout << "s_y_prec_W   = " << 10000 * f_gauss_fits_vtx_prec_y_W[0]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_prec_y_W[0]->GetParError(2) << endl;
	cout << "s_y_prec_E-W = " << 10000 * f_gauss_fits_diff_prec_y[0]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_prec_synth_y[0]->GetParError(2) << endl;
	cout << "s_y_prec_Syn = " << 10000 * f_gauss_fits_vtx_prec_synth_y[0]->GetParameter(2) << " +/- " << f_gauss_fits_vtx_prec_synth_y[0]->GetParError(2) << endl << endl;

	cout << "s_y_lof     = " << 10000 * f_gauss_fits_vtx_lof_y[0]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_lof_y[0]->GetParError(2) << endl;
	cout << "s_y_lof_E   = " << 10000 * f_gauss_fits_vtx_lof_y_E[0]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_lof_y_E[0]->GetParError(2) << endl;
	cout << "s_y_lof_W   = " << 10000 * f_gauss_fits_vtx_lof_y_W[0]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_lof_y_W[0]->GetParError(2) << endl;
	cout << "s_y_lof_E-W = " << 10000 * f_gauss_fits_diff_lof_y[0]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_diff_lof_y[0]->GetParError(2) << endl;
	cout << "s_y_lof_Syn = " << 10000 * f_gauss_fits_vtx_lof_synth_y[0]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_lof_synth_y[0]->GetParError(2) << endl << endl;

	cout << "BS_x_prec        = " << 10000 * bspt_prec[0][0] << endl;
	cout << "BS_y_prec        = " << 10000 * bspt_prec[1][0] << endl;
	cout << "RE_x_prec        = " << 10000 * resol_prec[0][0] << endl;
	cout << "RE_y_prec        = " << 10000 * resol_prec[1][0] << endl << endl;

	cout << "BS_x_lof        = " << 10000 * bspt_lof[0][0] << endl;
	cout << "BS_y_lof        = " << 10000 * bspt_lof[1][0] << endl;
	cout << "RE_x_lof        = " << 10000 * resol_lof[0][0] << endl;
	cout << "RE_y_lof        = " << 10000 * resol_lof[1][0] << endl << endl;

	cout << "BC_x_prec       = " << 10000 * f_gauss_fits_vtx_prec_x[0]->GetParameter(1) << endl;
	cout << "BC_y_prec       = " << 10000 * f_gauss_fits_vtx_prec_y[0]->GetParameter(1) << endl;
	cout << "BC_x_lof       = " << 10000 * f_gauss_fits_vtx_lof_x[0]->GetParameter(1) << endl;
	cout << "BC_y_lof       = " << 10000 * f_gauss_fits_vtx_lof_y[0]->GetParameter(1) << endl << endl;

	cout << "-----------------------------------------" << endl;
	cout << " TRACK REQUIREMENT = 2" << endl;
	cout << "-----------------------------------------" << endl << endl;

	cout << "**BC_prec_x  = " << 10000 * f_gauss_fits_vtx_prec_x[1]->GetParameter(1) << " +/- " << 10000 * f_gauss_fits_vtx_prec_x[1]->GetParError(1) << endl;
	cout << "**BC_prec_y  = " << 10000 * f_gauss_fits_vtx_prec_y[1]->GetParameter(1) << " +/- " << 10000 * f_gauss_fits_vtx_prec_y[1]->GetParError(1) << endl << endl;

	cout << "**BC_lof_x  = " << 10000 * f_gauss_fits_vtx_lof_x[1]->GetParameter(1) << " +/- " << 10000 * f_gauss_fits_vtx_lof_x[1]->GetParError(1) << endl;
	cout << "**BC_lof_y  = " << 10000 * f_gauss_fits_vtx_lof_y[1]->GetParameter(1) << " +/- " << 10000 * f_gauss_fits_vtx_lof_y[1]->GetParError(1) << endl << endl;

	cout << "s_x_prec     = " << 10000 * f_gauss_fits_vtx_prec_x[1]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_prec_x[1]->GetParError(2) << endl;
	cout << "s_x_prec_E   = " << 10000 * f_gauss_fits_vtx_prec_x_E[1]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_prec_x_E[1]->GetParError(2) << endl;
	cout << "s_x_prec_W   = " << 10000 * f_gauss_fits_vtx_prec_x_W[1]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_prec_x_W[1]->GetParError(2) << endl;
	cout << "s_x_prec_E-W = " << 10000 * f_gauss_fits_diff_prec_x[1]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_prec_synth_x[1]->GetParError(2) << endl;
	cout << "s_x_prec_Syn = " << 10000 * f_gauss_fits_vtx_prec_synth_x[1]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_prec_synth_x[1]->GetParError(2) << endl << endl;

	cout << "s_x_lof     = " << 10000 * f_gauss_fits_vtx_lof_x[1]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_lof_x[1]->GetParError(2) << endl;
	cout << "s_x_lof_E   = " << 10000 * f_gauss_fits_vtx_lof_x_E[1]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_lof_x_E[1]->GetParError(2) << endl;
	cout << "s_x_lof_W   = " << 10000 * f_gauss_fits_vtx_lof_x_W[1]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_lof_x_W[1]->GetParError(2) << endl;
	cout << "s_x_lof_E-W = " << 10000 * f_gauss_fits_diff_lof_x[1]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_diff_lof_x[1]->GetParError(2) << endl;
	cout << "s_x_lof_Syn = " << 10000 * f_gauss_fits_vtx_lof_synth_x[1]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_lof_synth_x[1]->GetParError(2) << endl << endl;

	cout << "s_y_prec     = " << 10000 * f_gauss_fits_vtx_prec_y[1]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_prec_y[1]->GetParError(2) << endl;
	cout << "s_y_prec_E   = " << 10000 * f_gauss_fits_vtx_prec_y_E[1]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_prec_y_E[1]->GetParError(2) << endl;
	cout << "s_y_prec_W   = " << 10000 * f_gauss_fits_vtx_prec_y_W[1]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_prec_y_W[1]->GetParError(2) << endl;
	cout << "s_y_prec_E-W = " << 10000 * f_gauss_fits_diff_prec_y[1]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_prec_synth_y[1]->GetParError(2) << endl;
	cout << "s_y_prec_Syn = " << 10000 * f_gauss_fits_vtx_prec_synth_y[1]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_prec_synth_y[1]->GetParError(2) << endl << endl;

	cout << "s_y_lof     = " << 10000 * f_gauss_fits_vtx_lof_y[1]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_lof_y[1]->GetParError(2) << endl;
	cout << "s_y_lof_E   = " << 10000 * f_gauss_fits_vtx_lof_y_E[1]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_lof_y_E[1]->GetParError(2) << endl;
	cout << "s_y_lof_W   = " << 10000 * f_gauss_fits_vtx_lof_y_W[1]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_lof_y_W[1]->GetParError(2) << endl;
	cout << "s_y_lof_E-W = " << 10000 * f_gauss_fits_diff_lof_y[1]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_diff_lof_y[1]->GetParError(2) << endl;
	cout << "s_y_lof_Syn = " << 10000 * f_gauss_fits_vtx_lof_synth_y[1]->GetParameter(2) << " +/- " << 10000 * f_gauss_fits_vtx_lof_synth_y[1]->GetParError(2) << endl << endl;

	cout << "BS_x_prec        = " << 10000 * bspt_prec[0][1] << endl;
	cout << "BS_y_prec        = " << 10000 * bspt_prec[1][1] << endl;
	cout << "RE_x_prec        = " << 10000 * resol_prec[0][1] << endl;
	cout << "RE_y_prec        = " << 10000 * resol_prec[1][1] << endl << endl;

	cout << "BS_x_lof        = " << 10000 * bspt_lof[0][1] << endl;
	cout << "BS_y_lof        = " << 10000 * bspt_lof[1][1] << endl;
	cout << "RE_x_lof        = " << 10000 * resol_lof[0][1] << endl;
	cout << "RE_y_lof        = " << 10000 * resol_lof[1][1] << endl << endl;

	cout << "BC_x_prec       = " << 10000 * f_gauss_fits_vtx_prec_x[1]->GetParameter(1) << endl;
	cout << "BC_y_prec       = " << 10000 * f_gauss_fits_vtx_prec_y[1]->GetParameter(1) << endl;
	cout << "BC_x_lof       = " << 10000 * f_gauss_fits_vtx_lof_x[1]->GetParameter(1) << endl;
	cout << "BC_y_lof       = " << 10000 * f_gauss_fits_vtx_lof_y[1]->GetParameter(1) << endl << endl;
}

void savePlots()
{
	TFile *fout = new TFile("output_data.root","RECREATE");
	g_track_frequency->SetName("g_track_frequency_prec");
	g_resol_prec_total_x->SetName("g_resol_prec_total_x");
	g_resol_prec_total_y->SetName("g_resol_prec_total_y");
	g_track_frequency->Write();
	g_resol_prec_total_x->Write();
	g_resol_prec_total_y->Write();
	fout->Close();
}

void ComputeBeamSpot()
{
	readHistograms();             //Create histograms from TTree
	fitHistograms();              //Fit histograms with a Gaussian
	calculateResolution();        //Calculate the resolution and beam spot size for selected E/W track combinations
	calculateTotalResolution();   //Calculate the resolution as a function of total number of beam tracks, with a fixed beam spread
	plotHistograms();             //Plot gaussian histograms and fits
	plotResolution();             //Plot resolution for E/W track combinations
	plotBeamSpot();               //Plot beam spread for E/W track combinations
	plotTrackFrequency();
	printParameters();
	savePlots();
}