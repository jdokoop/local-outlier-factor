//--------------------------------------------
// Macro to test the LOF vertexing algorithm
// by generating the following set of plots
// discriminating by number of segments,
// namely 2, 3, 4, 5-8, and > 8
//
// 1) x_lof - x_true
// 2) y_lof - y_true
// 3) z_lof - z_true
// 4) (x_lof - x_true) / sigma_x
// 5) (y_lof - y_true) / sigma_y
// 6) (z_lof - z_true) / sigma_z
//--------------------------------------------

#include <iostream>

using namespace std;

//--------------------------------------------
// Variables
//--------------------------------------------

//Number of segment categories
const int NCUTS = 6;

//Tree with event-level vertex and BBC multiplicity information
TTree *ntp_event;

//Distribution of x,y,z-vertices from pisa (true vertices)
TH1F *h_pisavtx_x[NCUTS];
TH1F *h_pisavtx_y[NCUTS];
TH1F *h_pisavtx_z[NCUTS];

//Distribution of x,y,z-vertices from LOF algorithm
TH1F *h_lofvtx_x[NCUTS];
TH1F *h_lofvtx_y[NCUTS];
TH1F *h_lofvtx_z[NCUTS];

//Distribution of x,y,z-vertices from precise vertex algorithm
TH1F *h_precvtx_x[NCUTS];
TH1F *h_precvtx_y[NCUTS];
TH1F *h_precvtx_z[NCUTS];

//Distribution of the difference between LOF and PISA x,y,z-vertices
TH1F *h_lofvtx_diff_x[NCUTS];
TH1F *h_lofvtx_diff_y[NCUTS];
TH1F *h_lofvtx_diff_z[NCUTS];

//Distribution of the difference between precise and PISA x,y,z-vertices
TH1F *h_precvtx_diff_x[NCUTS];
TH1F *h_precvtx_diff_y[NCUTS];
TH1F *h_precvtx_diff_z[NCUTS];

//Cut on the number of segment variable in ntp_event
string segmentCut[NCUTS] = {"nseg > 1", "nseg == 2", "nseg == 3", "nseg == 4", "nseg > 4 && nseg <= 8", "nseg >= 9"};

//Labels for segment cuts
string segmentCutLabel[NCUTS] = {"SEG > 1", "SEG = 2", "SEG = 3", "SEG = 4", "4 < SEG < 9", "SEG > 9"};

//Cut to restrict to events with reconstructed LOF and PRECISE vertex
string cutvtxLOFx = "TRUE";//"vtx_lof[0] > 0.05 && vtx_lof[0] < 0.25";
string cutvtxLOFy = "TRUE";//"vtx_lof[1] > -0.025 && vtx_lof[1] < 0.2";
string cutvtxPRECx = "TRUE";//"vtx_prec[0] > 0.05 && vtx_prec[0] < 0.25";
string cutvtxPRECy = "TRUE";//"vtx_prec[1] > -0.025 && vtx_prec[1] < 0.2";

//--------------------------------------------
// Functions
//--------------------------------------------

void fitResolutionHistograms()
{
	//Fit the difference between the reconstructed and the pisa vertices with a Gaussian
	double r1, r2;
	double p0, p1, p2;
	TF1 *fSignal;

	for (int i = 0; i < NCUTS; i++)
	{
		// --> LOF X
		p0 = h_lofvtx_diff_x[i]->GetMaximum();
		p1 = h_lofvtx_diff_x[i]->GetBinCenter(h_lofvtx_diff_x[i]->GetMaximumBin());
		p2 = 0.9 * h_lofvtx_diff_x[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		fSignal = new TF1(Form("f_lofvtx_diff_x_%i", i), "gaus", r1, r2);
		fSignal->SetParameters(p0, p1, p2);

		h_lofvtx_diff_x[i]->Fit(Form("f_lofvtx_diff_x_%i", i), "Q0R");

		// --> LOF Y
		p0 = h_lofvtx_diff_y[i]->GetMaximum();
		p1 = h_lofvtx_diff_y[i]->GetBinCenter(h_lofvtx_diff_y[i]->GetMaximumBin());
		p2 = 0.9 * h_lofvtx_diff_y[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		fSignal = new TF1(Form("f_lofvtx_diff_y_%i", i), "gaus", r1, r2);
		fSignal->SetParameters(p0, p1, p2);

		h_lofvtx_diff_y[i]->Fit(Form("f_lofvtx_diff_y_%i", i), "Q0R");

		// --> LOF Z
		p0 = h_lofvtx_diff_z[i]->GetMaximum();
		p1 = h_lofvtx_diff_z[i]->GetBinCenter(h_lofvtx_diff_z[i]->GetMaximumBin());
		p2 = 0.9 * h_lofvtx_diff_z[i]->GetRMS();

		r1 = p1 - p2;
		r2 = p1 + p2;

		fSignal = new TF1(Form("f_lofvtx_diff_z_%i", i), "gaus", r1, r2);
		fSignal->SetParameters(p0, p1, p2);

		h_lofvtx_diff_z[i]->Fit(Form("f_lofvtx_diff_z_%i", i), "Q0R");
	}
}

void normalizeHistograms()
{
	for (int i = 0; i < NCUTS; i++)
	{
		h_lofvtx_diff_x[i]->Scale(1.0 / h_lofvtx_diff_x[i]->GetMaximum());//h_lofvtx_diff_x[i]->Integral());
		h_precvtx_diff_x[i]->Scale(1.0 / h_precvtx_diff_x[i]->GetMaximum());//h_precvtx_diff_x[i]->Integral());

		h_lofvtx_diff_y[i]->Scale(1.0 / h_lofvtx_diff_y[i]->GetMaximum());//h_lofvtx_diff_y[i]->Integral());
		h_precvtx_diff_y[i]->Scale(1.0 / h_precvtx_diff_y[i]->GetMaximum());//h_precvtx_diff_y[i]->Integral());

		h_lofvtx_diff_z[i]->Scale(1.0 / h_lofvtx_diff_z[i]->GetMaximum());//h_lofvtx_diff_z[i]->Integral());
		h_precvtx_diff_z[i]->Scale(1.0 / h_precvtx_diff_z[i]->GetMaximum());//h_precvtx_diff_z[i]->Integral());
	}
}

void drawVertexDifferenceDistributions()
{
	gStyle->SetOptStat(0);

	TCanvas *cDiff[NCUTS];
	TLegend *legDiff[NCUTS];

	TLatex *tlatMeanLOFX[NCUTS];
	TLatex *tlatMeanLOFY[NCUTS];
	TLatex *tlatMeanLOFZ[NCUTS];

	TLatex *tlatMeanPRECX[NCUTS];
	TLatex *tlatMeanPRECY[NCUTS];
	TLatex *tlatMeanPRECZ[NCUTS];

	TLatex *tlatFracLOFX[NCUTS];
	TLatex *tlatFracLOFY[NCUTS];
	TLatex *tlatFracLOFZ[NCUTS];

	TLatex *tlatFracPRECX[NCUTS];
	TLatex *tlatFracPRECY[NCUTS];
	TLatex *tlatFracPRECZ[NCUTS];

	TLatex *tlatRMSLOFX[NCUTS];
	TLatex *tlatRMSLOFY[NCUTS];
	TLatex *tlatRMSLOFZ[NCUTS];

	TLatex *tlatRMSPRECX[NCUTS];
	TLatex *tlatRMSPRECY[NCUTS];
	TLatex *tlatRMSPRECZ[NCUTS];

	TF1 *fLOFDiff[NCUTS];

	for (int i = 0; i < NCUTS; i++)
	{
		cDiff[i] = new TCanvas(Form("c_vertex_%i", i), Form("Vertex Difference %s", segmentCutLabel[i].c_str()), 1100, 500);

		cDiff[i]->Divide(3, 1);
		cDiff[i]->cd(1);
		//gPad->SetLogy();

		h_lofvtx_diff_x[i]->GetXaxis()->SetTitle("(reco - pisa)_x [cm]");
		h_lofvtx_diff_x[i]->GetXaxis()->SetTitleFont(62);
		h_lofvtx_diff_x[i]->GetXaxis()->SetLabelFont(62);
		h_lofvtx_diff_x[i]->GetYaxis()->SetRangeUser(0, 1.2);
		h_lofvtx_diff_x[i]->GetYaxis()->SetTitle("AU");
		h_lofvtx_diff_x[i]->GetYaxis()->SetTitleFont(62);
		h_lofvtx_diff_x[i]->GetYaxis()->SetLabelFont(62);
		h_lofvtx_diff_x[i]->Draw();
		h_lofvtx_diff_x[i]->SetTitle("");
		h_precvtx_diff_x[i]->SetLineColor(kRed);
		h_precvtx_diff_x[i]->Draw("same");

		fLOFDiff[i] = (TF1*) h_lofvtx_diff_x[i]->GetFunction(Form("f_lofvtx_diff_x_%i", i));
		fLOFDiff[i]->SetLineColor(kBlack);
		//fLOFDiff[i]->SetLineWidth(2);
		fLOFDiff[i]->Draw("same");
		cout << "LOF RMS " << i << ": " << fLOFDiff[i]->GetParameter(2) << endl;

		tlatMeanLOFX[i] = new TLatex(0.15, 0.8, Form("Mean = %.3g", h_lofvtx_diff_x[i]->GetMean()));
		tlatMeanLOFX[i]->SetTextColor(kBlue);
		tlatMeanLOFX[i]->SetTextSize(0.035);
		tlatMeanLOFX[i]->SetNDC(kTRUE);
		tlatMeanLOFX[i]->Draw("same");

		tlatMeanPRECX[i] = new TLatex(0.15, 0.77, Form("Mean = %.3g", h_precvtx_diff_x[i]->GetMean()));
		tlatMeanPRECX[i]->SetTextColor(kRed);
		tlatMeanPRECX[i]->SetTextSize(0.035);
		tlatMeanPRECX[i]->SetNDC(kTRUE);
		tlatMeanPRECX[i]->Draw("same");

		tlatRMSLOFX[i] = new TLatex(0.15, 0.6, Form("RMS = %.3g", h_lofvtx_diff_x[i]->GetRMS()));
		tlatRMSLOFX[i]->SetTextColor(kBlue);
		tlatRMSLOFX[i]->SetTextSize(0.035);
		tlatRMSLOFX[i]->SetNDC(kTRUE);
		tlatRMSLOFX[i]->Draw("same");

		tlatRMSPRECX[i] = new TLatex(0.15, 0.57, Form("RMS = %.3g", h_precvtx_diff_x[i]->GetRMS()));
		tlatRMSPRECX[i]->SetTextColor(kRed);
		tlatRMSPRECX[i]->SetTextSize(0.035);
		tlatRMSPRECX[i]->SetNDC(kTRUE);
		tlatRMSPRECX[i]->Draw("same");

		tlatFracLOFX[i] = new TLatex(0.15, 0.2, Form("Frac Within +/- 0.05 = %.3g", h_lofvtx_diff_x[i]->Integral(h_lofvtx_diff_x[i]->FindBin(-0.05), h_lofvtx_diff_x[i]->FindBin(0.05)) / h_lofvtx_diff_x[i]->Integral()));
		tlatFracLOFX[i]->SetTextColor(kBlue);
		tlatFracLOFX[i]->SetTextSize(0.035);
		tlatFracLOFX[i]->SetNDC(kTRUE);
		tlatFracLOFX[i]->Draw("same");

		tlatFracPRECX[i] = new TLatex(0.15, 0.17, Form("Frac Within +/- 0.05 = %.3g", h_precvtx_diff_x[i]->Integral(h_precvtx_diff_x[i]->FindBin(-0.05), h_precvtx_diff_x[i]->FindBin(0.05)) / h_precvtx_diff_x[i]->Integral()));
		tlatFracPRECX[i]->SetTextColor(kRed);
		tlatFracPRECX[i]->SetTextSize(0.035);
		tlatFracPRECX[i]->SetNDC(kTRUE);
		tlatFracPRECX[i]->Draw("same");

		legDiff[i] = new TLegend(0.6, 0.7, 0.85, 0.8);
		legDiff[i]->AddEntry(h_lofvtx_diff_x[i], "LOF", "L");
		legDiff[i]->AddEntry(h_precvtx_diff_x[i], "PRECISE", "L");
		legDiff[i]->SetLineColor(kWhite);
		legDiff[i]->SetTextSize(0.05);
		legDiff[i]->Draw("same");

		cDiff[i]->cd(2);
		//gPad->SetLogy();

		h_lofvtx_diff_y[i]->GetXaxis()->SetTitle("(reco - pisa)_y [cm]");
		h_lofvtx_diff_y[i]->GetXaxis()->SetTitleFont(62);
		h_lofvtx_diff_y[i]->GetXaxis()->SetLabelFont(62);
		h_lofvtx_diff_y[i]->GetYaxis()->SetRangeUser(0, 1.2);
		h_lofvtx_diff_y[i]->GetYaxis()->SetTitle("AU");
		h_lofvtx_diff_y[i]->GetYaxis()->SetTitleFont(62);
		h_lofvtx_diff_y[i]->GetYaxis()->SetLabelFont(62);
		h_lofvtx_diff_y[i]->Draw();
		h_lofvtx_diff_y[i]->SetTitle("");
		h_precvtx_diff_y[i]->SetLineColor(kRed);
		h_precvtx_diff_y[i]->Draw("same");

		tlatMeanLOFY[i] = new TLatex(0.15, 0.8, Form("Mean = %.3g", h_lofvtx_diff_y[i]->GetMean()));
		tlatMeanLOFY[i]->SetTextColor(kBlue);
		tlatMeanLOFY[i]->SetTextSize(0.035);
		tlatMeanLOFY[i]->SetNDC(kTRUE);
		tlatMeanLOFY[i]->Draw("same");

		tlatMeanPRECY[i] = new TLatex(0.15, 0.77, Form("Mean = %.3g", h_precvtx_diff_y[i]->GetMean()));
		tlatMeanPRECY[i]->SetTextColor(kRed);
		tlatMeanPRECY[i]->SetTextSize(0.035);
		tlatMeanPRECY[i]->SetNDC(kTRUE);
		tlatMeanPRECY[i]->Draw("same");

		tlatRMSLOFY[i] = new TLatex(0.15, 0.6, Form("RMS = %.3g", h_lofvtx_diff_y[i]->GetRMS()));
		tlatRMSLOFY[i]->SetTextColor(kBlue);
		tlatRMSLOFY[i]->SetTextSize(0.035);
		tlatRMSLOFY[i]->SetNDC(kTRUE);
		tlatRMSLOFY[i]->Draw("same");

		tlatRMSPRECY[i] = new TLatex(0.15, 0.57, Form("RMS = %.3g", h_precvtx_diff_y[i]->GetRMS()));
		tlatRMSPRECY[i]->SetTextColor(kRed);
		tlatRMSPRECY[i]->SetTextSize(0.035);
		tlatRMSPRECY[i]->SetNDC(kTRUE);
		tlatRMSPRECY[i]->Draw("same");

		tlatFracLOFX[i] = new TLatex(0.15, 0.2, Form("Frac Within +/- 0.05 = %.3g", h_lofvtx_diff_y[i]->Integral(h_lofvtx_diff_y[i]->FindBin(-0.05), h_lofvtx_diff_y[i]->FindBin(0.05)) / h_lofvtx_diff_y[i]->Integral()));
		tlatFracLOFX[i]->SetTextColor(kBlue);
		tlatFracLOFX[i]->SetTextSize(0.035);
		tlatFracLOFX[i]->SetNDC(kTRUE);
		tlatFracLOFX[i]->Draw("same");

		tlatFracPRECX[i] = new TLatex(0.15, 0.17, Form("Frac Within +/- 0.05 = %.3g", h_precvtx_diff_y[i]->Integral(h_precvtx_diff_y[i]->FindBin(-0.05), h_precvtx_diff_y[i]->FindBin(0.05)) / h_precvtx_diff_y[i]->Integral()));
		tlatFracPRECX[i]->SetTextColor(kRed);
		tlatFracPRECX[i]->SetTextSize(0.035);
		tlatFracPRECX[i]->SetNDC(kTRUE);
		tlatFracPRECX[i]->Draw("same");

		cDiff[i]->cd(3);
		//gPad->SetLogy();

		h_lofvtx_diff_z[i]->GetXaxis()->SetTitle("(reco - pisa)_z [cm]");
		h_lofvtx_diff_z[i]->GetXaxis()->SetTitleFont(62);
		h_lofvtx_diff_z[i]->GetXaxis()->SetLabelFont(62);
		h_lofvtx_diff_z[i]->GetYaxis()->SetRangeUser(0, 1.2);
		h_lofvtx_diff_z[i]->GetYaxis()->SetTitle("AU");
		h_lofvtx_diff_z[i]->GetYaxis()->SetTitleFont(62);
		h_lofvtx_diff_z[i]->GetYaxis()->SetLabelFont(62);
		h_lofvtx_diff_z[i]->Draw();
		h_lofvtx_diff_z[i]->SetTitle("");
		h_precvtx_diff_z[i]->SetLineColor(kRed);
		h_precvtx_diff_z[i]->Draw("same");

		tlatMeanLOFZ[i] = new TLatex(0.15, 0.8, Form("Mean = %.3g", h_lofvtx_diff_z[i]->GetMean()));
		tlatMeanLOFZ[i]->SetTextColor(kBlue);
		tlatMeanLOFZ[i]->SetTextSize(0.035);
		tlatMeanLOFZ[i]->SetNDC(kTRUE);
		tlatMeanLOFZ[i]->Draw("same");

		tlatMeanPRECZ[i] = new TLatex(0.15, 0.77, Form("Mean = %.3g", h_precvtx_diff_z[i]->GetMean()));
		tlatMeanPRECZ[i]->SetTextColor(kRed);
		tlatMeanPRECZ[i]->SetTextSize(0.035);
		tlatMeanPRECZ[i]->SetNDC(kTRUE);
		tlatMeanPRECZ[i]->Draw("same");

		tlatRMSLOFZ[i] = new TLatex(0.15, 0.6, Form("RMS = %.3g", h_lofvtx_diff_z[i]->GetRMS()));
		tlatRMSLOFZ[i]->SetTextColor(kBlue);
		tlatRMSLOFZ[i]->SetTextSize(0.035);
		tlatRMSLOFZ[i]->SetNDC(kTRUE);
		tlatRMSLOFZ[i]->Draw("same");

		tlatRMSPRECZ[i] = new TLatex(0.15, 0.57, Form("RMS = %.3g", h_precvtx_diff_z[i]->GetRMS()));
		tlatRMSPRECZ[i]->SetTextColor(kRed);
		tlatRMSPRECZ[i]->SetTextSize(0.035);
		tlatRMSPRECZ[i]->SetNDC(kTRUE);
		tlatRMSPRECZ[i]->Draw("same");

		tlatFracLOFX[i] = new TLatex(0.15, 0.2, Form("Frac Within +/- 0.05 = %.3g", h_lofvtx_diff_z[i]->Integral(h_lofvtx_diff_z[i]->FindBin(-0.05), h_lofvtx_diff_z[i]->FindBin(0.05)) / h_lofvtx_diff_z[i]->Integral()));
		tlatFracLOFX[i]->SetTextColor(kBlue);
		tlatFracLOFX[i]->SetTextSize(0.035);
		tlatFracLOFX[i]->SetNDC(kTRUE);
		tlatFracLOFX[i]->Draw("same");

		tlatFracPRECX[i] = new TLatex(0.15, 0.17, Form("Frac Within +/- 0.05 = %.3g", h_precvtx_diff_z[i]->Integral(h_precvtx_diff_z[i]->FindBin(-0.05), h_precvtx_diff_z[i]->FindBin(0.05)) / h_precvtx_diff_z[i]->Integral()));
		tlatFracPRECX[i]->SetTextColor(kRed);
		tlatFracPRECX[i]->SetTextSize(0.035);
		tlatFracPRECX[i]->SetNDC(kTRUE);
		tlatFracPRECX[i]->Draw("same");
	}
}

void drawVertexDistributions()
{
	gStyle->SetOptStat(0);

	//Plot distribution of x,y,z vertices
	TCanvas *c_vertex_0 = new TCanvas("c_vertex_0", Form("Vertex %s", segmentCutLabel[0].c_str()), 1200, 500);
	c_vertex_0->Divide(3, 1);
	c_vertex_0->cd(1);
	gPad->SetLogy();

	h_lofvtx_x[0]->Scale(1.0 / h_lofvtx_x[0]->GetMaximum());
	h_precvtx_x[0]->Scale(1.0 / h_precvtx_x[0]->GetMaximum());
	h_pisavtx_x[0]->Scale(1.0 / h_pisavtx_x[0]->GetMaximum());

	h_lofvtx_x[0]->GetXaxis()->SetTitle("x-vertex [cm]");
	h_lofvtx_x[0]->GetXaxis()->SetTitleFont(62);
	h_lofvtx_x[0]->GetXaxis()->SetLabelFont(62);
	h_lofvtx_x[0]->GetYaxis()->SetTitle("AU");
	h_lofvtx_x[0]->GetYaxis()->SetTitleFont(62);
	h_lofvtx_x[0]->GetYaxis()->SetLabelFont(62);
	h_lofvtx_x[0]->Draw();
	h_lofvtx_x[0]->SetTitle("");
	h_precvtx_x[0]->SetLineColor(kRed);
	h_precvtx_x[0]->Draw("same");
	h_pisavtx_x[0]->SetLineColor(kBlack);
	h_pisavtx_x[0]->Draw("same");

	TLegend *leg_nseg0_vertex = new TLegend(0.15, 0.7, 0.4, 0.85);
	leg_nseg0_vertex->AddEntry(h_lofvtx_x[0], "LOF", "L");
	leg_nseg0_vertex->AddEntry(h_precvtx_x[0], "PRECISE", "L");
	leg_nseg0_vertex->AddEntry(h_pisavtx_x[0], "AMPT", "L");
	leg_nseg0_vertex->SetLineColor(kWhite);
	leg_nseg0_vertex->SetTextSize(0.05);
	leg_nseg0_vertex->Draw("same");

	c_vertex_0->cd(2);
	gPad->SetLogy();

	h_lofvtx_y[0]->Scale(1.0 / h_lofvtx_y[0]->GetMaximum());
	h_precvtx_y[0]->Scale(1.0 / h_precvtx_y[0]->GetMaximum());
	h_pisavtx_y[0]->Scale(1.0 / h_pisavtx_y[0]->GetMaximum());

	h_lofvtx_y[0]->GetXaxis()->SetTitle("y-vertex [cm]");
	h_lofvtx_y[0]->GetXaxis()->SetTitleFont(62);
	h_lofvtx_y[0]->GetXaxis()->SetLabelFont(62);
	h_lofvtx_y[0]->GetYaxis()->SetTitle("AU");
	h_lofvtx_y[0]->GetYaxis()->SetTitleFont(62);
	h_lofvtx_y[0]->GetYaxis()->SetLabelFont(62);
	h_lofvtx_y[0]->Draw();
	h_lofvtx_y[0]->SetTitle("");
	h_precvtx_y[0]->SetLineColor(kRed);
	h_precvtx_y[0]->Draw("same");
	h_pisavtx_y[0]->SetLineColor(kBlack);
	h_pisavtx_y[0]->Draw("same");

	c_vertex_0->cd(3);

	h_lofvtx_z[0]->Scale(1.0 / h_lofvtx_z[0]->GetMaximum());
	h_precvtx_z[0]->Scale(1.0 / h_precvtx_z[0]->GetMaximum());
	h_pisavtx_z[0]->Scale(1.0 / h_pisavtx_z[0]->GetMaximum());

	h_lofvtx_z[0]->GetXaxis()->SetTitle("z-vertex [cm]");
	h_lofvtx_z[0]->GetXaxis()->SetTitleFont(62);
	h_lofvtx_z[0]->GetXaxis()->SetLabelFont(62);
	h_lofvtx_z[0]->GetYaxis()->SetTitle("AU");
	h_lofvtx_z[0]->GetYaxis()->SetTitleFont(62);
	h_lofvtx_z[0]->GetYaxis()->SetLabelFont(62);
	h_lofvtx_z[0]->Draw();
	h_lofvtx_z[0]->SetTitle("");
	h_precvtx_z[0]->SetLineColor(kRed);
	h_precvtx_z[0]->Draw("same");
	h_pisavtx_z[0]->SetLineColor(kBlack);
	h_pisavtx_z[0]->Draw("same");


	TCanvas *c_vertex_1 = new TCanvas("c_vertex_1", Form("Vertex %s", segmentCutLabel[1].c_str()), 1200, 500);
	c_vertex_1->Divide(3, 1);
	c_vertex_1->cd(1);
	gPad->SetLogy();

	h_lofvtx_x[1]->Scale(1.0 / h_lofvtx_x[1]->Integral());

	h_lofvtx_x[1]->GetXaxis()->SetTitle("x-vertex [cm]");
	h_lofvtx_x[1]->GetXaxis()->SetTitleFont(62);
	h_lofvtx_x[1]->GetXaxis()->SetLabelFont(62);
	h_lofvtx_x[1]->GetYaxis()->SetTitle("Counts");
	h_lofvtx_x[1]->GetYaxis()->SetTitleFont(62);
	h_lofvtx_x[1]->GetYaxis()->SetLabelFont(62);
	h_lofvtx_x[1]->Draw();
	h_lofvtx_x[1]->SetTitle("");
	h_precvtx_x[1]->SetLineColor(kRed);
	h_precvtx_x[1]->Draw("same");
	h_pisavtx_x[1]->SetLineColor(kBlack);
	h_pisavtx_x[1]->Draw("same");

	TLegend *leg_nseg1_vertex = new TLegend(0.15, 0.7, 0.4, 0.8);
	leg_nseg1_vertex->AddEntry(h_lofvtx_x[1], "LOF", "L");
	leg_nseg1_vertex->AddEntry(h_precvtx_x[1], "PRECISE", "L");
	leg_nseg1_vertex->SetLineColor(kWhite);
	leg_nseg1_vertex->Draw("same");

	c_vertex_1->cd(2);
	gPad->SetLogy();

	h_lofvtx_y[1]->Scale(1.0 / h_lofvtx_y[1]->Integral());

	h_lofvtx_y[1]->GetXaxis()->SetTitle("y-vertex [cm]");
	h_lofvtx_y[1]->GetXaxis()->SetTitleFont(62);
	h_lofvtx_y[1]->GetXaxis()->SetLabelFont(62);
	h_lofvtx_y[1]->GetYaxis()->SetTitle("Counts");
	h_lofvtx_y[1]->GetYaxis()->SetTitleFont(62);
	h_lofvtx_y[1]->GetYaxis()->SetLabelFont(62);
	h_lofvtx_y[1]->Draw();
	h_lofvtx_y[1]->SetTitle("");
	h_precvtx_y[1]->SetLineColor(kRed);
	h_precvtx_y[1]->Draw("same");
	h_pisavtx_y[1]->SetLineColor(kBlack);
	h_pisavtx_y[1]->Draw("same");
	c_vertex_1->cd(3);

	h_lofvtx_z[1]->GetXaxis()->SetTitle("z-vertex [cm]");
	h_lofvtx_z[1]->GetXaxis()->SetTitleFont(62);
	h_lofvtx_z[1]->GetXaxis()->SetLabelFont(62);
	h_lofvtx_z[1]->GetYaxis()->SetTitle("Counts");
	h_lofvtx_z[1]->GetYaxis()->SetTitleFont(62);
	h_lofvtx_z[1]->GetYaxis()->SetLabelFont(62);
	h_lofvtx_z[1]->Draw();
	h_lofvtx_z[1]->SetTitle("");
	h_precvtx_z[1]->SetLineColor(kRed);
	h_precvtx_z[1]->Draw("same");
	h_pisavtx_z[1]->SetLineColor(kBlack);
	h_pisavtx_z[1]->Draw("same");

	TCanvas *c_vertex2 = new TCanvas("c_vertex2", Form("Vertex %s", segmentCutLabel[2].c_str()), 1200, 500);
	c_vertex2->Divide(3, 1);
	c_vertex2->cd(1);
	gPad->SetLogy();
	h_lofvtx_x[2]->GetXaxis()->SetTitle("x-vertex [cm]");
	h_lofvtx_x[2]->GetXaxis()->SetTitleFont(62);
	h_lofvtx_x[2]->GetXaxis()->SetLabelFont(62);
	h_lofvtx_x[2]->GetYaxis()->SetTitle("Counts");
	h_lofvtx_x[2]->GetYaxis()->SetTitleFont(62);
	h_lofvtx_x[2]->GetYaxis()->SetLabelFont(62);
	h_lofvtx_x[2]->Draw();
	h_lofvtx_x[2]->SetTitle("");
	h_precvtx_x[2]->SetLineColor(kRed);
	h_precvtx_x[2]->Draw("same");
	h_pisavtx_x[2]->SetLineColor(kBlack);
	h_pisavtx_x[2]->Draw("same");

	TLegend *leg_nseg2_vertex = new TLegend(0.15, 0.7, 0.4, 0.8);
	leg_nseg2_vertex->AddEntry(h_lofvtx_x[2], "LOF", "L");
	leg_nseg2_vertex->AddEntry(h_precvtx_x[2], "PRECISE", "L");
	leg_nseg2_vertex->SetLineColor(kWhite);
	leg_nseg2_vertex->Draw("same");

	c_vertex2->cd(2);
	gPad->SetLogy();
	h_lofvtx_y[2]->GetXaxis()->SetTitle("y-vertex [cm]");
	h_lofvtx_y[2]->GetXaxis()->SetTitleFont(62);
	h_lofvtx_y[2]->GetXaxis()->SetLabelFont(62);
	h_lofvtx_y[2]->GetYaxis()->SetTitle("Counts");
	h_lofvtx_y[2]->GetYaxis()->SetTitleFont(62);
	h_lofvtx_y[2]->GetYaxis()->SetLabelFont(62);
	h_lofvtx_y[2]->Draw();
	h_lofvtx_y[2]->SetTitle("");
	h_precvtx_y[2]->SetLineColor(kRed);
	h_precvtx_y[2]->Draw("same");
	h_pisavtx_y[2]->SetLineColor(kBlack);
	h_pisavtx_y[2]->Draw("same");

	c_vertex2->cd(3);
	h_lofvtx_z[2]->GetXaxis()->SetTitle("z-vertex [cm]");
	h_lofvtx_z[2]->GetXaxis()->SetTitleFont(62);
	h_lofvtx_z[2]->GetXaxis()->SetLabelFont(62);
	h_lofvtx_z[2]->GetYaxis()->SetTitle("Counts");
	h_lofvtx_z[2]->GetYaxis()->SetTitleFont(62);
	h_lofvtx_z[2]->GetYaxis()->SetLabelFont(62);
	h_lofvtx_z[2]->Draw();
	h_lofvtx_z[2]->SetTitle("");
	h_precvtx_z[2]->SetLineColor(kRed);
	h_precvtx_z[2]->Draw("same");
	h_pisavtx_z[2]->SetLineColor(kBlack);
	h_pisavtx_z[2]->Draw("same");

	TCanvas *c_vertex3 = new TCanvas("c_vertex3", Form("Vertex %s", segmentCutLabel[3].c_str()), 1200, 500);
	c_vertex3->Divide(3, 1);
	c_vertex3->cd(1);
	gPad->SetLogy();
	h_lofvtx_x[3]->GetXaxis()->SetTitle("x-vertex [cm]");
	h_lofvtx_x[3]->GetXaxis()->SetTitleFont(62);
	h_lofvtx_x[3]->GetXaxis()->SetLabelFont(62);
	h_lofvtx_x[3]->GetYaxis()->SetTitle("Counts");
	h_lofvtx_x[3]->GetYaxis()->SetTitleFont(62);
	h_lofvtx_x[3]->GetYaxis()->SetLabelFont(62);
	h_lofvtx_x[3]->Draw();
	h_lofvtx_x[3]->SetTitle("");
	h_precvtx_x[3]->SetLineColor(kRed);
	h_precvtx_x[3]->Draw("same");
	h_pisavtx_x[3]->SetLineColor(kBlack);
	h_pisavtx_x[3]->Draw("same");

	TLegend *leg_nseg3_vertex = new TLegend(0.15, 0.7, 0.4, 0.8);
	leg_nseg3_vertex->AddEntry(h_lofvtx_x[3], "LOF", "L");
	leg_nseg3_vertex->AddEntry(h_precvtx_x[3], "PRECISE", "L");
	leg_nseg3_vertex->SetLineColor(kWhite);
	leg_nseg3_vertex->Draw("same");

	c_vertex3->cd(2);
	gPad->SetLogy();
	h_lofvtx_y[3]->GetXaxis()->SetTitle("y-vertex [cm]");
	h_lofvtx_y[3]->GetXaxis()->SetTitleFont(62);
	h_lofvtx_y[3]->GetXaxis()->SetLabelFont(62);
	h_lofvtx_y[3]->GetYaxis()->SetTitle("Counts");
	h_lofvtx_y[3]->GetYaxis()->SetTitleFont(62);
	h_lofvtx_y[3]->GetYaxis()->SetLabelFont(62);
	h_lofvtx_y[3]->Draw();
	h_lofvtx_y[3]->SetTitle("");
	h_precvtx_y[3]->SetLineColor(kRed);
	h_precvtx_y[3]->Draw("same");
	h_pisavtx_y[3]->SetLineColor(kBlack);
	h_pisavtx_y[3]->Draw("same");

	c_vertex3->cd(3);
	h_lofvtx_z[3]->GetXaxis()->SetTitle("z-vertex [cm]");
	h_lofvtx_z[3]->GetXaxis()->SetTitleFont(62);
	h_lofvtx_z[3]->GetXaxis()->SetLabelFont(62);
	h_lofvtx_z[3]->GetYaxis()->SetTitle("Counts");
	h_lofvtx_z[3]->GetYaxis()->SetTitleFont(62);
	h_lofvtx_z[3]->GetYaxis()->SetLabelFont(62);
	h_lofvtx_z[3]->Draw();
	h_lofvtx_z[3]->SetTitle("");
	h_precvtx_z[3]->SetLineColor(kRed);
	h_precvtx_z[3]->Draw("same");
	h_pisavtx_z[3]->SetLineColor(kBlack);
	h_pisavtx_z[3]->Draw("same");

	TCanvas *c_vertex4 = new TCanvas("c_vertex4", Form("Vertex %s", segmentCutLabel[4].c_str()), 1200, 500);
	c_vertex4->Divide(3, 1);
	c_vertex4->cd(1);
	gPad->SetLogy();
	h_lofvtx_x[4]->GetXaxis()->SetTitle("x-vertex [cm]");
	h_lofvtx_x[4]->GetXaxis()->SetTitleFont(62);
	h_lofvtx_x[4]->GetXaxis()->SetLabelFont(62);
	h_lofvtx_x[4]->GetYaxis()->SetTitle("Counts");
	h_lofvtx_x[4]->GetYaxis()->SetTitleFont(62);
	h_lofvtx_x[4]->GetYaxis()->SetLabelFont(62);
	h_lofvtx_x[4]->Draw();
	h_lofvtx_x[4]->SetTitle("");
	h_precvtx_x[4]->SetLineColor(kRed);
	h_precvtx_x[4]->Draw("same");
	h_pisavtx_x[4]->SetLineColor(kBlack);
	h_pisavtx_x[4]->Draw("same");

	TLegend *leg_nseg4_vertex = new TLegend(0.15, 0.7, 0.4, 0.8);
	leg_nseg4_vertex->AddEntry(h_lofvtx_x[4], "LOF", "L");
	leg_nseg4_vertex->AddEntry(h_precvtx_x[4], "PRECISE", "L");
	leg_nseg4_vertex->SetLineColor(kWhite);
	leg_nseg4_vertex->Draw("same");

	c_vertex4->cd(2);
	gPad->SetLogy();
	h_lofvtx_y[4]->GetXaxis()->SetTitle("y-vertex [cm]");
	h_lofvtx_y[4]->GetXaxis()->SetTitleFont(62);
	h_lofvtx_y[4]->GetXaxis()->SetLabelFont(62);
	h_lofvtx_y[4]->GetYaxis()->SetTitle("Counts");
	h_lofvtx_y[4]->GetYaxis()->SetTitleFont(62);
	h_lofvtx_y[4]->GetYaxis()->SetLabelFont(62);
	h_lofvtx_y[4]->Draw();
	h_lofvtx_y[4]->SetTitle("");
	h_precvtx_y[4]->SetLineColor(kRed);
	h_precvtx_y[4]->Draw("same");
	h_pisavtx_y[4]->SetLineColor(kBlack);
	h_pisavtx_y[4]->Draw("same");

	c_vertex4->cd(3);
	h_lofvtx_z[4]->GetXaxis()->SetTitle("x-vertex [cm]");
	h_lofvtx_z[4]->GetXaxis()->SetTitleFont(62);
	h_lofvtx_z[4]->GetXaxis()->SetLabelFont(62);
	h_lofvtx_z[4]->GetYaxis()->SetTitle("Counts");
	h_lofvtx_z[4]->GetYaxis()->SetTitleFont(62);
	h_lofvtx_z[4]->GetYaxis()->SetLabelFont(62);
	h_lofvtx_z[4]->Draw();
	h_lofvtx_z[4]->SetTitle("");
	h_precvtx_z[4]->SetLineColor(kRed);
	h_precvtx_z[4]->Draw("same");
	h_pisavtx_z[4]->SetLineColor(kBlack);
	h_pisavtx_z[4]->Draw("same");

	TCanvas *c_vertex5 = new TCanvas("c_vertex5", Form("Vertex %s", segmentCutLabel[5].c_str()), 1200, 500);
	c_vertex5->Divide(3, 1);
	c_vertex5->cd(1);
	gPad->SetLogy();
	h_lofvtx_x[5]->GetXaxis()->SetTitle("x-vertex [cm]");
	h_lofvtx_x[5]->GetXaxis()->SetTitleFont(62);
	h_lofvtx_x[5]->GetXaxis()->SetLabelFont(62);
	h_lofvtx_x[5]->GetYaxis()->SetTitle("Counts");
	h_lofvtx_x[5]->GetYaxis()->SetTitleFont(62);
	h_lofvtx_x[5]->GetYaxis()->SetLabelFont(62);
	h_lofvtx_x[5]->Draw();
	h_lofvtx_x[5]->SetTitle("");
	h_precvtx_x[5]->SetLineColor(kRed);
	h_precvtx_x[5]->Draw("same");
	h_pisavtx_x[5]->SetLineColor(kBlack);
	h_pisavtx_x[5]->Draw("same");

	TLegend *leg_nseg5_vertex = new TLegend(0.15, 0.7, 0.4, 0.8);
	leg_nseg5_vertex->AddEntry(h_lofvtx_x[5], "LOF", "L");
	leg_nseg5_vertex->AddEntry(h_precvtx_x[5], "PRECISE", "L");
	leg_nseg5_vertex->SetLineColor(kWhite);
	leg_nseg5_vertex->Draw("same");

	c_vertex5->cd(2);
	gPad->SetLogy();
	h_lofvtx_y[5]->GetXaxis()->SetTitle("x-vertex [cm]");
	h_lofvtx_y[5]->GetXaxis()->SetTitleFont(62);
	h_lofvtx_y[5]->GetXaxis()->SetLabelFont(62);
	h_lofvtx_y[5]->GetYaxis()->SetTitle("Counts");
	h_lofvtx_y[5]->GetYaxis()->SetTitleFont(62);
	h_lofvtx_y[5]->GetYaxis()->SetLabelFont(62);
	h_lofvtx_y[5]->Draw();
	h_lofvtx_y[5]->SetTitle("");
	h_precvtx_y[5]->SetLineColor(kRed);
	h_precvtx_y[5]->Draw("same");
	h_pisavtx_y[5]->SetLineColor(kBlack);
	h_pisavtx_y[5]->Draw("same");

	c_vertex5->cd(3);
	h_lofvtx_z[5]->GetXaxis()->SetTitle("x-vertex [cm]");
	h_lofvtx_z[5]->GetXaxis()->SetTitleFont(62);
	h_lofvtx_z[5]->GetXaxis()->SetLabelFont(62);
	h_lofvtx_z[5]->GetYaxis()->SetTitle("Counts");
	h_lofvtx_z[5]->GetYaxis()->SetTitleFont(62);
	h_lofvtx_z[5]->GetYaxis()->SetLabelFont(62);
	h_lofvtx_z[5]->Draw();
	h_lofvtx_z[5]->SetTitle("");
	h_precvtx_z[5]->SetLineColor(kRed);
	h_precvtx_z[5]->Draw("same");
	h_pisavtx_z[5]->SetLineColor(kBlack);
	h_pisavtx_z[5]->Draw("same");
}

void getVertexDistributions()
{
	for (int i = 0; i < NCUTS; i++)
	{
		//PISA vertices
		ntp_event->Draw(Form("vtx_pisa[0]>>h_pisavtx_x_%i(200,-0.04,0.36)", i), segmentCut[i].c_str(), "goff");
		h_pisavtx_x[i] = (TH1F*) gDirectory->FindObject(Form("h_pisavtx_x_%i", i));

		ntp_event->Draw(Form("vtx_pisa[1]>>h_pisavtx_y_%i(200,-0.125,0.275)", i), segmentCut[i].c_str(), "goff");
		h_pisavtx_y[i] = (TH1F*) gDirectory->FindObject(Form("h_pisavtx_y_%i", i));

		ntp_event->Draw(Form("vtx_pisa[2]>>h_pisavtx_z_%i(400,-20,20)", i), segmentCut[i].c_str(), "goff");
		h_pisavtx_z[i] = (TH1F*) gDirectory->FindObject(Form("h_pisavtx_z_%i", i));

		//LOF vertices
		ntp_event->Draw(Form("vtx_lof[0]>>h_lofvtx_x_%i(200,-0.04,0.36)", i), segmentCut[i].c_str(), "goff");
		h_lofvtx_x[i] = (TH1F*) gDirectory->FindObject(Form("h_lofvtx_x_%i", i));

		ntp_event->Draw(Form("vtx_lof[1]>>h_lofvtx_y_%i(200,-0.125,0.275)", i), segmentCut[i].c_str(), "goff");
		h_lofvtx_y[i] = (TH1F*) gDirectory->FindObject(Form("h_lofvtx_y_%i", i));

		ntp_event->Draw(Form("vtx_lof[2]>>h_lofvtx_z_%i(400,-20,20)", i), segmentCut[i].c_str(), "goff");
		h_lofvtx_z[i] = (TH1F*) gDirectory->FindObject(Form("h_lofvtx_z_%i", i));

		//Precise vertices
		ntp_event->Draw(Form("vtx_prec[0]>>h_precvtx_x_%i(200,-0.5,0.5)", i), segmentCut[i].c_str(), "goff");
		h_precvtx_x[i] = (TH1F*) gDirectory->FindObject(Form("h_precvtx_x_%i", i));

		ntp_event->Draw(Form("vtx_prec[1]>>h_precvtx_y_%i(200,-0.125,0.275)", i), segmentCut[i].c_str(), "goff");
		h_precvtx_y[i] = (TH1F*) gDirectory->FindObject(Form("h_precvtx_y_%i", i));

		ntp_event->Draw(Form("vtx_prec[2]>>h_precvtx_z_%i(400,-20,20)", i), segmentCut[i].c_str(), "goff");
		h_precvtx_z[i] = (TH1F*) gDirectory->FindObject(Form("h_precvtx_z_%i", i));

		//Difference between LOF and PISA vertices
		ntp_event->Draw(Form("vtx_lof[0]-vtx_pisa[0]>>h_lofvtx_diff_x_%i(200,-0.2,0.2)", i), (segmentCut[i]).c_str(), "goff");
		h_lofvtx_diff_x[i] = (TH1F*) gDirectory->FindObject(Form("h_lofvtx_diff_x_%i", i));

		ntp_event->Draw(Form("vtx_lof[1]-vtx_pisa[1]>>h_lofvtx_diff_y_%i(200,-0.2,0.2)", i), (segmentCut[i]).c_str(), "goff");
		h_lofvtx_diff_y[i] = (TH1F*) gDirectory->FindObject(Form("h_lofvtx_diff_y_%i", i));

		ntp_event->Draw(Form("vtx_lof[2]-vtx_pisa[2]>>h_lofvtx_diff_z_%i(200,-0.2,0.2)", i), segmentCut[i].c_str(), "goff");
		h_lofvtx_diff_z[i] = (TH1F*) gDirectory->FindObject(Form("h_lofvtx_diff_z_%i", i));

		//Difference between precise and PISA vertices
		ntp_event->Draw(Form("vtx_prec[0]-vtx_pisa[0]>>h_precvtx_diff_x_%i(200,-0.2,0.2)", i), (segmentCut[i]).c_str(), "goff");
		h_precvtx_diff_x[i] = (TH1F*) gDirectory->FindObject(Form("h_precvtx_diff_x_%i", i));

		ntp_event->Draw(Form("vtx_prec[1]-vtx_pisa[1]>>h_precvtx_diff_y_%i(200,-0.2,0.2)", i), (segmentCut[i]).c_str(), "goff");
		h_precvtx_diff_y[i] = (TH1F*) gDirectory->FindObject(Form("h_precvtx_diff_y_%i", i));

		ntp_event->Draw(Form("vtx_prec[2]-vtx_pisa[2]>>h_precvtx_diff_z_%i(200,-0.2,0.2)", i), segmentCut[i].c_str(), "goff");
		h_precvtx_diff_z[i] = (TH1F*) gDirectory->FindObject(Form("h_precvtx_diff_z_%i", i));
	}
}

void getEventFractionNarrowVtx()
{
	for (int i = 0; i < NCUTS; i++)
	{
		cout << "****  " << segmentCutLabel[i] << "  ****" << endl;
		cout << "  Fraction of Events Within +/- 0.05 cm" << endl;
		cout << "       --->In Current Algorithm = " << h_precvtx_diff_z[i]->Integral(h_precvtx_diff_z[i]->FindBin(-0.05), h_precvtx_diff_z[i]->FindBin(0.05)) / h_precvtx_diff_z[i]->Integral() << endl;
		cout << "       --->In LOF     Algorithm = " << h_lofvtx_diff_z[i]->Integral(h_lofvtx_diff_z[i]->FindBin(-0.05), h_lofvtx_diff_z[i]->FindBin(0.05)) / h_lofvtx_diff_z[i]->Integral() << endl;

	}
}

void drawEventFractionWithVertex()
{
	float xvals[5] = {2, 3, 4, 5, 6};
	float fracLOF[5] = {0};
	float fracPREC[5] = {0};

	//Fraction of MB, narrow vtx events with nseg = xvals[i]
	fracLOF[0] = (float) ntp_event->GetEntries("TMath::Abs(vtx_bbc[2]) < 10 && pmtbbcn > 0 && pmtbbcs > 0 && vtx_lof[2] == vtx_lof[2] && nseg == 2") / ntp_event->GetEntries("TMath::Abs(vtx_bbc[2]) < 10 && pmtbbcn > 0 && pmtbbcs > 0 && nseg == 2");
	fracLOF[1] = (float) ntp_event->GetEntries("TMath::Abs(vtx_bbc[2]) < 10 && pmtbbcn > 0 && pmtbbcs > 0 && vtx_lof[2] == vtx_lof[2] && nseg == 3") / ntp_event->GetEntries("TMath::Abs(vtx_bbc[2]) < 10 && pmtbbcn > 0 && pmtbbcs > 0 && nseg == 3");
	fracLOF[2] = (float) ntp_event->GetEntries("TMath::Abs(vtx_bbc[2]) < 10 && pmtbbcn > 0 && pmtbbcs > 0 && vtx_lof[2] == vtx_lof[2] && nseg == 4") / ntp_event->GetEntries("TMath::Abs(vtx_bbc[2]) < 10 && pmtbbcn > 0 && pmtbbcs > 0 && nseg == 4");
	fracLOF[3] = (float) ntp_event->GetEntries("TMath::Abs(vtx_bbc[2]) < 10 && pmtbbcn > 0 && pmtbbcs > 0 && vtx_lof[2] == vtx_lof[2] && nseg == 5") / ntp_event->GetEntries("TMath::Abs(vtx_bbc[2]) < 10 && pmtbbcn > 0 && pmtbbcs > 0 && nseg == 5");
	fracLOF[4] = (float) ntp_event->GetEntries("TMath::Abs(vtx_bbc[2]) < 10 && pmtbbcn > 0 && pmtbbcs > 0 && vtx_lof[2] == vtx_lof[2] && nseg > 5") / ntp_event->GetEntries("TMath::Abs(vtx_bbc[2]) < 10 && pmtbbcn > 0 && pmtbbcs > 0 && nseg > 5");

	fracPREC[0] = (float) ntp_event->GetEntries("TMath::Abs(vtx_bbc[2]) < 10 && pmtbbcn > 0 && pmtbbcs > 0 && vtx_prec[2] == vtx_prec[2] && nseg == 2") / ntp_event->GetEntries("TMath::Abs(vtx_bbc[2]) < 10 && pmtbbcn > 0 && pmtbbcs > 0 && nseg == 2");
	fracPREC[1] = (float) ntp_event->GetEntries("TMath::Abs(vtx_bbc[2]) < 10 && pmtbbcn > 0 && pmtbbcs > 0 && vtx_prec[2] == vtx_prec[2] && nseg == 3") / ntp_event->GetEntries("TMath::Abs(vtx_bbc[2]) < 10 && pmtbbcn > 0 && pmtbbcs > 0 && nseg == 3");
	fracPREC[2] = (float) ntp_event->GetEntries("TMath::Abs(vtx_bbc[2]) < 10 && pmtbbcn > 0 && pmtbbcs > 0 && vtx_prec[2] == vtx_prec[2] && nseg == 4") / ntp_event->GetEntries("TMath::Abs(vtx_bbc[2]) < 10 && pmtbbcn > 0 && pmtbbcs > 0 && nseg == 4");
	fracPREC[3] = (float) ntp_event->GetEntries("TMath::Abs(vtx_bbc[2]) < 10 && pmtbbcn > 0 && pmtbbcs > 0 && vtx_prec[2] == vtx_prec[2] && nseg == 5") / ntp_event->GetEntries("TMath::Abs(vtx_bbc[2]) < 10 && pmtbbcn > 0 && pmtbbcs > 0 && nseg == 5");
	fracPREC[4] = (float) ntp_event->GetEntries("TMath::Abs(vtx_bbc[2]) < 10 && pmtbbcn > 0 && pmtbbcs > 0 && vtx_prec[2] == vtx_prec[2] && nseg > 5") / ntp_event->GetEntries("TMath::Abs(vtx_bbc[2]) < 10 && pmtbbcn > 0 && pmtbbcs > 0 && nseg > 5");

	TGraph *gEventFraction = new TGraph(5, xvals, fracLOF);
	TGraph *gEventFractionPrec = new TGraph(5, xvals, fracPREC);

	TCanvas *cEventFrac = new TCanvas("cEventFrac", "cEventFrac", 600, 600);

	gEventFraction->GetXaxis()->SetBinLabel(gEventFraction->GetXaxis()->FindBin(2), "2");
	gEventFraction->GetXaxis()->SetBinLabel(gEventFraction->GetXaxis()->FindBin(3), "3");
	gEventFraction->GetXaxis()->SetBinLabel(gEventFraction->GetXaxis()->FindBin(4), "4");
	gEventFraction->GetXaxis()->SetBinLabel(gEventFraction->GetXaxis()->FindBin(5), "5");
	gEventFraction->GetXaxis()->SetBinLabel(gEventFraction->GetXaxis()->FindBin(6), ">5");

	gEventFraction->GetXaxis()->SetTitleFont(62);
	gEventFraction->GetXaxis()->SetLabelFont(62);
	gEventFraction->GetYaxis()->SetTitleFont(62);
	gEventFraction->GetYaxis()->SetLabelFont(62);

	gEventFraction->GetXaxis()->SetLabelSize(0.045);

	gEventFraction->GetXaxis()->SetTitleOffset(1.3);
	gEventFraction->GetYaxis()->SetTitleOffset(1.3);

	gEventFraction->SetLineColor(kBlue);
	gEventFraction->SetMarkerColor(kBlue);

	gEventFraction->SetTitle("");
	gEventFraction->SetMarkerStyle(21);
	gEventFraction->GetXaxis()->SetTitle("Number of Segments");
	gEventFraction->GetYaxis()->SetTitle("Fraction of Events");
	gEventFraction->Draw("ALP");

	gEventFractionPrec->SetTitle("");
	gEventFractionPrec->SetMarkerStyle(21);
	gEventFractionPrec->SetLineColor(kRed);
	gEventFractionPrec->SetMarkerColor(kRed);
	gEventFractionPrec->Draw("LP,same");

	TLegend *tlEventFraction = new TLegend(0.6, 0.2, 0.85, 0.4);
	tlEventFraction->AddEntry(gEventFraction, "LOF", "LP");
	tlEventFraction->AddEntry(gEventFractionPrec, "PRECISE", "LP");
	tlEventFraction->SetLineColor(kWhite);
	tlEventFraction->SetTextSize(0.05);
	tlEventFraction->Draw("same");
}

void AnalyzeResolution()
{
	//Read file
	TFile *fin = new TFile("Data/423844_reco_lof.root");
	ntp_event  = (TTree*) fin->Get("ntp_event");

	getVertexDistributions();
	normalizeHistograms();
	fitResolutionHistograms();
	//drawVertexDistributions();
	drawVertexDifferenceDistributions();
	drawEventFractionWithVertex();
}