#include <iostream>

using namespace std;

//-------------------------------
// Variables
//-------------------------------

TH1F *h_total;
TH1F *h_resol;
TH1F *h_beam;

//Vertexing resolution
const float SIGMA_RES_X = 201.6;
const float SIGMA_RES_Y = 156.9;

//Beam spot spread
const float SIGMA_BEAM_X = 131.8;
const float SIGMA_BEAM_Y = 114.6;

//Beam spot spread
const float BEAM_CTR_X = 1612.0;
const float BEAM_CTR_Y = 723.0;

//Number of events to simulate
const int NPOINTS = 1e6;

//Distribution of true vertices
TH1F *hTrueVertexX;
TH1F *hTrueVertexY;

//Distribution of reconstructed vertices
TH1F *hRecoVertexX;
TH1F *hRecoVertexY;

//Distribution of distances between the true and the reconstructed vertex
TH1F *hResTrueReco;

//Distribution of distances between the true vertex and the beam center
TH1F *hResTrueBeamCenter;

//-------------------------------
// Functions
//-------------------------------

/*
 * Plot relevant things
 */
void plot()
{
	gStyle->SetOptStat(0);

	TCanvas *cRecoVertex = new TCanvas("cRecoVertex", "cRecoVertex", 900, 400);
	cRecoVertex->Divide(2, 1);

	cRecoVertex->cd(1);
	hRecoVertexX->Scale(1.0 / hRecoVertexX->GetMaximum());
	hTrueVertexX->Scale(1.0 / hTrueVertexX->GetMaximum());
	hRecoVertexX->SetTitle("");
	hRecoVertexX->GetXaxis()->SetTitle("x-vertex [#mum]");
	hRecoVertexX->GetXaxis()->SetTitleFont(62);
	hRecoVertexX->GetXaxis()->SetLabelFont(62);
	hRecoVertexX->GetYaxis()->SetTitle("AU [max. norm.]");
	hRecoVertexX->GetYaxis()->SetTitleFont(62);
	hRecoVertexX->GetYaxis()->SetLabelFont(62);
	hTrueVertexX->SetLineColor(kRed);
	hRecoVertexX->Draw();
	hTrueVertexX->Draw("same");

	TLegend *tlVertices = new TLegend(0.6, 0.7, 0.85, 0.8);
	tlVertices->AddEntry(hRecoVertexX, "RECO", "L");
	tlVertices->AddEntry(hTrueVertexX, "TRUE", "L");
	tlVertices->SetLineColor(kWhite);
	tlVertices->SetTextSize(0.05);
	tlVertices->Draw("same");

	cRecoVertex->cd(2);
	hRecoVertexY->Scale(1.0 / hRecoVertexY->GetMaximum());
	hTrueVertexY->Scale(1.0 / hTrueVertexY->GetMaximum());
	hRecoVertexY->SetTitle("");
	hRecoVertexY->GetXaxis()->SetTitle("y-vertex [#mum]");
	hRecoVertexY->GetXaxis()->SetTitleFont(62);
	hRecoVertexY->GetXaxis()->SetLabelFont(62);
	hRecoVertexY->GetYaxis()->SetTitle("AU [max. norm.]");
	hRecoVertexY->GetYaxis()->SetTitleFont(62);
	hRecoVertexY->GetYaxis()->SetLabelFont(62);
	hTrueVertexY->SetLineColor(kRed);
	hRecoVertexY->Draw();
	hTrueVertexY->Draw("same");

	TCanvas *cResiduals = new TCanvas("cResiduals", "cResiduals", 600, 400);
	hResTrueReco->Scale(1.0/hResTrueReco->Integral());
	hResTrueBeamCenter->Scale(1.0/hResTrueBeamCenter->Integral());
	cResiduals->SetLeftMargin(0.2);
	hResTrueReco->SetTitle("");
	hResTrueReco->GetXaxis()->SetTitle("Distance [#mum]");
	hResTrueReco->GetXaxis()->SetTitleFont(62);
	hResTrueReco->GetXaxis()->SetLabelFont(62);
	hResTrueReco->GetYaxis()->SetTitle("Probability");
	hResTrueReco->GetYaxis()->SetTitleOffset(1.8);
	hResTrueReco->GetYaxis()->SetTitleFont(62);
	hResTrueReco->GetYaxis()->SetLabelFont(62);
	hResTrueReco->GetYaxis()->SetRangeUser(0, hResTrueBeamCenter->GetMaximum() + 0.1*hResTrueBeamCenter->GetMaximum());
	hResTrueReco->SetLineColor(kBlue);
	hResTrueReco->Draw();
	hResTrueBeamCenter->SetLineColor(kRed);
	hResTrueBeamCenter->Draw("same");

	TLegend *tlResiduals = new TLegend(0.6, 0.7, 0.85, 0.8);
	tlResiduals->AddEntry(hResTrueReco, "RECO - TRUE", "L");
	tlResiduals->AddEntry(hResTrueBeamCenter, "BC - TRUE", "L");
	tlResiduals->SetLineColor(kWhite);
	tlResiduals->SetTextSize(0.05);
	tlResiduals->Draw("same");

	float areaBelow200BC = hResTrueBeamCenter->Integral(1, hResTrueBeamCenter->FindBin(200.0)) / hResTrueBeamCenter->Integral();
	float areaBelow200Reco = hResTrueReco->Integral(1, hResTrueReco->FindBin(200.0)) / hResTrueReco->Integral();

	TLatex *tlatAreaBC = new TLatex(0.6, 0.4, Form("Area Below 200 = %.3g", areaBelow200BC));
	tlatAreaBC->SetTextColor(kRed);
	tlatAreaBC->SetTextSize(0.035);
	tlatAreaBC->SetNDC(kTRUE);
	tlatAreaBC->Draw("same");

	TLatex *tlatAreaReco = new TLatex(0.6, 0.35, Form("Area Below 200 = %.3g", areaBelow200Reco));
	tlatAreaReco->SetTextColor(kBlue);
	tlatAreaReco->SetTextSize(0.035);
	tlatAreaReco->SetNDC(kTRUE);
	tlatAreaReco->Draw("same");

	TLine *tline200 = new TLine(200, 0, 200, hResTrueBeamCenter->GetMaximum() + 0.1*hResTrueBeamCenter->GetMaximum());
	tline200->SetLineStyle(2);
	tline200->Draw("same");
}

void MonteCarloResolution()
{
	//Initialize variables
	hResTrueReco       = new TH1F("hResTrueReco", "hResTrueReco", 200, 0, 700);
	hResTrueBeamCenter = new TH1F("hResTrueBeamCenter", "hResTrueBeamCenter", 200, 0, 700);
	hRecoVertexX       = new TH1F("hRecoVertexX", "hRecoVertexX", 150, BEAM_CTR_X - 5 * SIGMA_BEAM_X, BEAM_CTR_X + 5 * SIGMA_BEAM_X);
	hRecoVertexY       = new TH1F("hRecoVertexY", "hRecoVertexY", 150, BEAM_CTR_Y - 5 * SIGMA_BEAM_Y, BEAM_CTR_Y + 5 * SIGMA_BEAM_Y);
	hTrueVertexX       = new TH1F("hTrueVertexX", "hTrueVertexX", 150, BEAM_CTR_X - 5 * SIGMA_BEAM_X, BEAM_CTR_X + 5 * SIGMA_BEAM_X);
	hTrueVertexY       = new TH1F("hTrueVertexY", "hTrueVertexY", 150, BEAM_CTR_Y - 5 * SIGMA_BEAM_Y, BEAM_CTR_Y + 5 * SIGMA_BEAM_Y);

	TRandom rndm(0);

	//Generate events
	for (int i = 0; i < NPOINTS; i++)
	{
		//Generate true vertex by sampling around the beam center
		float vtx_true_x = rndm.Gaus(BEAM_CTR_X, SIGMA_BEAM_X);
		float vtx_true_y = rndm.Gaus(BEAM_CTR_Y, SIGMA_BEAM_Y);

		//Generate the reconstructed vertex by sampling around the collision point
		float vtx_reco_x = rndm.Gaus(vtx_true_x, SIGMA_RES_X);
		float vtx_reco_y = rndm.Gaus(vtx_true_y, SIGMA_RES_Y);

		hRecoVertexX->Fill(vtx_reco_x);
		hRecoVertexY->Fill(vtx_reco_y);

		hTrueVertexX->Fill(vtx_true_x);
		hTrueVertexY->Fill(vtx_true_y);

		//Determine the difference between the reconstructed and true vertices
		float resTrueReco = TMath::Sqrt(pow(vtx_true_x - vtx_reco_x, 2) + pow(vtx_true_y - vtx_reco_y, 2));
		hResTrueReco->Fill(resTrueReco);

		//Determine the difference between the true vertex and the beam center
		float resTrueBC = TMath::Sqrt(pow(vtx_true_x - BEAM_CTR_X, 2) + pow(vtx_true_y - BEAM_CTR_Y, 2));
		hResTrueBeamCenter->Fill(resTrueBC);
	}

	plot();
}