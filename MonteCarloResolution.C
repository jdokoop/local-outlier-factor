#include <iostream>

using namespace std;

//-------------------------------
// Variables
//-------------------------------

TH1F *h_total;
TH1F *h_resol;
TH1F *h_beam;

//Vertexing resolution
const float SIGMA_RES_X = 1.2 * 129.5;
const float SIGMA_RES_Y = 1.2 * 107.4;

//Beam spot spread
const float SIGMA_BEAM_X = 129.5;
const float SIGMA_BEAM_Y = 107.4;

//Beam spot center
const float BEAM_CTR_X = 1612.1;
const float BEAM_CTR_Y = 722.5;

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

//DCA distribution wrt the reconstructed vertex
TH1F *hDCAReco;

//DCA distribution wrt the beam center
TH1F *hDCABC;

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

	//Plot summary as a function of the number of tracks in each arm
	/*
	float x[4] = {1, 2, 3, 4};
	float y_reco_true[4] = {0.62, 0.72, 0.82, 0.93};
	float y_bc_true[4] = {0.77, 0.73, 0.79, 0.79};

	TGraph *g_reco_true = new TGraph(4, x, y_reco_true);
	TGraph *g_bc_true = new TGraph(4, x, y_bc_true);

	TCanvas *cSummary = new TCanvas("cSummary","cSummary",600,600);

	g_reco_true->SetTitle("");
	g_reco_true->SetMarkerStyle(20);
	g_reco_true->SetMarkerColor(kBlue);
	g_reco_true->SetLineColor(kBlue);
	g_reco_true->GetXaxis()->SetTitle("Number of Vertex Tracks in Each Arm");
	g_reco_true->GetXaxis()->SetTitleFont(62);
	g_reco_true->GetXaxis()->SetLabelFont(62);
	g_reco_true->GetYaxis()->SetTitle("Area Below 200");
	g_reco_true->GetYaxis()->SetTitleOffset(1.85);
	//g_reco_true->GetYaxis()->SetRangeUser(0, 0.025);
	g_reco_true->GetYaxis()->SetTitleFont(62);
	g_reco_true->GetYaxis()->SetLabelFont(62);
	g_reco_true->Draw("AP");
	g_bc_true->SetTitle("");
	g_bc_true->SetMarkerStyle(20);
	g_bc_true->SetMarkerColor(kRed);
	g_bc_true->SetLineColor(kRed);
	g_bc_true->GetXaxis()->SetTitleFont(62);
	g_bc_true->GetXaxis()->SetLabelFont(62);
	//g_bc_true->GetYaxis()->SetRangeUser(0, 0.025);
	g_bc_true->GetYaxis()->SetTitleFont(62);
	g_bc_true->GetYaxis()->SetLabelFont(62);
	g_bc_true->Draw("P,same");
	*/

	TCanvas *cResiduals = new TCanvas("cResiduals", "cResiduals", 600, 400);
	hDCAReco->Scale(1.0 / hDCAReco->Integral());
	hDCABC->Scale(1.0 / hDCABC->Integral());
	cResiduals->SetLeftMargin(0.2);
	hDCAReco->SetTitle("");
	hDCAReco->GetXaxis()->SetTitle("|DCA| [#mum]");
	hDCAReco->GetXaxis()->SetTitleFont(62);
	hDCAReco->GetXaxis()->SetLabelFont(62);
	hDCAReco->GetYaxis()->SetTitle("Probability");
	hDCAReco->GetYaxis()->SetTitleOffset(1.8);
	hDCAReco->GetYaxis()->SetTitleFont(62);
	hDCAReco->GetYaxis()->SetLabelFont(62);

	float max = (hDCAReco->GetMaximum() > hDCABC->GetMaximum()) ? hDCAReco->GetMaximum() : hDCABC->GetMaximum();
	hDCAReco->GetYaxis()->SetRangeUser(0, max + 0.1 * max);
	hDCAReco->SetLineColor(kBlue);
	hDCAReco->Draw();
	hDCABC->SetLineColor(kRed);
	hDCABC->Draw("same");

	TLegend *tlResiduals = new TLegend(0.6, 0.7, 0.85, 0.8);
	tlResiduals->AddEntry(hDCAReco, "RECO", "L");
	tlResiduals->AddEntry(hDCABC, "BC", "L");
	tlResiduals->SetLineColor(kWhite);
	tlResiduals->SetTextSize(0.05);
	tlResiduals->Draw("same");
}

/*
 * Get the distance of closest approach between the line defined by angle theta and point (ax, ay) and point (px, py)
 */
float getDCA(float theta, float px, float py, float ax, float ay)
{
	//Line parameter defining the point of closest approach to (x,y)
	float t = ((px * TMath::Cos(theta) + py * TMath::Sin(theta)) - (ax * TMath::Cos(theta) + ay * TMath::Sin(theta))) / (TMath::Cos(theta) * TMath::Cos(theta) + TMath::Sin(theta) * TMath::Sin(theta));

	//Point of closest approach
	float dcax = ax + t * TMath::Cos(theta);
	float dcay = ay + t * TMath::Sin(theta);

	//Check that the lines are perpendicular
	float testNorm = TMath::Sqrt((px - dcax) * (px - dcax) + (py - dcay) * (py - dcay));
	float ang = (1.0 / testNorm) * ((px - dcax) * TMath::Cos(theta) + (py - dcay) * TMath::Sin(theta));

	//Distance between (dcax, dcay) and (px,py)
	float dca = TMath::Sqrt((dcax - px) * (dcax - px) + (dcay - py) * (dcay - py));

	return dca;
}

/*
 * Run N simulated events by sampling a true and a reconstructed vertex, with a straight track from the true vertex
 * Compute the distance of closest approach between the track and the beam center, and the track and the reconstructed vertex
 */
void MonteCarloResolution()
{
	//Initialize variables
	hResTrueReco       = new TH1F("hResTrueReco", "hResTrueReco", 200, 0, 700);
	hResTrueBeamCenter = new TH1F("hResTrueBeamCenter", "hResTrueBeamCenter", 200, 0, 700);
	hRecoVertexX       = new TH1F("hRecoVertexX", "hRecoVertexX", 150, BEAM_CTR_X - 5 * SIGMA_BEAM_X, BEAM_CTR_X + 5 * SIGMA_BEAM_X);
	hRecoVertexY       = new TH1F("hRecoVertexY", "hRecoVertexY", 150, BEAM_CTR_Y - 5 * SIGMA_BEAM_Y, BEAM_CTR_Y + 5 * SIGMA_BEAM_Y);
	hTrueVertexX       = new TH1F("hTrueVertexX", "hTrueVertexX", 150, BEAM_CTR_X - 5 * SIGMA_BEAM_X, BEAM_CTR_X + 5 * SIGMA_BEAM_X);
	hTrueVertexY       = new TH1F("hTrueVertexY", "hTrueVertexY", 150, BEAM_CTR_Y - 5 * SIGMA_BEAM_Y, BEAM_CTR_Y + 5 * SIGMA_BEAM_Y);
	hDCAReco           = new TH1F("hDCAReco", "hDCAReco", 500, 0, 600);
	hDCABC             = new TH1F("hDCABC", "hDCABC", 500, 0, 600);

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

		//Generate straight track originating from the true vertex, parametrized by polar angle
		float trackAngle = rndm.Uniform(2 * TMath::Pi());
		float dcaRecoVtx = getDCA(trackAngle, vtx_reco_x, vtx_reco_y, vtx_true_x, vtx_true_y);
		float dcaBC      = getDCA(trackAngle, BEAM_CTR_X, BEAM_CTR_Y, vtx_true_x, vtx_true_y);

		hDCAReco->Fill(dcaRecoVtx);
		hDCABC->Fill(dcaBC);

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