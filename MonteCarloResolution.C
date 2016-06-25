#include <iostream>

using namespace std;

//-------------------------------
// Variables
//-------------------------------

TH1F *h_total;
TH1F *h_resol;
TH1F *h_beam;

const float SIGMA_RES_X = 2.0;
const float SIGMA_RES_Y = 2.0;

const float SIGMA_BEAM_X = 3.0;
const float SIGMA_BEAM_Y = 3.0;

const float BEAM_CTR_X = 0.0;
const float BEAM_CTR_Y = 0.0;

const int NPOINTS = 100000;

//-------------------------------
// Functions
//-------------------------------

void MonteCarloResolution()
{
	TRandom rndm(0);

	for(int i=0; i<NPOINTS; i++)
	{
		//Generate true vertex by sampling around the beam center
		float vtx_true_x = rndm.Gaus(BEAM_CTR_X,SIGMA_BEAM_X);
		float vtx_true_y = rndm.Gaus(BEAM_CTR_Y,SIGMA_BEAM_Y);

		//Generate the reconstructed vertex by sampling around the collision point
		float vtx_reco_x = rndm.Gaus(vtx_true_x,SIGMA_RES_X);
		float vtx_reco_y = rndm.Gaus(vtx_true_y,SIGMA_RES_Y);
	}
}