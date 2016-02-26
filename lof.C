//--------------------------------------------------
// Test implementation of the local outlier factor
// algorithm for outlier detection, as published by
// Breuning et al.
//--------------------------------------------------

#include <iostream>
#include <vector>

#include "TRandom.h"
#include "TMath.h"

using namespace std;

//------------------------------------------
// Variables
//------------------------------------------

//Structure to represent a 2-dimensional point in Cartesian coordinates
//Each point
struct point
{
	float x;
	float y;

	float kDistance;
	float lrd;
	float lof;
	vector<int> minPtsNeighbors;
};

//The 'k' in k-nearest-neighbors
const int K = 5;

//The number of points to compute reachability density
const int MINPTS = 10;

//Mean and width of 2D Gaussian from which points are sampled
const float MEAN = 0.0;
const float SIGMA = 1.0;
const float SIGMA_OUTLIERS = 4.0;

//Number of points to be generated in the cluster
const int NPOINTS = 80;

//Number of outliers to be generated with 4*sigma
const int NOUTLIERS = 10;

//Set of points to test the algorithm
vector<point> testPoints;

//------------------------------------------
// Functions
//------------------------------------------

/*
 * Plot the LOF score of each point as a 2D histogram
 */
void plotLOF()
{
	TProfile2D *hLOF2D = new TProfile2D("hLOF2D", "hLOF2D", 50, -8, 8, 50, -8, 8);
	TH2F *hPoints = new TH2F("hPoints", "hPoints", 50, -8, 8, 50, -8, 8);
	TH1F *hLOF1D = new TH1F("hLOF1D", "hLOF1D", 100, 0, 5);

	for (int i = 0; i < testPoints.size(); i++)
	{
		point p = testPoints[i];
		float lof = p.lof;
		float x = p.x;
		float y = p.y;

		hLOF2D->Fill(x, y, lof);
		hPoints->Fill(x, y);
		hLOF1D->Fill(lof);
	}

	TCanvas *c1 = new TCanvas("c1", "c1", 500, 500);
	hLOF1D->SetTitle("LOF Score Distribution");
	hLOF1D->GetXaxis()->SetTitle("LOF");
	hLOF1D->Draw();

	TCanvas *c2 = new TCanvas("c2", "c2", 1000, 500);
	c2->Divide(2, 1);

	c2->cd(1);
	hPoints->SetTitle("Spatial Distribution of Points");
	hPoints->GetXaxis()->SetTitle("x");
	hPoints->GetXaxis()->SetTitleOffset(1.6);
	hPoints->GetYaxis()->SetTitle("y");
	hPoints->GetYaxis()->SetTitleOffset(1.6);
	hPoints->GetZaxis()->SetTitle("N_{points}");
	hPoints->GetZaxis()->SetTitleOffset(1.4);
	hPoints->Draw("LEGO20");

	c2->cd(2);
	hLOF2D->SetTitle("Spatial Distribution of LOF Scores");
	hLOF2D->GetXaxis()->SetTitle("x");
	hLOF2D->GetXaxis()->SetTitleOffset(1.6);
	hLOF2D->GetYaxis()->SetTitle("y");
	hLOF2D->GetYaxis()->SetTitleOffset(1.6);
	hLOF2D->GetZaxis()->SetTitle("< LOF >");
	hLOF2D->GetZaxis()->SetTitleOffset(1.4);
	hLOF2D->Draw("LEGO20");
}

/*
 * Compute euclidean distance between two points
 */
float computeDistance(point p1, point p2)
{
	return TMath::Sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

/*
 * Compute the reachability density for every point by averaging the reach distance over MINPTS neighbors
 */
void computeReachDensity()
{
	for (int i = 0; i < testPoints.size(); i++)
	{
		point p1 = testPoints[i];
		vector<int> neighbors = p1.minPtsNeighbors;
		float summedReachDist = 0;

		//Iterate over the MINPTS neighbors to ith point
		for (int j = 0; j < neighbors.size(); j++)
		{
			point p2 = testPoints[j];

			float kDist      = p2.kDistance;
			float euclidDist = computeDistance(p1, p2);
			float reachDist  = max(kDist, euclidDist);

			summedReachDist += reachDist;
		}

		testPoints[i].lrd = (float) neighbors.size() / summedReachDist;
	}
}

/*
 * Compute the LOF score for each point
 */
void computeLOF()
{
	for (int i = 0; i < testPoints.size(); i++)
	{
		point p1 = testPoints[i];
		vector<int> neighbors = p1.minPtsNeighbors;
		float lrd1 = p1.lrd;
		float summedLRDRatio = 0;

		for (int j = 0; j < neighbors.size(); j++)
		{
			point p2 = testPoints[j];
			float lrd2 = p2.lrd;
			summedLRDRatio += (lrd2 / lrd1);
		}

		testPoints[i].lof = (float) summedLRDRatio / neighbors.size();
	}
}

/*
 * Find the k-nearest neighbors to each point
 * Use brute force at first (small number of points)
 * Eventually implement something better...
 */
void findNearestNeighbors()
{
	vector<point> neighbors;

	//Iterate over points
	//Use brute force to compute the distance to every other point
	for (int i = 0; i < testPoints.size(); i++)
	{
		point p1 = testPoints[i];

		//Vector with distance from current point to every other point in the array
		vector<float> distances;
		//Vector with the indices of the points arranged in increasing distance to the current point
		vector<int> indices;

		for (int j = 0; j < testPoints.size(); j++)
		{
			if (i == j) continue;

			point p2 = testPoints[j];

			float dist = computeDistance(p1, p2);

			//Insert first two elements in order
			if (distances.size() == 0)
			{
				distances.push_back(dist);
				indices.push_back(j);
			}
			else if (distances.size() == 1)
			{
				if (dist > distances[0])
				{
					distances.push_back(dist);
					indices.push_back(j);
				}
				else
				{
					distances.insert(distances.begin(), dist);
					indices.insert(indices.begin(), j);
				}
			}

			//Insert distance into array such that it is always ordered
			if (dist < distances[0])
			{
				distances.insert(distances.begin(), dist);
				indices.insert(indices.begin(), j);
			}
			else if (dist > distances[distances.size() - 1])
			{
				distances.insert(distances.begin() + (distances.size()), dist);
				indices.insert(indices.begin() + (indices.size()), j);
			}
			else
			{
				for (int i = 0; i < distances.size() - 1; i++)
				{
					if (dist > distances[i] && dist < distances[i + 1])
					{
						distances.insert(distances.begin() + (i + 1), dist);
						indices.insert(indices.begin() + (i + 1), j);
					}
				}
			}
		}

		//Having computed all distances, take the max distance to the k-nearest neighbors as the k-distance for the point
		testPoints[i].kDistance = distances[K - 1];

		//Store the indices of the closest MINPTS neighbors in the testPoints vector
		for (int k = 0; k < MINPTS; k++)
		{
			testPoints[i].minPtsNeighbors.push_back(indices[k]);
		}
	}
}

/*
 * Generate a synthetic test cluster with a Gaussian profile
 */
void generateGaussianPoints()
{
	TRandom *rand = new TRandom();

	//Generate points in cluster
	for (int i = 0; i < NPOINTS; i++)
	{
		point pt;
		pt.x = rand->Gaus(MEAN, SIGMA);
		pt.y = rand->Gaus(MEAN, SIGMA);

		testPoints.push_back(pt);
	}

	//Generate outliers sampled from Gaussian of width SIGMA_OUTLIERS > SIGMA
	for (int i = 0; i < NOUTLIERS; i++)
	{
		point pt;
		pt.x = rand->Gaus(MEAN, SIGMA_OUTLIERS);
		pt.y = rand->Gaus(MEAN, SIGMA_OUTLIERS);

		testPoints.push_back(pt);
	}
}

/*
 * Generate a synthetic test cluster sampled uniformly in a disk
 */
void generateUniformPoints()
{
	TRandom *rand = new TRandom(0);

	//Generate points in cluster
	for (int i = 0; i < NPOINTS; i++)
	{
		float r = rand->Uniform(0.0,3);
		float phi = rand->Uniform(0,2*TMath::Pi());

		point pt;
		pt.x = r*TMath::Cos(phi);
		pt.y = r*TMath::Sin(phi);

		testPoints.push_back(pt);
	}

	//Generate outliers sampled from Gaussian of width SIGMA_OUTLIERS > SIGMA
	for (int i = 0; i < NOUTLIERS; i++)
	{
		float r = rand->Uniform(4,8);
		float phi = rand->Uniform(0,2*TMath::Pi());

		point pt;
		pt.x = r*TMath::Cos(phi);
		pt.y = r*TMath::Sin(phi);

		testPoints.push_back(pt);
	}
}

void lof()
{
	gStyle->SetOptStat(0);
	generateUniformPoints();
	findNearestNeighbors();
	computeReachDensity();
	computeLOF();
	plotLOF();
}