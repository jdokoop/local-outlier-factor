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
	float reachDistance;
	vector<float> minPtsDistance;
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
const int NPOINTS = 200;

//Number of outliers to be generated with 4*sigma
const int NOUTLIERS = 10;

//Set of points to test the algorithm
vector<point> testPoints;

//------------------------------------------
// Functions
//------------------------------------------

/*
 * Compute euclidean distance between two points
 */
float computeDistance(point p1, point p2)
{
	return TMath::Sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
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

		for (int j = 0; j < testPoints.size(); j++)
		{
			if (i == j) continue;

			point p2 = testPoints[j];

			float dist = computeDistance(p1, p2);

			//Insert first two elements in order
			if(distances.size() == 0)
			{
				distances.push_back(dist);
			}
			else if(distances.size() == 1)
			{
				if(dist > distances[0]) 
				{
					distances.push_back(dist);
				}
				else
				{
					distances.insert(distances.begin(), dist);
				}
			}

			//Insert distance into array such that it is always ordered
			if (dist < distances[0])
			{
				distances.insert(distances.begin(), dist);
			}
			else if (dist > distances[distances.size() - 1])
			{
				distances.insert(distances.begin() + (distances.size()), dist);
			}
			else
			{
				for (int i = 0; i < distances.size() - 1; i++)
				{
					if (dist > distances[i] && dist < distances[i + 1])
					{
						distances.insert(distances.begin() + (i + 1), dist);
					}
				}
			}
		}

		//Having computed all distances, take the max distance to the k-nearest neighbors as the k-distance for the point
		testPoints[i].kDistance = distances[K-1];

		//Store the distance from a given point to its closest MINPTS neighbors
		for(int k=0; k<MINPTS; k++)
		{
			testPoints[i].minPtsDistance.push_back(distances[k]);
		}
	}
}

/*
 * Generate a synthetic test cluster with a Gaussian profile
 */
void generatePoints()
{
	TRandom *rand = new TRandom();
	TH2F *h = new TH2F("h", "h", 50, -8, 8, 50, -8, 8);

	//Generate points in cluster
	for (int i = 0; i < NPOINTS; i++)
	{
		point pt;
		pt.x = rand->Gaus(MEAN, SIGMA);
		pt.y = rand->Gaus(MEAN, SIGMA);

		testPoints.push_back(pt);
		h->Fill(pt.x, pt.y);
	}

	//Generate outliers sampled from Gaussian of width SIGMA_OUTLIERS > SIGMA
	for (int i = 0; i < NOUTLIERS; i++)
	{
		point pt;
		pt.x = rand->Gaus(MEAN, SIGMA_OUTLIERS);
		pt.y = rand->Gaus(MEAN, SIGMA_OUTLIERS);

		testPoints.push_back(pt);
		h->Fill(pt.x, pt.y);
	}

	TCanvas *c = new TCanvas("c", "c", 500, 500);
	gStyle->SetOptStat(0);
	h->Draw("LEGO");
}

void lof()
{
	generatePoints();
	findNearestNeighbors();
}