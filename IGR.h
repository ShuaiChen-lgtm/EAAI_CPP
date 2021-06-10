#pragma once
#include "LSMethods.h"
class IGR:
	public LSMethods
{
public:

	void SetParameter(int d, int Limp, float T, int time_factor);
	int IGR_Run(float lamda, float mu);  //IGA
	void IGR_Run2(float lamda, float mu, vector<int>& _results);
private:

	void ChooseExtJob(vector<int>& Des_Seq);

	int Destruction(const vector<int>& Des_Seq,
		vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij);
	int Reconstruction(const vector<int>& Des_Seq, 
		vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij);
	int Perturbation(vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij);

	void VNS(vector<vector<int>>& sol, vector<int>& factft, int& totaltft, vector<vector<vector<int>>>& perEij);

	void Reassignment(vector<vector<int>>& sol, vector<int>& factft, int& totaltft, vector<vector<vector<int>>>& perEij);

	void Permutation(vector<vector<int>>& sol, vector<int>& factft, int& totaltft, vector<vector<vector<int>>>& perEij);

	int m_d;
	int m_limp;
	float m_T;
	long m_CPUtime;

	vector < vector<vector<int>>> m_perEij;
};

