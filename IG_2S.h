#pragma once
#include "LSMethods.h"
class IG_2S:
	public LSMethods
{
public:
	void SetParameter_IG2S(int S_d1, int S_d2, double S_T, double S_ratio, long time_factor);
	int IG2S_Run(); //IG2S
private:
	void Destruction(vector<int>& des_seq,
		vector<vector<int>>& sol, vector<int>& facspan, vector<int>& factft, vector<vector<vector<int>>>& perEij);

	int Reconstruction(vector<int>& des_seq,
		vector<vector<int>>& sol, vector<int>& facspan, vector<int>& factft, vector<vector<vector<int>>>& perEij);

	void Destruction_2S(vector<int>& des_seq, vector<int>& seq);

	void Reconstruction_2S(vector<int>& des_seq, vector<int>& seq, vector<vector<int>>& Eij, int& mspan, int& tft);

	int m_d1;
	int m_d2;
	int m_T;
	double m_ratio;
	long m_CPU_time;

	vector<int> m_BestSol_FacSpan;

	vector<vector<vector<int>>> m_perEij;
};

