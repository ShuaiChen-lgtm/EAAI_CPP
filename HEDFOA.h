#pragma once
#include "LSMethods.h"
class HEDFOA:
	public LSMethods
{
public:
	void SetParameter_HEDFOA(int SetNP, int SetSN, int Setd, float SetT, int SetC);
	int HEDFOA_Run();
	void HEDFOA_Run(vector<int>& Results);
private:

	int Smell_based_foraging(vector<vector<int>>& Solution, vector<int>& FacSpan, vector<int> &FacTFT,
		vector<vector<vector<int>>>& perEij);

	int m_NP;
	int m_SN;
	int m_d;
	int m_C;
	float m_T;

	vector<vector<vector<int>>> m_SwarmSol;
	vector<vector<int>> m_SwarmFacSpan;
	vector<vector<int>> m_SwarmFacTft;
	vector<int> m_SwarmTotaltft;

	vector<int> m_BestSol_FacSpan;

	long m_CPUtime;

};

