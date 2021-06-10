#pragma once
#include "ConstructiveHeuristics.h"
class HIG:
	public ConstructiveHeuristics
{
public:
	int HIG_Run(); //HIG1
	void SetParameters(float value_T, int Iters_perT, float Value_Cooling,
		float TLTlow, float TLThigh, int Minnum_Des, int Maxnum_Des, long time_fator);
private:

	void Destruction(vector<vector<int>>& Sol, vector<int>& Des_Seq, vector<int>& FacTFT, vector<int>& FacSpan,
		vector<vector<vector<int>>>& perEij);
	int Reconstruction(const vector<int>& Des_Seq, vector<vector<int>>& Sol, vector<int>& FacTFT, vector<int>& FacSpan,
		vector<vector<vector<int>>>& perEij);

	bool TestExsitJobnotinTL(vector<int> Seq);

	bool RemoveOneJob(vector<int>& Seq, vector<int>& Des_Seq, int TLT);

	void UpdateBestSol(const vector<vector<int>>& Sol, const vector<int>& FacTFT, int TotalTFT);
	void UpdateSol(const vector<vector<int>>& Sol, const vector<int>& FacTFT, const vector<int>& FacSpan, int TotalTFT,
		const vector<vector<vector<int>>>& perEij);

	queue<int> m_queueTL;
	vector<int> m_JobFlag;

	float m_Value_T;
	int m_Iters_perT;
	float m_Value_Cooling;
	float m_Temperature;

	float m_TLTlow;
	float m_TLThigh;
	int m_Minnum_Des;
	int m_Maxnum_Des;

	long m_limittime;

	vector<vector<int>> m_CurSol;
	vector<int> m_CurFacTFT;
	vector<int> m_CurFacSpan;
	int m_CurTotalTFT;
	vector<vector<vector<int>>> m_CurperEij;

};

