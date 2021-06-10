#pragma once
#include <vector>
using namespace std;

class Problem
{
public:
	Problem();
	~Problem();
	void SetInstance(int factories, int jobs, int machines, const vector<vector<int>>& pTime);
protected:

	void GetTotalPTime();
	void GetTotalLoad();

	int m_Factories;
	int m_Machines;
	int m_Jobs;

	vector <vector<int>> m_Ptime;
	vector<int> m_TotalPTime;
	int m_TotalLoad;

	vector<int> m_FullSeq;

	vector<vector<int>> m_BestSol_Seq;
	vector<int> m_BestSol_FacTFT;
	int m_BestSol_TotalTFT;
};

