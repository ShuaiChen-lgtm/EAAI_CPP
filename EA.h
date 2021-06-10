#pragma once
#include "LSMethods.h"
class EA:
	public LSMethods
{
public:
	EA();
	~EA();
	void SetParameters_EA(int num_pop, long CPUtime);
	int EA_Run();
	int EA_Run2(vector<int>& _result);
protected:

	int m_num_pop;
	vector<vector<vector<int>>> m_popu_sol;
	vector<vector<int>> m_popu_factft;
	vector<int> m_popu_totaltft;

	vector<vector<vector<vector<int>>>> m_popu_perEij;

	bool m_flag_improve;

	long m_CPUtime;

	int InitialPopulation();
	void PopulationUpdate(const vector<vector<int>>& NewSol, const vector<int>& NewFacTFT, int NewTotalTFT,
		const vector<vector<vector<int>>>& perEij);

	bool JudgeSameFromPopulation(const vector<vector<int>>& Sol, const vector<int>& FacTFT, int TotalTFT);
	void Mutation(vector<vector<int>>& evo_sol, vector<int>& evo_factft, int& evo_totaltft,
				vector<vector<vector<int>>>& perEij);

};

