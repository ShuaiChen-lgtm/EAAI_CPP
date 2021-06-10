#pragma once
#include "LSMethods.h"
class MyPBIG:
	public LSMethods
{
public:
	MyPBIG();
	void SetParameters_PBIG(int num_pop, long CPUtime, int d);
	void SetParameters2_PBIG(int num_pop, long CPUtime, int d, float lamda);
	void SetParameters3_PBIG(int num_pop, int Type_Evolution, long CPUtime, int d, float lamda);
	int MyPBIG_Run(string typeInitial);
protected:

	int InitialPopulationDLR();
	int InitialPopulationDPFE();
	int InitialPopulationDNEH();
	void PopulationUpdate(const vector<vector<int>>& NewSol, const vector<int>& NewFacTFT, int NewTotalTFT,
		const vector<vector<vector<int>>>& perEij);
	void PopulationUpdate2(const vector<vector<int>>& NewSol, const vector<int>& NewFacTFT, int NewTotalTFT, const vector<vector<vector<int>>>& perEij);
	void PopulationUpdate22(int parent, const vector<vector<int>>& NewSol, const vector<int>& NewFacTFT, int NewTotalTFT, const vector<vector<vector<int>>>& perEij);
	void PopulationUpdate3(const vector<vector<int>>& NewSol, const vector<int>& NewFacTFT, int NewTotalTFT, const vector<vector<vector<int>>>& perEij);
	void PopulationUpdate4(const vector<vector<int>>& NewSol, const vector<int>& NewFacTFT, int NewTotalTFT, const vector<vector<vector<int>>>& perEij);
	bool JudgeSameFromPopulation(const vector<vector<int>>& Sol, const vector<int>& FacTFT, int TotalTFT);
	bool Judgesame();
	void ChooseExtJob(vector<int>& Des_Seq);

	void EvaluationDR1(vector<vector<int>>& evo_sol, vector<int>& evo_factft, int& evo_totaltft, vector<vector<vector<int>>>& perEij);
	void EvaluationDR2(vector<vector<int>>& evo_sol, vector<int>& evo_factft, int& evo_totaltft,
		vector<vector<vector<int>>>& perEij);
	void EvaluationDR3(vector<vector<int>>& evo_sol, vector<int>& evo_factft, int& evo_totaltft, 
		vector<vector<vector<int>>>& perEij);

	int Destruction_Unknown_PBIG(vector<int>& Des_Seq, vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij);

	int Destruction_Known_PBIG(const vector<int>& Des_Seq, vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij);

	int Reconstruction_PBIG(vector<int>& Des_Seq, vector<vector<int>>& Sol, vector<int>& FacTFT,
		vector<vector<vector<int>>>& perEij);
private:

	int m_d;
	int m_num_pop;
	float m_lamda;
	vector<vector<vector<int>>> m_popu_sol;
	vector<vector<int>> m_popu_factft;
	vector<int> m_popu_totaltft;

	vector<vector<vector<vector<int>>>> m_popu_perEij;

	bool m_flag_improve;

	long m_CPUtime;
	int m_Evolution_Type;
};

