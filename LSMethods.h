#pragma once
#include "ConstructiveHeuristics.h"
class LSMethods:
	public ConstructiveHeuristics
{
public:
	LSMethods();
	void SetStartTime(long starttime);
	long GetStartTime();
	int Con_Heu_LS_Run(int type_conheu);
protected:
	
	//Swap
	void LS_SwapAmongAllFactories_Reference(const vector<vector<int>>& Lamda_Seq,
		vector<vector<int>>& Sol, vector<int>& FacTFT, int& TotalTFT, vector<vector<vector<int>>>& perEij);
	void LS_SwapAmongAllFactories_Reference(const vector<vector<int>>& Lamda_Seq,
		const vector<int>& RefFacTFT, vector<vector<int>>& Sol, vector<int>& FacTFT, int& TotalTFT, vector<vector<vector<int>>>& perEij);
	void LS_SwapAmongAllFactories_Reference(const vector<vector<int>>& Lamda_Seq,
		vector<vector<int>>& Sol, vector<int>& FacTFT, int& TotalTFT, vector<vector<vector<int>>>& perEij, int& cur_pos);
	
	void LS_SwapInFactory(vector<int>& Sol, int& FacTFT, int &Totaltft, vector<vector<int>>& Eij);

	//Insert
	void LS_InsertAmongAllFactories_Reference(const vector<vector<int>>& Lamda_Seq,
		vector<vector<int>>& Sol, vector<int>& FacTFT, int& TotalTFT, vector<vector<vector<int>>>& perEij);
	void LS_InsertAmongAllFactories_Reference(const vector<vector<int>>& Lamda_Seq,
		vector<vector<int>>& Sol, vector<int>& FacTFT, int& TotalTFT, vector<vector<vector<int>>>& perEij, int& cur_pos);

	void LS_InsertInFactory(vector<int>& Sol, int& FacTFT, int& Totaltft, vector<vector<int>>& Eij);

	//Hybrid
	void LS_HybridAmongAllFactories_Reference(const vector<vector<int>>& Lamda_Seq, float lamda,
		vector<vector<int>>& Sol, vector<int>& FacTFT, int& TotalTFT, vector<vector<vector<int>>>& perEij);

	void LS3_IG2S(vector<vector<int>>& sol, vector<int>& factft, vector<int>& facspan, int& totaltft, vector<vector<vector<int>>>& perEij);

	void LS1_IG2S(vector<int>& seq, int& tft, int& span, vector<vector<int>>& Eij);


	void LocalSearch_HEDFOA(vector<vector<int>>& sol, vector<int>& factft, vector<int>& facspan, int& totaltft, vector<vector<vector<int>>>& perEij);

	long m_elapCPUtime;

private:
	long m_starttime;
};

