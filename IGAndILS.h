#pragma once
#include "LSMethods.h"
class IGAndILS:
	public LSMethods
{
public:
	void SetParameters_IG(long C_of_CPUtime, int Setd, float SetT, int Type_Initial);
	void SetParametersSecond_ILS(long C_of_CPUtime, int SetPi, int Setd, float SetT, int Type_Initial);
	int Initialization(int Type_Init, vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij);
	int IG_Run();  //IGPG
	int ILS_Run(); //ILS
protected:

	void LocalSearchType(vector<vector<int>>& Sol, vector<int>& FacTFT, int& TotalTFT, vector<vector<vector<int>>>& perEij);
	void Destruction_IG(vector<vector<int>>& Sol, vector<int>& FacTFT, vector<int>& Des_Seq, vector<vector<vector<int>>>& perEij);
	int Reconstruction_IG(vector<int>& Des_Seq, vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij);

	void Disturbance_ILS(vector<vector<int>>& Sol, vector<int>& FacTFT, int& TotalTFT, vector<vector<vector<int>>>& perEij);


private:

	int m_d;
	float m_T;
	long m_CPUtime;
	int m_Pi;
	int m_Type_initial;
	vector < vector<vector<int>>> m_perEij;
};

