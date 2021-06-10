#pragma once
#include "Operation.h"
class ConstructiveHeuristics :
	virtual public Operation
{
public:

	int Con_Heu_Run(int type_conheu);
	float lamda_DPEF;
protected:

	int DLR(vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij);
	int DLR_DNEH(float Ratio, vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij);
	int NEH2_en(vector<vector<int>>& Sol, vector<int>& FacTFT, vector<int>& FacSpan, vector<vector<vector<int>>>& perEij);
	int DNEH(vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij);

	int RandDNEH(vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij);

	int RandDLR(vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij);

	int DPFE(vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij);
	int RandDPFE(vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij);

	int DPF(vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij);

	int HPF23(float lamda, float mu, vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij);

	void GenerateSeedSeq_EHPF2(vector<int>& SeedSeq, float lamda);

	int EHPF2(vector<vector<int>>& Sol, vector<int>& FacTFT);

	int EHPF2(vector<vector<int>>& Sol, vector<int>& FacTFT, vector<int>& FacSpan);

	int NEHR2A4(vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij);

	int RandNEHR2A4(vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij);

	int PWNEH(vector<vector<int>>& Sol, vector<int>& FacTFT, vector<int>& FacSpan, vector<vector<vector<int>>>& perEij);

private:
	int FindBestJobPosFor_DLR(const vector<int>& CurSeq, const vector<int>& UnscheduledSeq, int& ResultMakeSpan,
		const vector<vector<int>>& Eij);
	void GenerateSeedSeq_IF0(vector<int>& SeedSeq);

	void GenerateSeedSeqSingle_IF0(vector<int>& SeedSeq);

	void GenerateSeedSeqForDPFE(vector<int>& SeedSeq);
	int FindBestJobPosForDPFE(const vector<int>& CurSeq, const vector<int>& UnscheduledSeq, 
		const vector<vector<int>>& Eij, int& ResultMakeSpan);

	int PWIndexFun(vector<vector<int>> Eij, vector<int> Remain_Seq, int k, vector<int>& RemainPtime);
	void PW(vector<int>& Seq);

	int FindBestJobPosFor_DPWE(const vector<int>& CurSeq, const vector<int>& UnscheduledSeq, 
		const vector<vector<int>>& Eij, int& ResultMakeSpan);

	void Rule_HPF23(vector<int>& Seq, float lamda, float mu);

	void Rule_HPF23_FirstJob(vector<int>& Permutation, float lamda);

	int FindBestJobPosForHPF23(const vector<vector<int>>& Eij, const vector<int>& Remain_Seq, int k, float mu);

	int FindBestJobPosForEHPF2(const vector<vector<int>>& Eij, const vector<int>& RemainSeq, int NextPos, float mu);

};

