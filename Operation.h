#pragma once
#include "ComFun.h"
class Operation :
	public ComFun
{
public:
	Operation();
	~Operation();
protected:

	void PushAJobToSeq(int Job, vector<int>& Seq, vector<vector<int>>& Eij, int& totaltftj);
	void InsertAJobToSeq(int InsertJob, int InsertPos, vector<int>& Seq, vector<vector<int>>& Eij, int& tft);
	void InsertAJobToItsSeq(int movepos, int tergetpos, vector<int>& Seq, vector<vector<int>>& Eij, int& tft);

	void SwapTwoJobsInsideFactory(int pos1, int pos2, vector<int>& Seq, vector<vector<int>>& Eij, int& Tft, int& TotalTFT);
	void SwapTwoJobsInsideFactory(int pos1, int pos2, vector<int>& Seq, vector<vector<int>>& Eij);
	void SwapTwoJobsBetweenFactories(int fac1, int fac2, int pos1, int pos2,
		vector<vector<int>>& Sol, vector<vector<vector<int>>>& perEij, vector<int>& FacTFT, int& Totaltft);
	bool CarryBestSwapByMinGlobalTFT(int SwapFac, int SwapPos, vector<vector<int>>& Sol, vector<int>& FacTFT, int& TotalTFT, vector<vector<vector<int>>>& perEij);
	bool CarryBestInsertByMinGlobalTFT(int MoveFac, int MovePos, vector<vector<int>>& Sol, vector<int>& FacTFT, int& TotalTFT, vector<vector<vector<int>>>& perEij);
	bool CarryBestInsertByMinGlobalTFT(int MoveFac, int MovePos, vector<vector<int>>& Sol, vector<int>& FacTFT, vector<int>& FacSpan, int& TotalTFT, vector<vector<vector<int>>>& perEij);
	bool CarryBestInsertByMinGlobalTFT_randfactory(int MoveFac, int MovePos, vector<vector<int>>& Sol, vector<int>& FacTFT, int& TotalTFT, vector<vector<vector<int>>>& perEij);
	bool CarryBestPermutationByMinGlobalTFT_randfactory(int plimit, int MoveFac, int MovePos, vector<vector<int>>& Sol, vector<int>& FacTFT, int& TotalTFT, vector<vector<vector<int>>>& perEij);
};
