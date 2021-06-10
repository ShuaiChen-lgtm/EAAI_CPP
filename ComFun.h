#pragma once
#include "Base.h"
#include "Problem.h"
#include "stdafx.h"
class ComFun :
	virtual public Problem
{
public:
	ComFun();
	~ComFun();
protected:
	void GetOrderSeq(int num_jobs, vector<int>& OrderSeq);
	void GetRandomSeq(int num_jobs, vector<int>& RefSeq);

	bool JudgeExist(int Job, const vector<int>& Terget_Seq);

	bool JudgeSameForTwoSol(const vector<vector<int>>& _1st_Sol, const vector<vector<int>>& _2nd_Sol);

	bool JudgeSameForTwoSol(const vector<vector<int>>& _1st_Sol, const vector<int>& _1st_factft, int _1st_TotalTFT, const vector<vector<int>>& _2nd_Sol, const vector<int>& _2nd_factft, int _2nd_TotalTFT);

	float JudgeSimilarityForTwoSol(const vector<vector<int>>& _1st_Sol, const vector<vector<int>>& _2nd_Sol);

	void GenerateRowByReduceDimen(const vector<vector<int>>& PreSeq, vector<int>& RefSeq);

	bool LocateJobInTwoDimArray(int Job, const vector<vector<int>>& TwoDimArray, int& FirstDimPos, int& SecondDimPos);

	void Check(const vector<vector<int>>& Sol, const vector<int>& FacTFT, int TotalTFT);

	int GetTotalFlowtime(const vector<int>& Seq);
	int GetTotalFlowtime(const vector<int>& Seq, int& Span);
	int GetTotalFlowtime(const vector<vector<int>>& Sol, vector<int>& FacTFT, vector<int>& FacSpan);
	int GetTotalFlowtime(const vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij);
	int GetEij(const vector<int>& Seq, vector<vector<int>>& Eij);
	int GetEij(const vector<int> &Seq, vector<vector<int>>& Eij, int& MakeSpan);

	int EraseAJob(vector<int>& Seq, int ErasePos, vector<vector<int>>& Eij);

	int EraseAJob_2(const vector<vector<int>>& PreEij, vector<int>& Seq, int ErasePos, vector<vector<int>>& Eij);

	int GetTFTAfterPushAJob(const vector<int>& Seq, int PushJob, const vector<vector<int>>& Eij);

	int GetTFTAfterInsertAJob(const vector<int>& Seq, int InsertJob, int InsertPos, const vector<vector<int>>& PreEij);

	void GetEijAfterInsertAJob(const vector<int>& Seq, int InsertJob, int InsertPos, vector<vector<int>>& Eij, int& Tft);


	int GetTFTAfterSwapAJob(const vector<int>& Seq, int SwapPos, int SwapJob, const vector<vector<int>>& Eij);
	void GetEijAfterSwapAJob(const vector<int>& Seq, int SwapPos, int SwapJob, vector<vector<int>>& Eij, int& Tft);

	int GetTFTAfterSwapJobsInFac(vector<int>& Seq, int ActivePos, int PassivePos, const vector<vector<int>>& Eij);

	int GetMinTFTAfterInsertSingleFac(const vector<int>& PreSeq, int InJob, int& BestPos, const vector<vector<int>>& Eij);

	int GetTFAfterInsertAJobInsideFactory(const vector<int>& Seq, int MovePos, int InsertPos, const vector<vector<int>>& Eij);

	int FindBestInsertPosInsideFactory(vector<int>& Seq, int MovePos, const vector<vector<int>>& Eij, int& BestPos);

	int FindBestInsertPosFacByMinIncreasedTFT(const vector<vector<int>>& Sol, const vector<int>& FacTFT,
		const vector<vector<vector<int>>>& perEij, int InJob, int& BestFac, int& BestPos);

	int GetBestSwapPosInOriginalFactory(vector<int>& PreSeq, int SwapPos, vector<vector<int>> Eij, int PreTFT, int& BestSwapPos);

	int GetBestSwapPosAfterSwapAmongTwoFactory(const vector<int>& ActiveSeq, const vector<vector<int>>& ActiveEij, int PreTFT_ActiveSeq, int SwapPos,
		const vector<int>& PassiveSeq, const vector<vector<int>>& PassiveEij, int PreTFT_PassiveSeq, int& BestSwapPos);

	template <typename T>
	void DescendingOrder(int Number, const vector<T>& GistSeq, vector<int>& Seq);

	template <typename T>
	void IncreasedOrder(int Number, const vector<T>& GistSeq, vector<int>& Seq);
};


template<typename T>
inline void ComFun::DescendingOrder(int Number, const vector<T>& GistSeq, vector<int>& Seq)
{
	Base::Pair<T>* ch = new Base::Pair<T>[Number];
	for (int j = 0; j < Number; j++)
	{
		ch[j].dim = j;
		ch[j].value = GistSeq[j];
	}
	sort(ch, ch + Number, Base::PairGreater<T>());
	for (int j = 0; j < Number; j++)
		Seq[j] = ch[j].dim;
	delete[]ch;
}

template<typename T>
inline void ComFun::IncreasedOrder(int Number, const vector<T>& GistSeq, vector<int>& Seq)
{
	Base::Pair<T>* ch = new Base::Pair<T>[Number];
	for (int j = 0; j < Number; j++)
	{
		ch[j].dim = j;
		ch[j].value = GistSeq[j];
	}
	sort(ch, ch + Number, Base::PairLess<T>());
	for (int j = 0; j < Number; j++)
		Seq[j] = ch[j].dim;
	delete[]ch;
}

