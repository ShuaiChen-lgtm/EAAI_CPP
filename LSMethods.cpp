#include "LSMethods.h"

LSMethods::LSMethods()
{
	m_starttime = 0;
}

void LSMethods::SetStartTime(long starttime)
{
	m_starttime = starttime;
}

long LSMethods::GetStartTime()
{
	return m_starttime;
}

int LSMethods::Con_Heu_LS_Run(int type_conheu)
{
	GetTotalPTime();
	vector<vector<vector<int>>> perEij(m_Factories);
	switch (type_conheu)
	{
	case 1:this->m_BestSol_TotalTFT = NEHR2A4(this->m_BestSol_Seq, this->m_BestSol_FacTFT, perEij); break;
	case 2:this->m_BestSol_TotalTFT = DNEH(this->m_BestSol_Seq, this->m_BestSol_FacTFT, perEij); break;
	case 3:this->m_BestSol_TotalTFT = DLR(this->m_BestSol_Seq, this->m_BestSol_FacTFT, perEij); break;
	case 4:this->m_BestSol_TotalTFT = DLR_DNEH(0.3, this->m_BestSol_Seq, this->m_BestSol_FacTFT, perEij); break;//0.3
	case 5:this->m_BestSol_TotalTFT = DPFE(this->m_BestSol_Seq, this->m_BestSol_FacTFT, perEij); break;
	case 6:this->m_BestSol_TotalTFT = DPF(this->m_BestSol_Seq, this->m_BestSol_FacTFT, perEij); break;
	default:
		break;
	}
	LS_SwapAmongAllFactories_Reference(m_BestSol_Seq, m_BestSol_Seq, m_BestSol_FacTFT,
		m_BestSol_TotalTFT, perEij);
	Check(this->m_BestSol_Seq, this->m_BestSol_FacTFT, this->m_BestSol_TotalTFT);
	return this->m_BestSol_TotalTFT;
}

//Swap
void LSMethods::LS_SwapAmongAllFactories_Reference(const vector<vector<int>>& Lamda_Seq,
	vector<vector<int>>& Sol, vector<int>& FacTFT, int& TotalTFT, vector<vector<vector<int>>>& perEij)
{
	vector<int> RefSeq(this->m_Jobs);
	this->GenerateRowByReduceDimen(Lamda_Seq, RefSeq);

	int Cnt = 0;
	int j = 0;
	while (Cnt < this->m_Jobs)
	{
		// 此次移除RefSeq中第 Cnt 个工件，先找到在当前解中的位置
		int RemoveJob = RefSeq[j];
		int Factory = -1, Position = -1;
		this->LocateJobInTwoDimArray(RemoveJob, Sol, Factory, Position);

		bool ImpFlag = this->CarryBestSwapByMinGlobalTFT(Factory, Position, Sol, FacTFT, TotalTFT, perEij);

		if (ImpFlag)
		{
			Cnt = 0;
		}
		else
			Cnt++;
		j = (j + 1) % this->m_Jobs;
	}
}

void LSMethods::LS_SwapAmongAllFactories_Reference(const vector<vector<int>>& Lamda_Seq, const vector<int>& RefFacTFT,
	vector<vector<int>>& Sol, vector<int>& FacTFT, int& TotalTFT, vector<vector<vector<int>>>& perEij)
{
	vector<int> RefSeq(this->m_Jobs);
	this->GenerateRowByReduceDimen(Lamda_Seq, RefSeq);

	int Cnt = 0;
	int j = 0;
	while (Cnt < this->m_Jobs)
	{
		// 此次移除RefSeq中第 Cnt 个工件，先找到在当前解中的位置
		int RemoveJob = RefSeq[j];
		int Factory = -1, Position = -1;
		this->LocateJobInTwoDimArray(RemoveJob, Sol, Factory, Position);

		bool ImpFlag = this->CarryBestSwapByMinGlobalTFT(Factory, Position, Sol, FacTFT, TotalTFT, perEij);

		if (ImpFlag)
		{
			Cnt = 0;
		}
		else
			Cnt++;
		j = (j + 1) % this->m_Jobs;
	}
}

void LSMethods::LS_SwapAmongAllFactories_Reference(const vector<vector<int>>& Lamda_Seq,
	vector<vector<int>>& Sol, vector<int>& FacTFT, int& TotalTFT, vector<vector<vector<int>>>& perEij, int& cur_pos)
{
	vector<int> RefSeq(this->m_Jobs);
	this->GenerateRowByReduceDimen(Lamda_Seq, RefSeq);

	int Cnt = 0;
	while (Cnt < this->m_Jobs)
	{
		// 此次移除RefSeq中第 Cnt 个工件，先找到在当前解中的位置
		int RemoveJob = RefSeq[cur_pos];
		int Factory = -1, Position = -1;
		this->LocateJobInTwoDimArray(RemoveJob, Sol, Factory, Position);

		bool ImpFlag = this->CarryBestSwapByMinGlobalTFT(Factory, Position, Sol, FacTFT, TotalTFT, perEij);

		if (ImpFlag)
		{
			Cnt = 0;
		}
		else
			Cnt++;
		cur_pos++;
		if (cur_pos >= m_Jobs)
			cur_pos = 0;
	}
}

void LSMethods::LS_SwapInFactory(vector<int>& Sol, int& FacTFT, int& Totaltft, vector<vector<int>>& Eij)
{
	int nJobs = Sol.size();
	vector<int> RefSeq(nJobs);
	GetOrderSeq(nJobs, RefSeq);
	random_shuffle(RefSeq.begin(), RefSeq.end());
	int cnt = 0;
	while (cnt < nJobs)
	{
		int pos = RefSeq[cnt];
		int bestswappos = -1;
		int mindeceasetf = GetBestSwapPosInOriginalFactory(Sol, pos, Eij, FacTFT, bestswappos);
		if (mindeceasetf > 0)
		{
			SwapTwoJobsInsideFactory(pos, bestswappos, Sol, Eij, FacTFT, Totaltft);
			cnt = 0;
			random_shuffle(RefSeq.begin(), RefSeq.end());
		}
		else
			cnt++;
	}
}

//Insert
void LSMethods::LS_InsertAmongAllFactories_Reference(const vector<vector<int>> &Lamda_Seq, 
	vector<vector<int>>& Sol, vector<int>& FacTFT, int& TotalTFT, vector<vector<vector<int>>>& perEij)
{
	vector<int> RefSeq(this->m_Jobs);
	this->GenerateRowByReduceDimen(Lamda_Seq, RefSeq);
	int Cnt = 0;
	int j = 0;
	while (Cnt < this->m_Jobs)
	{
		// 此次移除RefSeq中第 j 个工件，先找到在当前解中的位置
		int RemoveJob = RefSeq[j];//--------------------------------------------------
		int Factory = -1, Position = -1;
		this->LocateJobInTwoDimArray(RemoveJob, Sol, Factory, Position);
		if (Sol[Factory].size() == 1)
		{
			Cnt++;
			continue;
		}
		bool ImpFlag = this->CarryBestInsertByMinGlobalTFT(Factory, Position, Sol, FacTFT, TotalTFT, perEij);
		if (ImpFlag)
		{
			Cnt = 0;
			if (Base::GetElapsedProcessTime() - GetStartTime() > m_elapCPUtime)
				break;
		}
		else
			Cnt++;
		j = (j + 1) % this->m_Jobs;
	}
}

void LSMethods::LS_InsertAmongAllFactories_Reference(const vector<vector<int>> &Lamda_Seq, 
	vector<vector<int>>& Sol, vector<int>& FacTFT, int& TotalTFT, vector<vector<vector<int>>>& perEij, int& cur_pos)
{
	vector<int> RefSeq(this->m_Jobs);
	this->GenerateRowByReduceDimen(Lamda_Seq, RefSeq);
	int Cnt = 0;
	while (Cnt < this->m_Jobs)
	{
		// 此次移除RefSeq中第 j 个工件，先找到在当前解中的位置
		int RemoveJob = RefSeq[cur_pos];//--------------------------------------------------
		int Factory = -1, Position = -1;
		this->LocateJobInTwoDimArray(RemoveJob, Sol, Factory, Position);
		if (Sol[Factory].size() == 1)
		{
			Cnt++;
			continue;
		}
		bool ImpFlag = this->CarryBestInsertByMinGlobalTFT(Factory, Position, Sol, FacTFT, TotalTFT, perEij);
		if (ImpFlag)
		{
			Cnt = 0;
		}
		else
			Cnt++;
		cur_pos++;
		if (cur_pos >= m_Jobs)
			cur_pos = 0;
	}
}

void LSMethods::LS_InsertInFactory(vector<int>& Sol, int& FacTFT, int& Totaltft, vector<vector<int>>& Eij)
{
	Totaltft -= FacTFT;
	int nJobs = Sol.size();
	vector<int> RefSeq(nJobs);
	GetOrderSeq(nJobs, RefSeq);
	random_shuffle(RefSeq.begin(), RefSeq.end());
	int cnt = 0;
	while (cnt < nJobs)
	{
		int pos = RefSeq[cnt];
		int job = Sol[pos];

		vector<vector<int>> newEij;
		Sol.erase(Sol.begin() + pos);
		GetEij(Sol, newEij);

		int bestpos = -1;
		int mintft = GetMinTFTAfterInsertSingleFac(Sol, job, bestpos, newEij);
		if (mintft< FacTFT)
		{
			InsertAJobToSeq(job, bestpos, Sol, newEij, FacTFT);
			Eij = newEij;
			cnt = 0;
			random_shuffle(RefSeq.begin(), RefSeq.end());
		}
		else
		{
			Sol.insert(Sol.begin() + pos, job);
			cnt++;
		}
	}
	Totaltft += FacTFT;
}

//Hybrid
void LSMethods::LS_HybridAmongAllFactories_Reference(const vector<vector<int>> &Lamda_Seq, float lamda, 
	vector<vector<int>>& Sol, vector<int>& FacTFT, int& TotalTFT, vector<vector<vector<int>>>& perEij)
{
	vector<int> RefSeq(this->m_Jobs);
	//this->GenerateRowByReduceDimen(Lamda_Seq, RefSeq);
	this->GenerateRowByReduceDimen(Sol, RefSeq);

	int Cnt = 0;
	int j = 0;
	while (Cnt < this->m_Jobs)
	{
		// 此次移除RefSeq中第 j 个工件，先找到在当前解中的位置
		int RemoveJob = RefSeq[j];//--------------------------------------------------test
		int Factory = -1, Position = -1;
		this->LocateJobInTwoDimArray(RemoveJob, Sol, Factory, Position);
		if (Sol[Factory].size() == 1)
		{
			Cnt++;
			continue;
		}
		bool ImpFlag;
		float _randnum = (float)rand() / RAND_MAX;
		if (_randnum < lamda)
			ImpFlag = this->CarryBestSwapByMinGlobalTFT(Factory, Position, Sol, FacTFT, TotalTFT, perEij);
		else
			ImpFlag = this->CarryBestInsertByMinGlobalTFT(Factory, Position, Sol, FacTFT, TotalTFT, perEij);
		if (ImpFlag)
		{
			Cnt = 0;
			if (Base::GetElapsedProcessTime() - GetStartTime() > m_elapCPUtime)
				break;
		}
		else
			Cnt++;
		j = (j + 1) % this->m_Jobs;
	}
}

// IG_2S
void LSMethods::LS3_IG2S(vector<vector<int>>& sol, vector<int>& factft, vector<int>& facspan, int& totaltft,
	vector<vector<vector<int>>>& perEij)
{
	//detect the factory that has the maximum span
	int extfac = max_element(facspan.begin(), facspan.end()) - facspan.begin();

	int num_extjobs = sol[extfac].size();
	vector<int> RefSeq(num_extjobs);
	GetOrderSeq(num_extjobs, RefSeq);
	random_shuffle(RefSeq.begin(), RefSeq.end());
	int Cnt = 0;
	int k = 0;
	while (Cnt < num_extjobs)
	{
		int pos = RefSeq[k];
		bool ImpFlag = this->CarryBestInsertByMinGlobalTFT(extfac, pos, sol, factft, facspan, totaltft, perEij);
		if (ImpFlag)
		{
			//detect the factory that has the maximum span
			extfac = max_element(facspan.begin(), facspan.end()) - facspan.begin();
			num_extjobs = sol[extfac].size();

			GetOrderSeq(num_extjobs, RefSeq);
			random_shuffle(RefSeq.begin(), RefSeq.end());
			Cnt = 0;
		}
		else
		{
			Cnt++;
		}
		k = (k + 1) % num_extjobs;
	}
}

void LSMethods::LS1_IG2S(vector<int>& seq, int& tft, int& span, vector<vector<int>>& Eij)
{
	int sizeofseq = seq.size();
	vector<int> RefSeq(sizeofseq);
	GetOrderSeq(sizeofseq, RefSeq);
	random_shuffle(RefSeq.begin(), RefSeq.end());
	for (int j = 0; j < sizeofseq; j++)
	{
		int pos = RefSeq[j];
		int job = seq[pos];

		int BPos = -1;
		int remintft = FindBestInsertPosInsideFactory(seq, pos, Eij, BPos);
		if (remintft < tft)
		{
			seq.erase(seq.begin() + pos);
			if (BPos > pos)
				seq.insert(seq.begin() + (BPos - 1), job);
			else
				seq.insert(seq.begin() + BPos, job);
			tft = remintft;
			GetEij(seq, Eij);
			span = Eij[sizeofseq - 1][m_Machines - 1];
		}
	}
}

void LSMethods::LocalSearch_HEDFOA(vector<vector<int>>& sol, vector<int>& factft, vector<int>& facspan, int& totaltft,
	vector<vector<vector<int>>>& perEij)
{
	//detect the factory that has the maximum span
	int extfac = max_element(facspan.begin(), facspan.end()) - facspan.begin();

	int num_extjobs = sol[extfac].size();
	vector<int> RefSeq(num_extjobs);
	GetOrderSeq(num_extjobs, RefSeq);
	random_shuffle(RefSeq.begin(), RefSeq.end());
	int Cnt = 0;
	while (Cnt < num_extjobs)
	{
		int pos = RefSeq[Cnt];
		bool ImpFlag = this->CarryBestInsertByMinGlobalTFT(extfac, pos, sol, factft, facspan, totaltft, perEij);
		if (ImpFlag)
		{
			//detect the factory that has the maximum span
			extfac = max_element(facspan.begin(), facspan.end()) - facspan.begin();
			num_extjobs = sol[extfac].size();

			GetOrderSeq(num_extjobs, RefSeq);
			random_shuffle(RefSeq.begin(), RefSeq.end());
			Cnt = 0;
		}
		else
		{
			Cnt++;
		}
	}
}


