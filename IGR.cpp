#include "IGR.h"

void IGR::SetParameter(int d, int Limp, float T, int time_factor)
{
	m_d=d;
	m_limp = Limp;
	m_T = T;
	m_CPUtime = time_factor * m_Jobs * m_Machines;

	m_perEij.resize(m_Factories);
}

int IGR::IGR_Run(float lamda,float mu)
{
	vector<int> ConvergentArray;
	SetStartTime(Base::GetElapsedProcessTime());
	this->GetTotalPTime();
	this->GetTotalLoad();
	float Temperature = this->m_T * this->m_TotalLoad / (10.0 * this->m_Jobs * this->m_Machines);
	this->m_BestSol_TotalTFT = HPF23(lamda, mu, this->m_BestSol_Seq, this->m_BestSol_FacTFT, m_perEij);
	VNS(this->m_BestSol_Seq, this->m_BestSol_FacTFT, m_BestSol_TotalTFT, m_perEij);

	vector<vector<int>> CurSol(this->m_BestSol_Seq);
	vector<int> CurFacTFT(this->m_BestSol_FacTFT);
	int CurTotalTFT = this->m_BestSol_TotalTFT;
	int Cnt = 0;
	while (Base::GetElapsedProcessTime() - GetStartTime() < m_CPUtime)
	{
		vector<vector<int>> Evo_Sol(CurSol);
		vector<int> Evo_FacTFT(CurFacTFT);
		int Evo_TotalTFT = CurTotalTFT;
		vector<vector<vector<int>>> Evo_perEij(this->m_perEij);

		Evo_TotalTFT = Perturbation(Evo_Sol, Evo_FacTFT, Evo_perEij);
		VNS(Evo_Sol, Evo_FacTFT, Evo_TotalTFT, Evo_perEij);
		Reassignment(Evo_Sol, Evo_FacTFT, Evo_TotalTFT, Evo_perEij);
		Permutation(Evo_Sol, Evo_FacTFT, Evo_TotalTFT, Evo_perEij);

		//更新解
		if (Evo_TotalTFT < CurTotalTFT)
		{
			CurSol = Evo_Sol;
			CurFacTFT = Evo_FacTFT;
			CurTotalTFT = Evo_TotalTFT;
			m_perEij = Evo_perEij;
			if (Evo_TotalTFT < this->m_BestSol_TotalTFT)
			{
				this->m_BestSol_Seq = Evo_Sol;
				this->m_BestSol_FacTFT = Evo_FacTFT;
				this->m_BestSol_TotalTFT = Evo_TotalTFT;
			}
		}
		else
		{
			float RandNumber = (float)rand() / RAND_MAX;
			float Index = 1.00 * exp((CurTotalTFT - Evo_TotalTFT) / Temperature);
			if (RandNumber < Index)
			{
				CurSol = Evo_Sol;
				CurFacTFT = Evo_FacTFT;
				CurTotalTFT = Evo_TotalTFT;
				m_perEij = Evo_perEij;
			}
		}
		Cnt++;
	}
	cout << Cnt << ends << (Base::GetElapsedProcessTime() - GetStartTime()) << ends;
	Check(m_BestSol_Seq, m_BestSol_FacTFT, m_BestSol_TotalTFT);
	return m_BestSol_TotalTFT;
}


void IGR::IGR_Run2(float lamda, float mu, vector<int>& _results)
{
	vector<int> ConvergentArray;
	SetStartTime(Base::GetElapsedProcessTime());
	this->GetTotalPTime();
	this->GetTotalLoad();
	float Temperature = this->m_T * this->m_TotalLoad / (10.0 * this->m_Jobs * this->m_Machines);
	this->m_BestSol_TotalTFT = HPF23(lamda, mu, this->m_BestSol_Seq, this->m_BestSol_FacTFT, m_perEij);
	VNS(this->m_BestSol_Seq, this->m_BestSol_FacTFT, m_BestSol_TotalTFT, m_perEij);

	vector<vector<int>> CurSol(this->m_BestSol_Seq);
	vector<int> CurFacTFT(this->m_BestSol_FacTFT);
	int CurTotalTFT = this->m_BestSol_TotalTFT;

	m_CPUtime = m_Jobs * m_Machines * 10;
	while (Base::GetElapsedProcessTime() - GetStartTime() < m_CPUtime)
	{
		vector<vector<int>> Evo_Sol(CurSol);
		vector<int> Evo_FacTFT(CurFacTFT);
		int Evo_TotalTFT = CurTotalTFT;
		vector<vector<vector<int>>> Evo_perEij(this->m_perEij);

		Evo_TotalTFT = Perturbation(Evo_Sol, Evo_FacTFT, Evo_perEij);
		VNS(Evo_Sol, Evo_FacTFT, Evo_TotalTFT, Evo_perEij);
		Reassignment(Evo_Sol, Evo_FacTFT, Evo_TotalTFT, Evo_perEij);
		Permutation(Evo_Sol, Evo_FacTFT, Evo_TotalTFT, Evo_perEij);

		//更新解
		if (Evo_TotalTFT < CurTotalTFT)
		{
			CurSol = Evo_Sol;
			CurFacTFT = Evo_FacTFT;
			CurTotalTFT = Evo_TotalTFT;
			m_perEij = Evo_perEij;
			if (Evo_TotalTFT < this->m_BestSol_TotalTFT)
			{
				this->m_BestSol_Seq = Evo_Sol;
				this->m_BestSol_FacTFT = Evo_FacTFT;
				this->m_BestSol_TotalTFT = Evo_TotalTFT;
			}
		}
		else
		{
			float RandNumber = (float)rand() / RAND_MAX;
			float Index = 1.00 * exp((CurTotalTFT - Evo_TotalTFT) / Temperature);
			if (RandNumber < Index)
			{
				CurSol = Evo_Sol;
				CurFacTFT = Evo_FacTFT;
				CurTotalTFT = Evo_TotalTFT;
				m_perEij = Evo_perEij;
			}
		}
	}
	_results[0] = m_BestSol_TotalTFT;

	m_CPUtime = m_Jobs * m_Machines * 20;
	while (Base::GetElapsedProcessTime() - GetStartTime() < m_CPUtime)
	{
		vector<vector<int>> Evo_Sol(CurSol);
		vector<int> Evo_FacTFT(CurFacTFT);
		int Evo_TotalTFT = CurTotalTFT;
		vector<vector<vector<int>>> Evo_perEij(this->m_perEij);

		Evo_TotalTFT = Perturbation(Evo_Sol, Evo_FacTFT, Evo_perEij);
		VNS(Evo_Sol, Evo_FacTFT, Evo_TotalTFT, Evo_perEij);
		Reassignment(Evo_Sol, Evo_FacTFT, Evo_TotalTFT, Evo_perEij);
		Permutation(Evo_Sol, Evo_FacTFT, Evo_TotalTFT, Evo_perEij);

		//更新解
		if (Evo_TotalTFT < CurTotalTFT)
		{
			CurSol = Evo_Sol;
			CurFacTFT = Evo_FacTFT;
			CurTotalTFT = Evo_TotalTFT;
			m_perEij = Evo_perEij;
			if (Evo_TotalTFT < this->m_BestSol_TotalTFT)
			{
				this->m_BestSol_Seq = Evo_Sol;
				this->m_BestSol_FacTFT = Evo_FacTFT;
				this->m_BestSol_TotalTFT = Evo_TotalTFT;
			}
		}
		else
		{
			float RandNumber = (float)rand() / RAND_MAX;
			float Index = 1.00 * exp((CurTotalTFT - Evo_TotalTFT) / Temperature);
			if (RandNumber < Index)
			{
				CurSol = Evo_Sol;
				CurFacTFT = Evo_FacTFT;
				CurTotalTFT = Evo_TotalTFT;
				m_perEij = Evo_perEij;
			}
		}
	}
	_results[1] = m_BestSol_TotalTFT;

	m_CPUtime = m_Jobs * m_Machines * 40;
	while (Base::GetElapsedProcessTime() - GetStartTime() < m_CPUtime)
	{
		vector<vector<int>> Evo_Sol(CurSol);
		vector<int> Evo_FacTFT(CurFacTFT);
		int Evo_TotalTFT = CurTotalTFT;
		vector<vector<vector<int>>> Evo_perEij(this->m_perEij);

		Evo_TotalTFT = Perturbation(Evo_Sol, Evo_FacTFT, Evo_perEij);
		VNS(Evo_Sol, Evo_FacTFT, Evo_TotalTFT, Evo_perEij);
		Reassignment(Evo_Sol, Evo_FacTFT, Evo_TotalTFT, Evo_perEij);
		Permutation(Evo_Sol, Evo_FacTFT, Evo_TotalTFT, Evo_perEij);

		//更新解
		if (Evo_TotalTFT < CurTotalTFT)
		{
			CurSol = Evo_Sol;
			CurFacTFT = Evo_FacTFT;
			CurTotalTFT = Evo_TotalTFT;
			m_perEij = Evo_perEij;
			if (Evo_TotalTFT < this->m_BestSol_TotalTFT)
			{
				this->m_BestSol_Seq = Evo_Sol;
				this->m_BestSol_FacTFT = Evo_FacTFT;
				this->m_BestSol_TotalTFT = Evo_TotalTFT;
			}
		}
		else
		{
			float RandNumber = (float)rand() / RAND_MAX;
			float Index = 1.00 * exp((CurTotalTFT - Evo_TotalTFT) / Temperature);
			if (RandNumber < Index)
			{
				CurSol = Evo_Sol;
				CurFacTFT = Evo_FacTFT;
				CurTotalTFT = Evo_TotalTFT;
				m_perEij = Evo_perEij;
			}
		}
	}
	_results[2] = m_BestSol_TotalTFT;

	m_CPUtime = m_Jobs * m_Machines * 80;
	while (Base::GetElapsedProcessTime() - GetStartTime() < m_CPUtime)
	{
		vector<vector<int>> Evo_Sol(CurSol);
		vector<int> Evo_FacTFT(CurFacTFT);
		int Evo_TotalTFT = CurTotalTFT;
		vector<vector<vector<int>>> Evo_perEij(this->m_perEij);

		Evo_TotalTFT = Perturbation(Evo_Sol, Evo_FacTFT, Evo_perEij);
		VNS(Evo_Sol, Evo_FacTFT, Evo_TotalTFT, Evo_perEij);
		Reassignment(Evo_Sol, Evo_FacTFT, Evo_TotalTFT, Evo_perEij);
		Permutation(Evo_Sol, Evo_FacTFT, Evo_TotalTFT, Evo_perEij);

		//更新解
		if (Evo_TotalTFT < CurTotalTFT)
		{
			CurSol = Evo_Sol;
			CurFacTFT = Evo_FacTFT;
			CurTotalTFT = Evo_TotalTFT;
			m_perEij = Evo_perEij;
			if (Evo_TotalTFT < this->m_BestSol_TotalTFT)
			{
				this->m_BestSol_Seq = Evo_Sol;
				this->m_BestSol_FacTFT = Evo_FacTFT;
				this->m_BestSol_TotalTFT = Evo_TotalTFT;
			}
		}
		else
		{
			float RandNumber = (float)rand() / RAND_MAX;
			float Index = 1.00 * exp((CurTotalTFT - Evo_TotalTFT) / Temperature);
			if (RandNumber < Index)
			{
				CurSol = Evo_Sol;
				CurFacTFT = Evo_FacTFT;
				CurTotalTFT = Evo_TotalTFT;
				m_perEij = Evo_perEij;
			}
		}
	}
	_results[3] = m_BestSol_TotalTFT;

	cout << (Base::GetElapsedProcessTime() - GetStartTime()) << ends;
	Check(m_BestSol_Seq, m_BestSol_FacTFT, m_BestSol_TotalTFT);
}

void IGR::ChooseExtJob(vector<int>& Des_Seq)
{
	for (int j = 0; j < m_d; j++)
	{
		int job = -1;
		do
		{
			job = rand() % m_Jobs;
		} while (JudgeExist(job, Des_Seq));
		Des_Seq.push_back(job);
	}
}

int IGR::Destruction(const vector<int>& Des_Seq, vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij)
{
	int sizeof_desseq = Des_Seq.size();
	vector<int> flag_change(m_Factories, false);
	for (int j = 0; j < sizeof_desseq; j++)
	{
		int job = Des_Seq[j];
		int factory = -1, position = -1;
		LocateJobInTwoDimArray(job, Sol, factory, position);
		Sol[factory].erase(Sol[factory].begin() + position);
		flag_change[factory] = true;
	}
	for (int f = 0; f < m_Factories; f++)
	{
		if (!flag_change[f])
			continue;
		FacTFT[f] = GetEij(Sol[f], perEij[f]);
	}
	return accumulate(FacTFT.begin(), FacTFT.end(), 0);
}

int IGR::Reconstruction(const vector<int>& Des_Seq, vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij)
{
	for (int j = 0; j < m_d; j++)
	{
		int Job = Des_Seq[j];
		int BestPos = -1, BestFac = -1;
		int Min_TFT = FindBestInsertPosFacByMinIncreasedTFT(Sol, FacTFT, perEij, Job, BestFac, BestPos);
		FacTFT[BestFac] = Min_TFT;
		InsertAJobToSeq(Job, BestPos, Sol[BestFac], perEij[BestFac], FacTFT[BestFac]);
	}
	return accumulate(FacTFT.begin(), FacTFT.end(), 0);
}

int IGR::Perturbation(vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij)
{
	vector<int> seq_des;
	ChooseExtJob(seq_des);
	Destruction(seq_des, Sol, FacTFT, perEij);
	int totaltft = Reconstruction(seq_des, Sol, FacTFT, perEij);
	return totaltft;
}

void IGR::VNS(vector<vector<int>>& sol, vector<int>& factft, int& totaltft, vector<vector<vector<int>>>& perEij)
{
	for (int f = 0; f < m_Factories; f++)
	{
		if (rand() % 2)
		{
			LS_InsertInFactory(sol[f], factft[f], totaltft, perEij[f]);
			LS_SwapInFactory(sol[f], factft[f], totaltft, perEij[f]);
		}
		else 
		{
			LS_SwapInFactory(sol[f], factft[f], totaltft, perEij[f]);
			LS_InsertInFactory(sol[f], factft[f], totaltft, perEij[f]);
		}
	}
}

void IGR::Reassignment(vector<vector<int>>& sol, vector<int>& factft, int& totaltft, vector<vector<vector<int>>>& perEij)
{
	//detect the factory that has the maximum span
	int extfac = 0;
	int maxspan = INT_MAX;
	for (int f = 0; f < m_Factories; f++) {
		if (perEij[f][sol[f].size() - 1][m_Machines - 1] < maxspan)
		{
			extfac = f;
			maxspan = perEij[f][sol[f].size() - 1][m_Machines - 1];
		}
	}

	int njobs_extfac = sol[extfac].size();
	vector<int> RefSeq(njobs_extfac);
	GetOrderSeq(njobs_extfac, RefSeq);
	random_shuffle(RefSeq.begin(), RefSeq.end());
	int limitjobs = min(njobs_extfac, m_limp);
	int Cnt = 0;
	while (Cnt < limitjobs)
	{
		int pos = RefSeq[Cnt];
		bool ImpFlag = this->CarryBestInsertByMinGlobalTFT_randfactory(extfac, pos, sol, factft, totaltft, perEij);
		if (ImpFlag)
		{
			//detect the factory that has the maximum span
			extfac = 0;
			maxspan = INT_MAX;
			for (int f = 0; f < m_Factories; f++) {
				if (perEij[f][sol[f].size() - 1][m_Machines - 1] < maxspan)
				{
					extfac = f;
					maxspan = perEij[f][sol[f].size() - 1][m_Machines - 1];
				}
			}

			njobs_extfac = sol[extfac].size();
			GetOrderSeq(njobs_extfac, RefSeq);
			random_shuffle(RefSeq.begin(), RefSeq.end());
			limitjobs = min(njobs_extfac, m_limp);
			Cnt = 0;

			if (Base::GetElapsedProcessTime() - GetStartTime() >= m_CPUtime)
				break;
		}
		else
		{
			Cnt++;
		}
	}
}

void IGR::Permutation(vector<vector<int>>& sol, vector<int>& factft, int& totaltft, vector<vector<vector<int>>>& perEij)
{
	//detect the factory that has the maximum span
	int extfac = 0;
	int maxspan = INT_MAX;
	for (int f = 0; f < m_Factories; f++) {
		if (perEij[f][sol[f].size() - 1][m_Machines - 1] < maxspan)
		{
			extfac = f;
			maxspan = perEij[f][sol[f].size() - 1][m_Machines - 1];
		}
	}

	int njobs_extfac = sol[extfac].size();
	vector<int> RefSeq(njobs_extfac);
	GetOrderSeq(njobs_extfac, RefSeq);
	random_shuffle(RefSeq.begin(), RefSeq.end());
	int limitjobs = min(njobs_extfac, m_limp);
	int Cnt = 0;
	while (Cnt < limitjobs)
	{
		int pos = RefSeq[Cnt];
		bool ImpFlag = this->CarryBestPermutationByMinGlobalTFT_randfactory(m_limp, extfac, pos, sol, factft, totaltft, perEij);
		if (ImpFlag)
		{
			//detect the factory that has the maximum span
			extfac = 0;
			maxspan = INT_MAX;
			for (int f = 0; f < m_Factories; f++) {
				if (perEij[f][sol[f].size() - 1][m_Machines - 1] < maxspan)
				{
					extfac = f;
					maxspan = perEij[f][sol[f].size() - 1][m_Machines - 1];
				}
			}

			njobs_extfac = sol[extfac].size();
			GetOrderSeq(njobs_extfac, RefSeq);
			random_shuffle(RefSeq.begin(), RefSeq.end());
			limitjobs = min(njobs_extfac, m_limp);
			Cnt = 0;

			if (Base::GetElapsedProcessTime() - GetStartTime() >= m_CPUtime)
				break;
		}
		else
		{
			Cnt++;
		}
	}
}

