#include "IGAndILS.h"
#include "Base.h"

void IGAndILS::SetParameters_IG(long C_of_CPUtime, int Setd, float SetT,int Type_Initial)
{
	this->m_CPUtime = C_of_CPUtime * this->m_Jobs * this->m_Machines;
	this->m_d = Setd;
	this->m_T = SetT;
	this->m_Type_initial = Type_Initial;
	this->m_elapCPUtime = this->m_CPUtime;

	this->m_BestSol_Seq.resize(this->m_Factories);
	this->m_BestSol_FacTFT.resize(this->m_Factories);
	this->m_BestSol_TotalTFT = -1;
	this->m_perEij.resize(this->m_Factories);
}

void IGAndILS::SetParametersSecond_ILS(long C_of_CPUtime, int SetPi, int Setd, float SetT, int Type_Initial)
{
	this->m_CPUtime = C_of_CPUtime * this->m_Jobs * this->m_Machines;
	this->m_d = Setd;
	this->m_T = SetT;
	this->m_Pi = SetPi;
	this->m_Type_initial = Type_Initial;
	this->m_elapCPUtime = this->m_CPUtime;

	this->m_BestSol_Seq.resize(this->m_Factories);
	this->m_BestSol_FacTFT.resize(this->m_Factories);
	this->m_BestSol_TotalTFT = -1;
	this->m_perEij.resize(this->m_Factories);
}

int IGAndILS::Initialization(int Type_Init, vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij)
{
	int TotalTFT = 0;
	switch (Type_Init)
	{
	case 0:TotalTFT = DLR(Sol, FacTFT, perEij); break;
	case 1:TotalTFT = DLR_DNEH(0.1, Sol, FacTFT, perEij); break;
	case 2:TotalTFT = DLR_DNEH(0.2, Sol, FacTFT, perEij); break;
	case 3:TotalTFT = DLR_DNEH(0.3, Sol, FacTFT, perEij); break;
	default:
		break;
	}
	return TotalTFT;
}

int IGAndILS::IG_Run()
{
	vector<int> ConvergentArray;
	SetStartTime(Base::GetElapsedProcessTime());
	this->GetTotalPTime();
	this->GetTotalLoad();
	float Temperature = 1.0 * this->m_T * this->m_TotalLoad / (10.0 * this->m_Jobs * this->m_Machines);
	this->m_BestSol_TotalTFT = Initialization(m_Type_initial, this->m_BestSol_Seq, this->m_BestSol_FacTFT, m_perEij);

	LocalSearchType(this->m_BestSol_Seq, this->m_BestSol_FacTFT, this->m_BestSol_TotalTFT, m_perEij);
	vector<vector<int>> CurSol(this->m_BestSol_Seq);
	vector<int> CurFacTFT(this->m_BestSol_FacTFT);
	int CurTotalTFT = this->m_BestSol_TotalTFT;
	int Cnt = 0;
	while (Base::GetElapsedProcessTime() - GetStartTime() < m_CPUtime)//Base::GetElapsedProcessTime()-StarTime<this->m_CPUtime
	{
		vector<vector<int>> Evo_Sol(CurSol);
		vector<int> Evo_FacTFT(CurFacTFT);
		int Evo_TotalTFT = CurTotalTFT;
		vector<vector<vector<int>>> Evo_perEij(this->m_perEij);

		vector<int> Des_Seq;
		Destruction_IG(Evo_Sol, Evo_FacTFT, Des_Seq, Evo_perEij);
		Evo_TotalTFT = Reconstruction_IG(Des_Seq, Evo_Sol, Evo_FacTFT, Evo_perEij);

		LocalSearchType(Evo_Sol, Evo_FacTFT, Evo_TotalTFT, Evo_perEij);
		Cnt++;
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
			float Index = 1.0f * exp((CurTotalTFT - Evo_TotalTFT) / Temperature);
			if (RandNumber < Index)
			{
				CurSol = Evo_Sol;
				CurFacTFT = Evo_FacTFT;
				CurTotalTFT = Evo_TotalTFT;
				m_perEij = Evo_perEij;
			}
		}
	}
	cout << Cnt << ends;
	cout << (Base::GetElapsedProcessTime() - GetStartTime()) << ends;
	this->Check(this->m_BestSol_Seq, this->m_BestSol_FacTFT, this->m_BestSol_TotalTFT);
	return this->m_BestSol_TotalTFT;
}

int IGAndILS::ILS_Run()
{
	SetStartTime(Base::GetElapsedProcessTime());
	this->GetTotalPTime();
	this->GetTotalLoad();
	float Temperature = 1.0 * this->m_T * this->m_TotalLoad / (10.0 * this->m_Jobs * this->m_Machines);
	this->m_BestSol_TotalTFT = Initialization(m_Type_initial, this->m_BestSol_Seq, this->m_BestSol_FacTFT, m_perEij);
	LocalSearchType(this->m_BestSol_Seq, this->m_BestSol_FacTFT, m_BestSol_TotalTFT, m_perEij);

	vector<vector<int>> CurSol(this->m_BestSol_Seq);
	vector<int> CurFacTFT(this->m_BestSol_FacTFT);
	int CurTotalTFT = this->m_BestSol_TotalTFT;
	int Cnt = 0;
	while (Base::GetElapsedProcessTime() - GetStartTime() < this->m_CPUtime)//Base::GetElapsedProcessTime()-StarTime<this->m_CPUtime
	{
		vector<vector<int>> Evo_Sol(CurSol);
		vector<int> Evo_FacTFT(CurFacTFT);
		int Evo_TotalTFT = CurTotalTFT;
		vector<vector<vector<int>>> Evo_perEij(m_perEij);
		Disturbance_ILS(Evo_Sol, Evo_FacTFT, Evo_TotalTFT, Evo_perEij);
		LocalSearchType(Evo_Sol, Evo_FacTFT, Evo_TotalTFT, Evo_perEij);	
		Cnt++;
		//更新解
		if (Evo_TotalTFT <= CurTotalTFT)
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
			float Index = 1.0f * exp((CurTotalTFT - Evo_TotalTFT) / Temperature);
			if (RandNumber < Index)
			{
				CurSol = Evo_Sol;
				CurFacTFT = Evo_FacTFT;
				CurTotalTFT = Evo_TotalTFT;
				m_perEij = Evo_perEij;
			}
		}
	}
	cout << Cnt << ends << (Base::GetElapsedProcessTime() - GetStartTime()) << ends;
	this->Check(this->m_BestSol_Seq, this->m_BestSol_FacTFT, this->m_BestSol_TotalTFT);
	return m_BestSol_TotalTFT;
}

void IGAndILS::LocalSearchType(vector<vector<int>>& Sol, vector<int>& FacTFT, int& TotalTFT, vector<vector<vector<int>>>& perEij)
{
	LS_InsertAmongAllFactories_Reference(m_BestSol_Seq, Sol, FacTFT, TotalTFT, perEij);
}

void IGAndILS::Destruction_IG(vector<vector<int>>& Sol, vector<int>& FacTFT, vector<int>& Des_Seq, vector<vector<vector<int>>>& perEij)
{
	vector<bool> flag_change(m_Factories, false);
	while (Des_Seq.size() < this->m_d)
	{
		int Factory = rand() % this->m_Factories;
		if (Sol[Factory].empty())
			continue;
		int Position = rand() % Sol[Factory].size();
		Des_Seq.push_back(Sol[Factory][Position]);
		Sol[Factory].erase(Sol[Factory].begin() + Position);
		flag_change[Factory] = true;
	}
	for (int f = 0; f < m_Factories; f++)
		if (flag_change[f])
			FacTFT[f] = GetEij(Sol[f], perEij[f]);
}

int IGAndILS::Reconstruction_IG(vector<int>& Des_Seq, vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij)
{
	while (Des_Seq.size())
	{
		int Pos = rand() % Des_Seq.size();
		int Job = Des_Seq[Pos];
		Des_Seq.erase(Des_Seq.begin() + Pos);
		int BestPos = -1, BestFac = -1;
		int Min_TFT = FindBestInsertPosFacByMinIncreasedTFT(Sol, FacTFT, perEij, Job, BestFac, BestPos);
		FacTFT[BestFac] = Min_TFT;
		InsertAJobToSeq(Job, BestPos, Sol[BestFac], perEij[BestFac], FacTFT[BestFac]);
	}
	int totaltft = accumulate(FacTFT.begin(), FacTFT.end(), 0);
	return totaltft;
}

void IGAndILS::Disturbance_ILS(vector<vector<int>>& Sol, vector<int>& FacTFT, int& TotalTFT, vector<vector<vector<int>>>& perEij)
{
	vector<vector<int>> bestdisturb_Sol;
	vector<int> bestdisturb_FacTFT;
	int bestdisturb_TotalTFT = INT_MAX;
	vector<vector<vector<int>>> bestdisturb_perEij;

	for (int w = 0; w < m_Pi; w++)
	{
		vector<vector<int>> CopySol(Sol);
		vector<int> CopySolFacTFT(FacTFT);
		vector<vector<vector<int>>> CopyperEij(perEij);
		int CopyTotalTFT = TotalTFT;

		for (int i = 0; i < m_d; i++)
		{
			//choose a job randomly and insert it into other position
			int ExtFac = -1;
			do
			{
				ExtFac = rand() % this->m_Factories;
			} while (CopySol[ExtFac].empty());
			int ExtPos = rand() % CopySol[ExtFac].size();
			int ExtJob = CopySol[ExtFac][ExtPos];
			int InsertFac = rand() % this->m_Factories;

			if (ExtFac != InsertFac)
			{
				int InsertPos = rand() % (CopySol[InsertFac].size() + 1);
				CopySol[ExtFac].erase(CopySol[ExtFac].begin() + ExtPos);
				CopySol[InsertFac].insert(CopySol[InsertFac].begin() + InsertPos, ExtJob);
			}
			else
			{
				if (CopySol[ExtFac].size() == 1)
				{
					i--;
					continue;
				}
				CopySol[ExtFac].erase(CopySol[ExtFac].begin() + ExtPos);
				int InsertPos = -1;
				do
				{
					InsertPos = rand() % (CopySol[InsertFac].size() + 1);
				} while (InsertPos == ExtPos);
				CopySol[InsertFac].insert(CopySol[InsertFac].begin() + InsertPos, ExtJob);
			}

		}

		int TesttotalTFT = GetTotalFlowtime(CopySol, CopySolFacTFT, CopyperEij);

		if (TesttotalTFT < bestdisturb_TotalTFT)
		{
			bestdisturb_Sol = CopySol;
			bestdisturb_FacTFT = CopySolFacTFT;
			bestdisturb_TotalTFT = TesttotalTFT;
			bestdisturb_perEij = CopyperEij;
			bestdisturb_TotalTFT = TesttotalTFT;
		}
	}
	Sol = bestdisturb_Sol;
	FacTFT = bestdisturb_FacTFT;
	TotalTFT = bestdisturb_TotalTFT;
	perEij = bestdisturb_perEij;
}

