#include "HIG.h"

int HIG::HIG_Run()
{
	long StartTime = Base::GetElapsedProcessTime();
	this->GetTotalPTime();
	this->GetTotalLoad();
	this->GetOrderSeq(this->m_Jobs, this->m_FullSeq);

	this->m_CurTotalTFT = PWNEH(this->m_CurSol, m_CurFacTFT, this->m_CurFacSpan, m_CurperEij);

	UpdateBestSol(m_CurSol, m_CurFacTFT, m_CurTotalTFT);

	m_Temperature = this->m_TotalLoad * this->m_Value_T / (m_Jobs * m_Machines * 10);
	int Cnt_perTemp = 0;//每个温度下迭代的次数

	while (Base::GetElapsedProcessTime() - StartTime < m_limittime)
	{
		vector<vector<int>> Evo_Sol(m_CurSol);
		vector<int> Evo_FacTFT(m_CurFacTFT);
		int Evo_TotalTFT = m_CurTotalTFT;
		vector<int>Evo_FacSpan(m_CurFacSpan);
		vector<vector<vector<int>>> Evo_perEij(m_CurperEij);

		//分解重构产生候选解
		vector<int> Des_Seq;
		Destruction(Evo_Sol, Des_Seq, Evo_FacTFT, Evo_FacSpan, Evo_perEij);
		Evo_TotalTFT = Reconstruction(Des_Seq, Evo_Sol, Evo_FacTFT, Evo_FacSpan, Evo_perEij);

		UpdateSol(Evo_Sol, Evo_FacTFT, Evo_FacSpan, Evo_TotalTFT, Evo_perEij);

		//更新温度
		Cnt_perTemp++;
		if (Cnt_perTemp % this->m_Iters_perT == 0)
		{
			m_Temperature = m_Temperature * this->m_Value_Cooling;
			Cnt_perTemp = 0;
		}
	}

	long EndTime = Base::GetElapsedProcessTime();
	cout << m_Temperature << ends << EndTime - StartTime << ends;
	this->Check(this->m_BestSol_Seq, this->m_BestSol_FacTFT, this->m_BestSol_TotalTFT);
	return this->m_BestSol_TotalTFT;
}

void HIG::SetParameters(float value_T, int Iters_perT, float Value_Cooling, float TLTlow, float TLThigh, int Minnum_Des, int Maxnum_Des, long time_fator)
{
	this->m_Value_T = value_T;
	this->m_Iters_perT = Iters_perT;
	this->m_Value_Cooling = Value_Cooling;
	this->m_TLTlow = TLTlow;
	this->m_TLThigh = TLThigh;
	this->m_Minnum_Des = Minnum_Des;
	this->m_Maxnum_Des = Maxnum_Des;

	this->m_limittime = m_Jobs * m_Machines * time_fator;
	
	m_CurperEij.resize(m_Factories);

	while (!this->m_queueTL.empty())
		this->m_queueTL.pop();
	this->m_JobFlag.resize(this->m_Jobs, false);
}

void HIG::Destruction(vector<vector<int>>& Sol, vector<int>& Des_Seq, vector<int>& FacTFT, vector<int>& FacSpan, 
	vector<vector<vector<int>>>& perEij)
{
	float u = (float)rand() / RAND_MAX;
	int TLT = this->m_TLTlow * this->m_Jobs + u * (this->m_TLThigh - this->m_TLTlow) * this->m_Jobs;
	while (this->m_queueTL.size() > TLT)//如果禁忌表容量超过TLT，则把前面的工件去掉
	{
		int Job = this->m_queueTL.front();
		this->m_JobFlag[Job] = false;
		this->m_queueTL.pop();
	}

	vector<bool> factoryflag_change(m_Factories, false);
	int curnum_extjobs = 0;

	vector<int> order_factories(m_Factories);
	IncreasedOrder(m_Factories, FacSpan, order_factories);
	//Extract a job from each factory with the maximum makespan
	int curfac = m_Factories - 1;
	do
	{
		factoryflag_change[order_factories[curfac]] = true;
		curnum_extjobs++;
		bool flag = RemoveOneJob(Sol[order_factories[curfac]], Des_Seq, TLT);
		curfac--;

	} while (curfac >= 0 && FacSpan[order_factories[curfac]] == FacSpan[order_factories[curfac + 1]]);

	//the case that all factories have the same makespan is not detailed in the literature,we extract a job from every factory here
	//Extract a job from each factory with the minimum makespan
	if (curfac != -1)
	{
		curfac = 0;
		do
		{
			factoryflag_change[order_factories[curfac]] = true;
			curnum_extjobs++;
			bool flag = RemoveOneJob(Sol[order_factories[curfac]], Des_Seq, TLT);
			curfac++;
		} while (curfac < m_Factories && FacSpan[order_factories[curfac]] == FacSpan[order_factories[curfac - 1]]);
	}

	int num_Des = rand() % (this->m_Maxnum_Des - this->m_Minnum_Des + 1) + this->m_Minnum_Des;
	//the remained （num_Des-curnum_extjobs） are extracted randomly
	while (curnum_extjobs < num_Des)
	{
		int factory = -1;
		do
		{
			factory = rand() % m_Factories;
		} while (Sol[factory].empty());
		bool flag = RemoveOneJob(Sol[factory], Des_Seq, TLT);
		if (flag)
		{
			factoryflag_change[factory] = true;
			curnum_extjobs++;
		}
	}

	for (int f = 0; f < m_Factories; f++)
		if (factoryflag_change[f])
			GetEij(Sol[f], perEij[f]);
}

int HIG::Reconstruction(const vector<int>& Des_Seq, vector<vector<int>>& Sol, vector<int>& FacTFT, vector<int>& FacSpan,
	vector<vector<vector<int>>>& perEij)
{
	unsigned int Size_DesSeq = Des_Seq.size();
	for (unsigned int j = 0; j < Size_DesSeq; j++)
	{
		int Job = Des_Seq[j];
		int BestPos = -1, BestFac = -1;
		int Min_TFT = FindBestInsertPosFacByMinIncreasedTFT(Sol, FacTFT, perEij, Job, BestFac, BestPos);
		FacTFT[BestFac] = Min_TFT;
		InsertAJobToSeq(Job, BestPos, Sol[BestFac], perEij[BestFac], FacTFT[BestFac]);
		FacSpan[BestFac] = perEij[BestFac][Sol[BestFac].size() - 1][m_Machines - 1];
	}
	return accumulate(FacTFT.begin(), FacTFT.end(), 0);
}


bool HIG::TestExsitJobnotinTL(vector<int> Seq)
{
	int Size_Seq = Seq.size();
	for (int j = 0; j < Size_Seq; j++)
	{
		if (!this->m_JobFlag[Seq[j]])
		{
			return true;
		}
	}
	return false;
}

bool HIG::RemoveOneJob(vector<int>& Seq, vector<int>& Des_Seq, int TLT)
{
	bool Flag = TestExsitJobnotinTL(Seq);
	if (Flag && !Seq.empty())
	{
		int extpos = -1;
		do
		{
			extpos = rand() % Seq.size();
		} while (this->m_JobFlag[Seq[extpos]]);

		int extjob = Seq[extpos];
		Des_Seq.push_back(extjob);

		if (this->m_queueTL.size() == TLT)
		{
			int ChangeJob = this->m_queueTL.front();
			this->m_JobFlag[ChangeJob] = false;
			this->m_queueTL.pop();
			this->m_queueTL.push(extjob);
			this->m_JobFlag[extjob] = true;
		}
		else
		{
			this->m_queueTL.push(extjob);
			this->m_JobFlag[extjob] = true;
		}

		Seq.erase(Seq.begin() + extpos);
		return true;
	}
	else
		return false;
}

void HIG::UpdateBestSol(const vector<vector<int>>& Sol, const vector<int>& FacTFT, int TotalTFT)
{
	if (TotalTFT < m_BestSol_TotalTFT)
	{
		m_BestSol_Seq = Sol;
		m_BestSol_FacTFT = FacTFT;
		m_BestSol_TotalTFT = TotalTFT;
	}
}

void HIG::UpdateSol(const vector<vector<int>>& Sol, const vector<int>& FacTFT, const vector<int>& FacSpan, int TotalTFT, const vector<vector<vector<int>>>& perEij)
{
	//更新解
	if (TotalTFT <= this->m_CurTotalTFT)
	{
		m_CurSol = Sol;
		m_CurFacTFT = FacTFT;
		m_CurFacSpan = FacSpan;
		m_CurTotalTFT = TotalTFT;
		m_CurperEij = perEij;
		UpdateBestSol(Sol, FacTFT, TotalTFT);
	}
	else
	{
		float RandNumber = (float)rand() / RAND_MAX;
		float Index = 1.00 * exp((m_CurTotalTFT - TotalTFT) / m_Temperature);
		if (RandNumber < Index)
		{
			m_CurSol = Sol;
			m_CurFacTFT = FacTFT;
			m_CurFacSpan = FacSpan;
			m_CurTotalTFT = TotalTFT;
			m_CurperEij = perEij;
		}
	}
}

