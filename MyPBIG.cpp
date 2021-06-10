#include "MyPBIG.h"

MyPBIG::MyPBIG()
{
	m_Evolution_Type = -1;
	m_lamda = -1;
}

void MyPBIG::SetParameters_PBIG(int num_pop, long CPUtime,int d)
{
	m_num_pop = num_pop;
	m_CPUtime = CPUtime * m_Jobs * m_Machines;
	m_elapCPUtime = m_CPUtime;
	m_d = d;

	m_flag_improve = false;
	m_popu_sol.resize(m_num_pop);
	m_popu_factft.resize(m_num_pop);
	m_popu_totaltft.resize(m_num_pop);
	m_popu_perEij.resize(m_num_pop);
	for (int p = 0; p < m_num_pop; p++)
		m_popu_perEij[p].resize(m_Factories);
	GetOrderSeq(m_Jobs, m_FullSeq);
}

void MyPBIG::SetParameters2_PBIG(int num_pop, long CPUtime, int d, float lamda)
{
	m_num_pop = num_pop;
	m_CPUtime = CPUtime * m_Jobs * m_Machines;
	m_elapCPUtime = m_CPUtime;
	m_d = d;
	m_lamda = lamda;

	m_flag_improve = false;
	m_popu_sol.resize(m_num_pop);
	m_popu_factft.resize(m_num_pop);
	m_popu_totaltft.resize(m_num_pop);
	m_popu_perEij.resize(m_num_pop);
	for (int p = 0; p < m_num_pop; p++)
		m_popu_perEij[p].resize(m_Factories);
	GetOrderSeq(m_Jobs, m_FullSeq);
}

void MyPBIG::SetParameters3_PBIG(int num_pop, int Type_Evolution, long CPUtime, int d, float lamda)
{
	m_num_pop = num_pop;
	m_CPUtime = CPUtime * m_Jobs * m_Machines;
	m_elapCPUtime = m_CPUtime;
	m_d = d;
	m_lamda = lamda;
	m_Evolution_Type = Type_Evolution;

	m_flag_improve = false;
	m_popu_sol.resize(m_num_pop);
	m_popu_factft.resize(m_num_pop);
	m_popu_totaltft.resize(m_num_pop);
	m_popu_perEij.resize(m_num_pop);
	for (int p = 0; p < m_num_pop; p++)
		m_popu_perEij[p].resize(m_Factories);
	GetOrderSeq(m_Jobs, m_FullSeq);
}

int MyPBIG::MyPBIG_Run(string typeInitial)
{
	SetStartTime(Base::GetElapsedProcessTime());
	this->GetTotalPTime();
	this->GetTotalLoad();
	int best_pop = -1;
	if (typeInitial == "DWPFE")
		best_pop = InitialPopulationDPFE();
	else if (typeInitial == "DLR")
		best_pop = InitialPopulationDLR();
	else if (typeInitial == "DNEH")
		best_pop = InitialPopulationDNEH();

	vector<vector<int>> evo_sol_1(m_popu_sol[best_pop]);
	vector<int> evo_factft_1(m_popu_factft[best_pop]);
	int evo_totaltft_1(m_popu_totaltft[best_pop]);
	vector<vector<vector<int>>> evo_perEij_1(m_popu_perEij[best_pop]);

	LS_HybridAmongAllFactories_Reference(m_BestSol_Seq, m_lamda, evo_sol_1, evo_factft_1, evo_totaltft_1, evo_perEij_1);
	if (evo_totaltft_1 < m_BestSol_TotalTFT)
	{
		m_BestSol_Seq = evo_sol_1;
		m_BestSol_FacTFT = evo_factft_1;
		m_BestSol_TotalTFT = evo_totaltft_1;
		m_flag_improve = true;
	}
	PopulationUpdate3(evo_sol_1, evo_factft_1, evo_totaltft_1, evo_perEij_1);

	while (Base::GetElapsedProcessTime() - GetStartTime() < m_CPUtime)
	{
		vector<vector<int>> evo_sol;
		vector<int> evo_factft;
		int evo_totaltft = 0;
		vector<vector<vector<int>>> evo_perEij;

		EvaluationDR3(evo_sol, evo_factft, evo_totaltft, evo_perEij);

		if (evo_totaltft < m_BestSol_TotalTFT)
		{
			m_BestSol_Seq = evo_sol;
			m_BestSol_FacTFT = evo_factft;
			m_BestSol_TotalTFT = evo_totaltft;
			m_flag_improve = true;
		}

		LS_HybridAmongAllFactories_Reference(m_BestSol_Seq, m_lamda, evo_sol, evo_factft, evo_totaltft, evo_perEij);

		if (evo_totaltft < m_BestSol_TotalTFT)
		{
			m_BestSol_Seq = evo_sol;
			m_BestSol_FacTFT = evo_factft;
			m_BestSol_TotalTFT = evo_totaltft;
			m_flag_improve = true;
		}

		PopulationUpdate3(evo_sol, evo_factft, evo_totaltft, evo_perEij);
	}
	cout << (Base::GetElapsedProcessTime() - GetStartTime()) << ends;
	Check(m_BestSol_Seq, m_BestSol_FacTFT, m_BestSol_TotalTFT);
	return m_BestSol_TotalTFT;
}

int MyPBIG::InitialPopulationDLR()
{
	m_popu_totaltft[0] = DLR(m_popu_sol[0], m_popu_factft[0], m_popu_perEij[0]);
	for (int p = 1; p < m_num_pop; p++)
	{
		m_popu_totaltft[p] = RandDLR(m_popu_sol[p], m_popu_factft[p], m_popu_perEij[p]);
	}
	int best_pop = min_element(m_popu_totaltft.begin(), m_popu_totaltft.end()) - m_popu_totaltft.begin();
	m_BestSol_Seq.assign(m_popu_sol[best_pop].begin(), m_popu_sol[best_pop].end());
	m_BestSol_FacTFT.assign(m_popu_factft[best_pop].begin(), m_popu_factft[best_pop].end());
	m_BestSol_TotalTFT = m_popu_totaltft[best_pop];
	return best_pop;
}

int MyPBIG::InitialPopulationDPFE()
{
	m_popu_totaltft[0] = DPFE(m_popu_sol[0], m_popu_factft[0], m_popu_perEij[0]);
	for (int p = 1; p < m_num_pop; p++)
	{
		m_popu_totaltft[p] = RandDPFE(m_popu_sol[p], m_popu_factft[p], m_popu_perEij[p]);
	}

	int best_pop = min_element(m_popu_totaltft.begin(), m_popu_totaltft.end()) - m_popu_totaltft.begin();
	m_BestSol_Seq.assign(m_popu_sol[best_pop].begin(), m_popu_sol[best_pop].end());
	m_BestSol_FacTFT.assign(m_popu_factft[best_pop].begin(), m_popu_factft[best_pop].end());
	m_BestSol_TotalTFT = m_popu_totaltft[best_pop];
	return best_pop;
} 

int MyPBIG::InitialPopulationDNEH()
{
	m_popu_totaltft[0] = DNEH(m_popu_sol[0], m_popu_factft[0], m_popu_perEij[0]);
	for (int p = 1; p < m_num_pop; p++)
	{
		m_popu_totaltft[p] = RandDNEH(m_popu_sol[p], m_popu_factft[p], m_popu_perEij[p]);
	}
	int best_pop = min_element(m_popu_totaltft.begin(), m_popu_totaltft.end()) - m_popu_totaltft.begin();
	m_BestSol_Seq.assign(m_popu_sol[best_pop].begin(), m_popu_sol[best_pop].end());
	m_BestSol_FacTFT.assign(m_popu_factft[best_pop].begin(), m_popu_factft[best_pop].end());
	m_BestSol_TotalTFT = m_popu_totaltft[best_pop];
	return best_pop;
}


void MyPBIG::PopulationUpdate3(const vector<vector<int>>& NewSol, const vector<int>& NewFacTFT, int NewTotalTFT,
	const vector<vector<vector<int>>>& perEij)
{
	if (m_flag_improve)//已知最优解有改进
	{
		m_popu_sol.clear();
		m_popu_factft.clear();
		m_popu_totaltft.clear();
		m_popu_perEij.clear();

		m_popu_sol.push_back(NewSol);
		m_popu_factft.push_back(NewFacTFT);
		m_popu_totaltft.push_back(NewTotalTFT);
		m_popu_perEij.push_back(perEij);

		m_flag_improve = false;
	}
	else
	{
		bool flag = JudgeSameFromPopulation(NewSol, NewFacTFT, NewTotalTFT);
		if (flag)
		{
			if (m_popu_totaltft.size() < m_num_pop)
			{
				m_popu_sol.push_back(NewSol);
				m_popu_factft.push_back(NewFacTFT);
				m_popu_totaltft.push_back(NewTotalTFT);
				m_popu_perEij.push_back(perEij);
			}
			else
			{
				int worst_pop = max_element(m_popu_totaltft.begin(), m_popu_totaltft.end()) - m_popu_totaltft.begin();
				if (NewTotalTFT <= m_popu_totaltft[worst_pop])
				{
					m_popu_sol[worst_pop] = NewSol;
					m_popu_factft[worst_pop] = NewFacTFT;
					m_popu_totaltft[worst_pop] = NewTotalTFT;
					m_popu_perEij[worst_pop] = perEij;
				}
			}
		}
	}
}

bool MyPBIG::JudgeSameFromPopulation(const vector<vector<int>>& Sol, const vector<int>& FacTFT, int TotalTFT)
{
	bool Flag = true;
	for (int p = 0; p < m_popu_totaltft.size(); p++)
	{
		Flag = JudgeSameForTwoSol(Sol, FacTFT, TotalTFT, m_popu_sol[p], m_popu_factft[p], m_popu_totaltft[p]);
		if (Flag == false)
			return false;
	}
	return true;
}

bool MyPBIG::Judgesame()
{
	for(int p=0;p<m_popu_totaltft.size();p++)
		for (int t = p + 1; t < m_popu_totaltft.size(); t++)
		{
			bool flag = JudgeSameForTwoSol(m_popu_sol[p], m_popu_factft[p], m_popu_totaltft[p],
				m_popu_sol[t], m_popu_factft[t], m_popu_totaltft[t]);
			if (!flag)
				return false;
		}
	return true;
}

void MyPBIG::ChooseExtJob(vector<int>& Des_Seq)
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

//随机选择的一个个体分解重构,evo---随机选择的
void MyPBIG::EvaluationDR1(vector<vector<int>>& evo_sol, vector<int>& evo_factft, int& evo_totaltft,
	vector<vector<vector<int>>>& perEij)
{
	vector<int> Des_Seq;
	ChooseExtJob(Des_Seq);
	evo_totaltft = Destruction_Known_PBIG(Des_Seq, evo_sol, evo_factft, perEij);
	evo_totaltft = Reconstruction_PBIG(Des_Seq, evo_sol, evo_factft, perEij);
}

//每个个体均分解相同的若干工件、重构，选择totaltft最小的新解
void MyPBIG::EvaluationDR2(vector<vector<int>>& evo_sol, vector<int>& evo_factft, int& evo_totaltft,
	vector<vector<vector<int>>>& perEij)
{
	vector<int> Des_Seq;
	ChooseExtJob(Des_Seq);
	int cursize_pop = m_popu_totaltft.size();
	int best_totaltft = INT_MAX;
	for (int p = 0; p < cursize_pop; p++)
	{
		vector<vector<int>> test_sol(m_popu_sol[p]);
		vector<int> test_factft(m_popu_factft[p]);
		int test_totaltft(m_popu_totaltft[p]);
		vector<vector<vector<int>>> test_perEij(m_popu_perEij[p]);

		vector<int> CopyDes_Seq(Des_Seq);
		test_totaltft = Destruction_Known_PBIG(CopyDes_Seq, test_sol, test_factft, test_perEij);
		test_totaltft = Reconstruction_PBIG(CopyDes_Seq, test_sol, test_factft, test_perEij);
		if (test_totaltft < best_totaltft)
		{
			evo_sol = test_sol;
			evo_factft = test_factft;
			evo_totaltft = test_totaltft;
			perEij = test_perEij;
			best_totaltft = test_totaltft;
		}
	}
}

//每个个体均抽出相同的若干工件，选择totaltft最小的部分解重构
void MyPBIG::EvaluationDR3(vector<vector<int>>& evo_sol, vector<int>& evo_factft, int& evo_totaltft,
	vector<vector<vector<int>>>& perEij)
{
	vector<int> Des_Seq;
	ChooseExtJob(Des_Seq);
	int cursize_pop = m_popu_totaltft.size();
	int partial_totaltft = INT_MAX;
	for (int p = 0; p < cursize_pop; p++)
	{
		vector<vector<int>> test_sol(m_popu_sol[p]);
		vector<int> test_factft(m_popu_factft[p]);
		int test_totaltft(m_popu_totaltft[p]);
		vector<vector<vector<int>>> test_perEij(m_popu_perEij[p]);

		vector<int> CopyDes_Seq(Des_Seq);
		test_totaltft = Destruction_Known_PBIG(CopyDes_Seq, test_sol, test_factft, test_perEij);
		if (test_totaltft < partial_totaltft)
		{
			evo_sol = test_sol;
			evo_factft = test_factft;
			evo_totaltft = test_totaltft;
			perEij = test_perEij;
			partial_totaltft = test_totaltft;
		}
	}
	evo_totaltft = Reconstruction_PBIG(Des_Seq, evo_sol, evo_factft, perEij);
}

//随机抽取出工件
int MyPBIG::Destruction_Unknown_PBIG(vector<int>& Des_Seq, vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij)
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
	return accumulate(FacTFT.begin(), FacTFT.end(), 0);
}

//抽取出已知工件
int MyPBIG::Destruction_Known_PBIG(const vector<int>& Des_Seq, vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij)
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

//重构-----随机
int MyPBIG::Reconstruction_PBIG(vector<int>& Des_Seq, vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij)
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
	return accumulate(FacTFT.begin(), FacTFT.end(), 0);
}
