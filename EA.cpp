#include "EA.h"

EA::EA()
{
}

EA::~EA()
{
}

void EA::SetParameters_EA(int num_pop, long CPUtime)
{
	m_num_pop = num_pop;
	m_CPUtime = CPUtime * m_Jobs * m_Machines;
	m_flag_improve = false;
}

int EA::EA_Run()
{
	SetStartTime(Base::GetElapsedProcessTime());
	this->GetTotalPTime();
	this->GetTotalLoad();
	int best_pop = InitialPopulation();
	int cur_pos_for_rls = 0;
	vector<vector<int>> evo_sol(m_popu_sol[best_pop]);
	vector<int> evo_factft(m_popu_factft[best_pop]);
	int evo_totaltft(m_popu_totaltft[best_pop]);
	vector<vector<vector<int>>> evo_perEij(m_popu_perEij[best_pop]);

	LS_InsertAmongAllFactories_Reference(m_BestSol_Seq, evo_sol, evo_factft, evo_totaltft, evo_perEij, cur_pos_for_rls);
	if (evo_totaltft < m_BestSol_TotalTFT)
		m_flag_improve = true;
	PopulationUpdate(evo_sol, evo_factft, evo_totaltft, evo_perEij);
	int Cnt = 0;
	while (Base::GetElapsedProcessTime() - GetStartTime() < m_CPUtime)
	{
		int choosed_pop = rand() % m_popu_totaltft.size();
		vector<vector<int>> evo_sol(m_popu_sol[choosed_pop]);
		vector<int> evo_factft(m_popu_factft[choosed_pop]);
		int evo_totaltft(m_popu_totaltft[choosed_pop]);
		vector<vector<vector<int>>> evo_perEij(m_popu_perEij[choosed_pop]);
		Mutation(evo_sol, evo_factft, evo_totaltft, evo_perEij);
		LS_InsertAmongAllFactories_Reference(m_BestSol_Seq, evo_sol, evo_factft, evo_totaltft, evo_perEij, cur_pos_for_rls);
		if (evo_totaltft < m_BestSol_TotalTFT)
		{
			m_BestSol_Seq = evo_sol;
			m_BestSol_FacTFT = evo_factft;
			m_BestSol_TotalTFT = evo_totaltft;
			m_flag_improve = true;
		}
		PopulationUpdate(evo_sol, evo_factft, evo_totaltft, evo_perEij);
		Cnt++;
	}
	cout << Cnt << ends << (Base::GetElapsedProcessTime() - GetStartTime()) << ends;
	Check(m_BestSol_Seq, m_BestSol_FacTFT, m_BestSol_TotalTFT);
	return m_BestSol_TotalTFT;
}

int EA::EA_Run2(vector<int>& _result)
{
	SetStartTime(Base::GetElapsedProcessTime());
	this->GetTotalPTime();
	this->GetTotalLoad();
	int best_pop = InitialPopulation();
	int cur_pos_for_rls = 0;
	vector<vector<int>> evo_sol(m_popu_sol[best_pop]);
	vector<int> evo_factft(m_popu_factft[best_pop]);
	int evo_totaltft(m_popu_totaltft[best_pop]);
	vector<vector<vector<int>>> evo_perEij(m_popu_perEij[best_pop]);

	LS_InsertAmongAllFactories_Reference(m_BestSol_Seq, evo_sol, evo_factft, evo_totaltft, evo_perEij, cur_pos_for_rls);
	if (evo_totaltft < m_BestSol_TotalTFT)
		m_flag_improve = true;
	PopulationUpdate(evo_sol, evo_factft, evo_totaltft, evo_perEij);

	long elapCPUtime = m_Jobs * m_Machines * 10;
	while (Base::GetElapsedProcessTime() - GetStartTime() < elapCPUtime)
	{
		int choosed_pop = rand() % m_popu_totaltft.size();
		vector<vector<int>> evo_sol(m_popu_sol[choosed_pop]);
		vector<int> evo_factft(m_popu_factft[choosed_pop]);
		int evo_totaltft(m_popu_totaltft[choosed_pop]);
		vector<vector<vector<int>>> evo_perEij(m_popu_perEij[choosed_pop]);
		Mutation(evo_sol, evo_factft, evo_totaltft, evo_perEij);
		LS_InsertAmongAllFactories_Reference(m_BestSol_Seq, evo_sol, evo_factft, evo_totaltft, evo_perEij, cur_pos_for_rls);
		if (evo_totaltft < m_BestSol_TotalTFT)
		{
			m_BestSol_Seq = evo_sol;
			m_BestSol_FacTFT = evo_factft;
			m_BestSol_TotalTFT = evo_totaltft;
			m_flag_improve = true;
		}
		PopulationUpdate(evo_sol, evo_factft, evo_totaltft, evo_perEij);
	}
	_result[0] = m_BestSol_TotalTFT;

	elapCPUtime = m_Jobs * m_Machines * 20;
	while (Base::GetElapsedProcessTime() - GetStartTime() < elapCPUtime)
	{
		int choosed_pop = rand() % m_popu_totaltft.size();
		vector<vector<int>> evo_sol(m_popu_sol[choosed_pop]);
		vector<int> evo_factft(m_popu_factft[choosed_pop]);
		int evo_totaltft(m_popu_totaltft[choosed_pop]);
		vector<vector<vector<int>>> evo_perEij(m_popu_perEij[choosed_pop]);
		Mutation(evo_sol, evo_factft, evo_totaltft, evo_perEij);
		LS_InsertAmongAllFactories_Reference(m_BestSol_Seq, evo_sol, evo_factft, evo_totaltft, evo_perEij, cur_pos_for_rls);
		if (evo_totaltft < m_BestSol_TotalTFT)
		{
			m_BestSol_Seq = evo_sol;
			m_BestSol_FacTFT = evo_factft;
			m_BestSol_TotalTFT = evo_totaltft;
			m_flag_improve = true;
		}
		PopulationUpdate(evo_sol, evo_factft, evo_totaltft, evo_perEij);
	}
	_result[1] = m_BestSol_TotalTFT;

	elapCPUtime = m_Jobs * m_Machines * 40;
	while (Base::GetElapsedProcessTime() - GetStartTime() < elapCPUtime)
	{
		int choosed_pop = rand() % m_popu_totaltft.size();
		vector<vector<int>> evo_sol(m_popu_sol[choosed_pop]);
		vector<int> evo_factft(m_popu_factft[choosed_pop]);
		int evo_totaltft(m_popu_totaltft[choosed_pop]);
		vector<vector<vector<int>>> evo_perEij(m_popu_perEij[choosed_pop]);
		Mutation(evo_sol, evo_factft, evo_totaltft, evo_perEij);
		LS_InsertAmongAllFactories_Reference(m_BestSol_Seq, evo_sol, evo_factft, evo_totaltft, evo_perEij, cur_pos_for_rls);
		if (evo_totaltft < m_BestSol_TotalTFT)
		{
			m_BestSol_Seq = evo_sol;
			m_BestSol_FacTFT = evo_factft;
			m_BestSol_TotalTFT = evo_totaltft;
			m_flag_improve = true;
		}
		PopulationUpdate(evo_sol, evo_factft, evo_totaltft, evo_perEij);
	}
	_result[2] = m_BestSol_TotalTFT;

	elapCPUtime = m_Jobs * m_Machines * 80;
	while (Base::GetElapsedProcessTime() - GetStartTime() < elapCPUtime)
	{
		int choosed_pop = rand() % m_popu_totaltft.size();
		vector<vector<int>> evo_sol(m_popu_sol[choosed_pop]);
		vector<int> evo_factft(m_popu_factft[choosed_pop]);
		int evo_totaltft(m_popu_totaltft[choosed_pop]);
		vector<vector<vector<int>>> evo_perEij(m_popu_perEij[choosed_pop]);
		Mutation(evo_sol, evo_factft, evo_totaltft, evo_perEij);
		LS_InsertAmongAllFactories_Reference(m_BestSol_Seq, evo_sol, evo_factft, evo_totaltft, evo_perEij, cur_pos_for_rls);
		if (evo_totaltft < m_BestSol_TotalTFT)
		{
			m_BestSol_Seq = evo_sol;
			m_BestSol_FacTFT = evo_factft;
			m_BestSol_TotalTFT = evo_totaltft;
			m_flag_improve = true;
		}
		PopulationUpdate(evo_sol, evo_factft, evo_totaltft, evo_perEij);
	}
	_result[3] = m_BestSol_TotalTFT;

	cout << (Base::GetElapsedProcessTime() - GetStartTime()) << ends;
	Check(m_BestSol_Seq, m_BestSol_FacTFT, m_BestSol_TotalTFT);
	return m_BestSol_TotalTFT;
}

int EA::InitialPopulation()
{
	m_popu_sol.resize(m_num_pop);
	m_popu_factft.resize(m_num_pop);
	m_popu_totaltft.resize(m_num_pop);
	m_popu_perEij.resize(m_num_pop);

	for (int p = 0; p < m_num_pop; p++)
		m_popu_perEij[p].resize(m_Factories);

	GetOrderSeq(m_Jobs, m_FullSeq);
	m_popu_totaltft[0] = NEHR2A4(m_popu_sol[0], m_popu_factft[0], m_popu_perEij[0]);
	for (int p = 1; p < m_num_pop; p++)
	{
		m_popu_totaltft[p] = RandNEHR2A4(m_popu_sol[p], m_popu_factft[p], m_popu_perEij[p]);
	}
	int best_pop = min_element(m_popu_totaltft.begin(), m_popu_totaltft.end()) - m_popu_totaltft.begin();
	m_BestSol_Seq.assign(m_popu_sol[best_pop].begin(), m_popu_sol[best_pop].end());
	m_BestSol_FacTFT.assign(m_popu_factft[best_pop].begin(), m_popu_factft[best_pop].end());
	m_BestSol_TotalTFT = m_popu_totaltft[best_pop];
	return best_pop;
}

void EA::PopulationUpdate(const vector<vector<int>>& NewSol, const vector<int>& NewFacTFT, int NewTotalTFT,
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
		m_BestSol_Seq = NewSol;
		m_BestSol_FacTFT = NewFacTFT;
		m_BestSol_TotalTFT = NewTotalTFT;
	}
	else
	{
		bool flag = JudgeSameFromPopulation(NewSol, NewFacTFT, NewTotalTFT);
		if (flag)
		{
			m_popu_sol.push_back(NewSol);
			m_popu_factft.push_back(NewFacTFT);
			m_popu_totaltft.push_back(NewTotalTFT);
			m_popu_perEij.push_back(perEij);

			if (m_popu_totaltft.size() > m_num_pop)//超出容量
			{
				int worst_pop = max_element(m_popu_totaltft.begin(), m_popu_totaltft.end()) - m_popu_totaltft.begin();
				m_popu_sol.erase(m_popu_sol.begin() + worst_pop);
				m_popu_factft.erase(m_popu_factft.begin() + worst_pop);
				m_popu_totaltft.erase(m_popu_totaltft.begin() + worst_pop);
				m_popu_perEij.erase(m_popu_perEij.begin() + worst_pop);
			}
		}
	}
}

bool EA::JudgeSameFromPopulation(const vector<vector<int>>& Sol, const vector<int>& FacTFT, int TotalTFT)
{
	bool Flag = true;
	for (int p = 0; p < m_popu_totaltft.size(); p++)
	{
		Flag = JudgeSameForTwoSol(Sol, FacTFT, TotalTFT, m_popu_sol[p], m_popu_factft[p], m_popu_totaltft[p]);
		if (!Flag)
			return Flag;
	}
	return Flag;
}

void EA::Mutation(vector<vector<int>>& evo_sol, vector<int>& evo_factft, int& evo_totaltft, 
	vector<vector<vector<int>>> &perEij)
{
	int fac1 = -1, fac2 = -1, pos1 = -1, pos2 = -1;
	do
	{
		fac1 = rand() % m_Factories;
		fac2 = rand() % m_Factories;
		pos1 = rand() % evo_sol[fac1].size();
		pos2 = rand() % evo_sol[fac2].size();
	} while (fac1 == fac2 && pos1 == pos2);

	if (fac1 == fac2)
		SwapTwoJobsInsideFactory(pos1, pos2, evo_sol[fac1], perEij[fac1], evo_factft[fac1], evo_totaltft);
	else
		SwapTwoJobsBetweenFactories(fac1, fac2, pos1, pos2, evo_sol, perEij, evo_factft, evo_totaltft);
}

