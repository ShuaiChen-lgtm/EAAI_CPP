#include "IG_2S.h"

void IG_2S::SetParameter_IG2S(int S_d1, int S_d2, double S_T, double S_ratio, long time_factor)
{
	this->m_d1 = S_d1;
	this->m_d2 = S_d2;
	this->m_T = S_T;
	this->m_ratio = S_ratio;
	this->m_CPU_time = time_factor * m_Machines * m_Jobs;

	this->m_BestSol_Seq.resize(this->m_Factories);
	this->m_BestSol_FacTFT.resize(this->m_Factories);
	this->m_BestSol_FacSpan.resize(this->m_Factories);
	this->m_BestSol_TotalTFT = -1;
	this->m_perEij.resize(this->m_Factories);
}

int IG_2S::IG2S_Run()
{
	vector<int> ConvergentArray;
	SetStartTime(Base::GetElapsedProcessTime());
	this->GetTotalPTime();
	this->GetTotalLoad();
	float Temperature = 1.0 * this->m_T * this->m_TotalLoad / (10.0f * this->m_Jobs * this->m_Machines);
	this->m_BestSol_TotalTFT = NEH2_en(this->m_BestSol_Seq, this->m_BestSol_FacTFT, this->m_BestSol_FacSpan, m_perEij);

	this->LS3_IG2S(this->m_BestSol_Seq, this->m_BestSol_FacTFT, this->m_BestSol_FacSpan, this->m_BestSol_TotalTFT, m_perEij);

	vector<vector<int>> CurSol(this->m_BestSol_Seq);
	vector<int> CurFacTFT(this->m_BestSol_FacTFT);
	vector<int> CurFacSpan(this->m_BestSol_FacSpan);
	int CurTotalTFT = this->m_BestSol_TotalTFT;

	//------------------the first stage-------------
	while (Base::GetElapsedProcessTime() - GetStartTime() < m_CPU_time * m_ratio)
	{
		vector<vector<int>> Evo_Sol(CurSol);
		vector<int> Evo_FacTFT(CurFacTFT);
		vector<int> Evo_FacSpan(CurFacSpan);
		int Evo_TotalTFT = CurTotalTFT;

		vector<vector<vector<int>>> Evo_perEij(this->m_perEij);
		vector<int> Des_Seq;
		Destruction(Des_Seq, Evo_Sol, Evo_FacSpan, Evo_FacTFT, Evo_perEij);

		Evo_TotalTFT = Reconstruction(Des_Seq, Evo_Sol, Evo_FacSpan, Evo_FacTFT, Evo_perEij);

		LS3_IG2S(Evo_Sol, Evo_FacTFT, Evo_FacSpan, Evo_TotalTFT, Evo_perEij);

		//¸üÐÂ½â
		if (Evo_TotalTFT < CurTotalTFT)
		{
			CurSol = Evo_Sol;
			CurFacTFT = Evo_FacTFT;
			CurTotalTFT = Evo_TotalTFT;
			CurFacSpan = Evo_FacSpan;
			m_perEij = Evo_perEij;
			if (Evo_TotalTFT < this->m_BestSol_TotalTFT)
			{
				this->m_BestSol_Seq = Evo_Sol;
				this->m_BestSol_FacTFT = Evo_FacTFT;
				this->m_BestSol_TotalTFT = Evo_TotalTFT;
				this->m_BestSol_FacSpan = Evo_FacSpan;
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
				CurFacSpan = Evo_FacSpan;
				m_perEij = Evo_perEij;
			}
		}
	}

	//------------------the second stage-------------
	while (Base::GetElapsedProcessTime() - GetStartTime() < m_CPU_time)
	{
		int critfac = max_element(m_BestSol_FacSpan.begin(), m_BestSol_FacSpan.end()) - m_BestSol_FacSpan.begin();

		vector<int> Evo_seq(m_BestSol_Seq[critfac]);
		int Evo_span = m_BestSol_FacSpan[critfac];
		vector<vector<int>> Eij;

		vector<int> Des_seq;
		Destruction_2S(Des_seq, Evo_seq);
		int Evo_tft = GetEij(Evo_seq, Eij);
		Reconstruction_2S(Des_seq, Evo_seq, Eij, Evo_span, Evo_tft);
		LS1_IG2S(Evo_seq, Evo_tft, Evo_span, Eij);

		if (Evo_tft < m_BestSol_FacTFT[critfac])
		{
			m_BestSol_TotalTFT -= m_BestSol_FacTFT[critfac];
			m_BestSol_Seq[critfac] = Evo_seq;
			m_BestSol_FacSpan[critfac] = Evo_span;
			m_BestSol_FacTFT[critfac] = Evo_tft;
			m_BestSol_TotalTFT += m_BestSol_FacTFT[critfac];
		}
	}
	cout << (Base::GetElapsedProcessTime() - GetStartTime()) << ends;
	this->Check(this->m_BestSol_Seq, this->m_BestSol_FacTFT, this->m_BestSol_TotalTFT);
	return this->m_BestSol_TotalTFT;
}

void IG_2S::Destruction(vector<int>& des_seq, vector<vector<int>>& sol, vector<int>& facspan, vector<int>& factft, vector<vector<vector<int>>>& perEij)
{
	des_seq.resize(m_d1);
	int cnt_jobs = 0;
	int critfac = max_element(facspan.begin(), facspan.end()) - facspan.begin();
	vector<bool> changeflag_fac(this->m_Factories, false);
	changeflag_fac[critfac] = true;
	while (cnt_jobs < m_d1/2 && !sol[critfac].empty())
	{
		int extpos = rand() % sol[critfac].size();
		des_seq[cnt_jobs] = sol[critfac][extpos];
		sol[critfac].erase(sol[critfac].begin() + extpos);
		cnt_jobs++;
	}
	while (cnt_jobs < m_d1)
	{
		int extfac = -1;
		do
		{
			extfac = rand() % m_Factories;
		} while (sol[extfac].empty());
		int extpos = rand() % sol[extfac].size();
		des_seq[cnt_jobs] = sol[extfac][extpos];
		sol[extfac].erase(sol[extfac].begin() + extpos);
		changeflag_fac[extfac] = true;
		cnt_jobs++;
	}
	for (int f = 0; f < m_Factories; f++)
		if (1)//changeflag_fac[f]
		{
			factft[f] = GetEij(sol[f], perEij[f]);
			if (sol[f].empty())
				facspan[f] = 0;
			else				 
				facspan[f] = perEij[f][sol[f].size() - 1][this->m_Machines - 1];
		}
}

int IG_2S::Reconstruction(vector<int>& des_seq, 
	vector<vector<int>>& sol, vector<int>& facspan, vector<int>& factft, vector<vector<vector<int>>>& perEij)
{
	for (int j = 0; j < m_d1; j++)
	{
		int Job = des_seq[j];
		int BestPos = -1, BestFac = -1;
		int Min_TFT = FindBestInsertPosFacByMinIncreasedTFT(sol, factft, perEij, Job, BestFac, BestPos);
		factft[BestFac] = Min_TFT;
		InsertAJobToSeq(Job, BestPos, sol[BestFac], perEij[BestFac], factft[BestFac]);
		facspan[BestFac] = perEij[BestFac][sol[BestFac].size() - 1][this->m_Machines - 1];

		//reinsert the prvious or follow job in original factory
		int RePos;
		if (sol[BestFac].size() > 2)
		{
			if (BestPos == 0)
				RePos = BestPos + 1;
			else if (BestPos == sol[BestFac].size() - 1)
				RePos = BestPos - 1;
			else if (rand() % 2)
				RePos = BestPos + 1;
			else
				RePos = BestPos - 1;

			int ReJob = sol[BestFac][RePos];
			int BPos = -1;
			Min_TFT = FindBestInsertPosInsideFactory(sol[BestFac], RePos, perEij[BestFac], BPos);
			if (Min_TFT < factft[BestFac])
			{
				sol[BestFac].erase(sol[BestFac].begin() + RePos);
				if (BPos > RePos)
					sol[BestFac].insert(sol[BestFac].begin() + (BPos - 1), ReJob);
				else
					sol[BestFac].insert(sol[BestFac].begin() + BPos, ReJob);
				factft[BestFac] = Min_TFT;
				GetEij(sol[BestFac], perEij[BestFac]);
				facspan[BestFac] = perEij[BestFac][sol[BestFac].size() - 1][this->m_Machines - 1];
			}
		}
	}
	int total_tft = accumulate(factft.begin(), factft.end(), 0);
	return total_tft;
}

void IG_2S::Destruction_2S(vector<int>& des_seq, vector<int>& seq)
{
	int limitjobs = min(seq.size(), m_d2);
	des_seq.resize(limitjobs);
	for (int j = 0; j < limitjobs; j++)
	{
		int pos = rand() % seq.size();
		des_seq[j] = seq[pos];
		seq.erase(seq.begin() + pos);
	}
}

void IG_2S::Reconstruction_2S(vector<int>& des_seq, vector<int>& seq, vector<vector<int>>& Eij, int& mspan,int &tft)
{
	while (!des_seq.empty())
	{
		int pos = rand() % des_seq.size();
		int job = des_seq[pos];
		des_seq.erase(des_seq.begin() + pos);
		int bestpos = -1;
		int resulting_tft = GetMinTFTAfterInsertSingleFac(seq, job, bestpos, Eij);
		InsertAJobToSeq(job, bestpos, seq, Eij, tft);

		int RePos;
		if (seq.size() > 3)
		{
			if (bestpos == 0)
				RePos = bestpos + 1;
			else if (bestpos == seq.size() - 1)
				RePos = bestpos - 1;
			else if (rand() % 2)
				RePos = bestpos + 1;
			else
				RePos = bestpos - 1;

			int ReJob = seq[RePos];
			int BPos = -1;
			int remintft = FindBestInsertPosInsideFactory(seq, RePos, Eij, BPos);
			if (remintft < tft)
			{
				seq.erase(seq.begin() + RePos);
				if (BPos > RePos)
					seq.insert(seq.begin() + (BPos - 1), ReJob);
				else
					seq.insert(seq.begin() + BPos, ReJob);
				tft = remintft;
				GetEij(seq, Eij);
			}
		}
	}
	mspan = Eij[seq.size() - 1][m_Machines - 1];
}



