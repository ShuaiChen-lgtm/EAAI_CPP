#include "ConstructiveHeuristics.h"

int ConstructiveHeuristics::Con_Heu_Run(int type_conheu)
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
	case 6:this->m_BestSol_TotalTFT = EHPF2(this->m_BestSol_Seq, this->m_BestSol_FacTFT); break;
	default:
		break;
	}
	Check(this->m_BestSol_Seq, this->m_BestSol_FacTFT, this->m_BestSol_TotalTFT);
	return this->m_BestSol_TotalTFT;
}

int ConstructiveHeuristics::DLR(vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij)
{
	vector<int> Unsche_Seq(this->m_Jobs);
	GenerateSeedSeq_IF0(Unsche_Seq);

	Sol.resize(this->m_Factories);
	FacTFT.resize(this->m_Factories);
	vector<int> FacSpan(this->m_Factories);
	for (int f = 0; f < this->m_Factories; f++)
	{
		int Job = Unsche_Seq[f];
		Sol[f].push_back(Job);
		FacSpan[f] = this->m_TotalPTime[Job];
		FacTFT[f] = GetEij(Sol[f], perEij[f]);
	}
	Unsche_Seq.erase(Unsche_Seq.begin(), Unsche_Seq.begin() + this->m_Factories);

	while (Unsche_Seq.size() > 1)
	{
		int SelectFac = min_element(FacSpan.begin(), FacSpan.end()) - FacSpan.begin();
		int BestJobPos = FindBestJobPosFor_DLR(Sol[SelectFac], Unsche_Seq, FacSpan[SelectFac], perEij[SelectFac]);
		PushAJobToSeq(Unsche_Seq[BestJobPos], Sol[SelectFac], perEij[SelectFac], FacTFT[SelectFac]);
		Unsche_Seq.erase(Unsche_Seq.begin() + BestJobPos);
	}
	int SelectFac = min_element(FacSpan.begin(), FacSpan.end()) - FacSpan.begin();
	PushAJobToSeq(Unsche_Seq[0], Sol[SelectFac], perEij[SelectFac], FacTFT[SelectFac]);

	int totaltft = accumulate(FacTFT.begin(), FacTFT.end(), 0);
	return totaltft;
}

int ConstructiveHeuristics::DLR_DNEH(float Ratio, vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij)
{
	Sol.resize(this->m_Factories);
	FacTFT.resize(this->m_Factories);
	vector<int> SeedPermutation(this->m_Jobs);
	GenerateSeedSeq_IF0(SeedPermutation);
	vector<int> FacSpan(this->m_Factories);
	for (int f = 0; f < m_Factories; f++)
	{
		int Job = SeedPermutation[f];
		Sol[f].push_back(Job);
		FacSpan[f] = this->m_TotalPTime[Job];
		FacTFT[f] = GetEij(Sol[f], perEij[f]);
	}
	SeedPermutation.erase(SeedPermutation.begin(), SeedPermutation.begin() + this->m_Factories);
	while (SeedPermutation.size() > m_Jobs * Ratio)
	{
		int SelectFac = min_element(FacSpan.begin(), FacSpan.end()) - FacSpan.begin();
		int BestJobPos = FindBestJobPosFor_DLR(Sol[SelectFac], SeedPermutation, FacSpan[SelectFac], perEij[SelectFac]);
		PushAJobToSeq(SeedPermutation[BestJobPos], Sol[SelectFac], perEij[SelectFac], FacTFT[SelectFac]);
		SeedPermutation.erase(SeedPermutation.begin() + BestJobPos);
	}
	int num_remainJob = SeedPermutation.size();
	for (int j = 0; j < num_remainJob; j++)
	{
		int Job = SeedPermutation[j];
		int BestPos = -1, BestFac = -1;
		int Min_TFT = FindBestInsertPosFacByMinIncreasedTFT(Sol, FacTFT, perEij, Job, BestFac, BestPos);
		this->InsertAJobToSeq(Job, BestPos, Sol[BestFac], perEij[BestFac], FacTFT[BestFac]);
	}
	int totaltft = accumulate(FacTFT.begin(), FacTFT.end(), 0);
	return totaltft;
}

int ConstructiveHeuristics::NEH2_en(vector<vector<int>>& Sol, vector<int>& FacTFT, vector<int>& FacSpan, vector<vector<vector<int>>>& perEij)
{
	Sol.resize(this->m_Factories);
	FacTFT.resize(this->m_Factories, 0);
	FacSpan.resize(this->m_Factories, 0);
	vector<int> SeedPermutation(this->m_Jobs);
	DescendingOrder(m_Jobs, m_TotalPTime, SeedPermutation);
	for (int f = 0; f < m_Factories; f++)
		this->PushAJobToSeq(SeedPermutation[f], Sol[f], perEij[f], FacTFT[f]);
	for (int j = this->m_Factories; j < this->m_Jobs; j++)
	{
		int Job = SeedPermutation[j];
		int BestPos = -1, BestFac = -1;
		int Min_TFT = FindBestInsertPosFacByMinIncreasedTFT(Sol, FacTFT, perEij, Job, BestFac, BestPos);
		this->InsertAJobToSeq(Job, BestPos, Sol[BestFac], perEij[BestFac], FacTFT[BestFac]);

		int RePos;
		if (Sol[BestFac].size() > 2)
		{
			if (BestPos == 0)
				RePos = BestPos + 1;
			else if (BestPos == Sol[BestFac].size() - 1)
				RePos = BestPos - 1;
			else if (rand() % 2)
				RePos = BestPos + 1;
			else
				RePos = BestPos - 1;

			int ReJob = Sol[BestFac][RePos];
			int BPos = -1;
			Min_TFT = FindBestInsertPosInsideFactory(Sol[BestFac], RePos, perEij[BestFac], BPos);
			if (Min_TFT < FacTFT[BestFac])
			{
				Sol[BestFac].erase(Sol[BestFac].begin() + RePos);
				if (BPos > RePos)
					Sol[BestFac].insert(Sol[BestFac].begin() + (BPos - 1), ReJob);
				else
					Sol[BestFac].insert(Sol[BestFac].begin() + BPos, ReJob);
				FacTFT[BestFac] = Min_TFT;
				GetEij(Sol[BestFac], perEij[BestFac]);
			}
		}
	}
	int TotalTFT = accumulate(FacTFT.begin(), FacTFT.end(), 0);
	for (int f = 0; f < m_Factories; f++)
		FacSpan[f] = perEij[f][Sol[f].size() - 1][this->m_Machines - 1];
	return TotalTFT;
}

int ConstructiveHeuristics::DNEH(vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij)
{
	Sol.resize(this->m_Factories);
	FacTFT.resize(this->m_Factories, 0);
	vector<int> SeedPermutation(this->m_Jobs);
	GenerateSeedSeq_IF0(SeedPermutation);
	for (int f = 0; f < m_Factories; f++)
		this->PushAJobToSeq(SeedPermutation[f], Sol[f], perEij[f], FacTFT[f]);
	for (int j = this->m_Factories; j < this->m_Jobs; j++)
	{
		int Job = SeedPermutation[j];
		int BestPos = -1, BestFac = -1;
		int Min_TFT = FindBestInsertPosFacByMinIncreasedTFT(Sol, FacTFT, perEij, Job, BestFac, BestPos);
		this->InsertAJobToSeq(Job, BestPos, Sol[BestFac], perEij[BestFac], FacTFT[BestFac]);
	}
	int TotalTFT = accumulate(FacTFT.begin(), FacTFT.end(), 0);
	return TotalTFT;
}

int ConstructiveHeuristics::RandDNEH(vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij)
{
	Sol.resize(this->m_Factories);
	FacTFT.resize(this->m_Factories, 0);
	vector<int> SeedPermutation(this->m_Jobs);
	GenerateSeedSeq_IF0(SeedPermutation);
	for (int f = 0; f < m_Factories; f++)
	{
		int r = rand() % SeedPermutation.size();
		this->PushAJobToSeq(SeedPermutation[r], Sol[f], perEij[f], FacTFT[f]);
		SeedPermutation.erase(SeedPermutation.begin() + r);
	}
	for (int j = 0; j < this->m_Jobs - m_Factories; j++)
	{
		int Job = SeedPermutation[j];
		int BestPos = -1, BestFac = -1;
		int Min_TFT = FindBestInsertPosFacByMinIncreasedTFT(Sol, FacTFT, perEij, Job, BestFac, BestPos);
		this->InsertAJobToSeq(Job, BestPos, Sol[BestFac], perEij[BestFac], FacTFT[BestFac]);
	}
	int TotalTFT = accumulate(FacTFT.begin(), FacTFT.end(), 0);
	return TotalTFT;
}

int ConstructiveHeuristics::RandDLR(vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij)
{
	Sol.clear();
	Sol.resize(m_Factories);
	FacTFT.clear();
	FacTFT.resize(m_Factories, 0);

	vector<int> Unsche_Seq(this->m_FullSeq);
	random_shuffle(Unsche_Seq.begin(), Unsche_Seq.end());

	Sol.resize(this->m_Factories);
	FacTFT.resize(this->m_Factories);
	vector<int> FacSpan(this->m_Factories);
	for (int f = 0; f < this->m_Factories; f++)
	{
		int Job = Unsche_Seq[f];
		Sol[f].push_back(Job);
		FacSpan[f] = this->m_TotalPTime[Job];
		FacTFT[f] = GetEij(Sol[f], perEij[f]);
	}
	Unsche_Seq.erase(Unsche_Seq.begin(), Unsche_Seq.begin() + this->m_Factories);

	while (Unsche_Seq.size() > 1)
	{
		int SelectFac = min_element(FacSpan.begin(), FacSpan.end()) - FacSpan.begin();
		int BestJobPos = FindBestJobPosFor_DLR(Sol[SelectFac], Unsche_Seq, FacSpan[SelectFac], perEij[SelectFac]);
		PushAJobToSeq(Unsche_Seq[BestJobPos], Sol[SelectFac], perEij[SelectFac], FacTFT[SelectFac]);
		Unsche_Seq.erase(Unsche_Seq.begin() + BestJobPos);
	}
	int SelectFac = min_element(FacSpan.begin(), FacSpan.end()) - FacSpan.begin();
	PushAJobToSeq(Unsche_Seq[0], Sol[SelectFac], perEij[SelectFac], FacTFT[SelectFac]);

	int totaltft = accumulate(FacTFT.begin(), FacTFT.end(), 0);
	return totaltft;
}

int ConstructiveHeuristics::DPFE(vector<vector<int>>& Sol,vector<int>& FacTFT, vector<vector<vector<int>>>& perEij)
{
	vector<int> Unsche_Seq(this->m_Jobs);
	GenerateSeedSeqForDPFE(Unsche_Seq);

	Sol.resize(this->m_Factories);
	FacTFT.resize(this->m_Factories);
	vector<int> FacSpan(this->m_Factories);
	for (int f = 0; f < this->m_Factories; f++)
	{
		int Job = Unsche_Seq[f];
		Sol[f].push_back(Job);
		FacSpan[f] = this->m_TotalPTime[Job];
		FacTFT[f] = GetEij(Sol[f], perEij[f]);
	}
	Unsche_Seq.erase(Unsche_Seq.begin(), Unsche_Seq.begin() + this->m_Factories);
	
	int limitjobs = m_Factories;
	while (Unsche_Seq.size() > limitjobs)
	{
		int SelectFac = min_element(FacSpan.begin(), FacSpan.end()) - FacSpan.begin();

		int BestJobPos = FindBestJobPosForDPFE(Sol[SelectFac], Unsche_Seq, perEij[SelectFac], FacSpan[SelectFac]);
		PushAJobToSeq(Unsche_Seq[BestJobPos], Sol[SelectFac], perEij[SelectFac], FacTFT[SelectFac]);
		Unsche_Seq.erase(Unsche_Seq.begin() + BestJobPos);
	}
	for (int j = 0; j < limitjobs; j++)
	{
		int Job = Unsche_Seq[j];
		int BestPos = -1, BestFac = -1;
		int Min_TFT = FindBestInsertPosFacByMinIncreasedTFT(Sol, FacTFT, perEij, Job, BestFac, BestPos);
		this->InsertAJobToSeq(Job, BestPos, Sol[BestFac], perEij[BestFac], FacTFT[BestFac]);
	}
	int totaltft = accumulate(FacTFT.begin(), FacTFT.end(), 0);
	return totaltft;
}

int ConstructiveHeuristics::RandDPFE(vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij)
{
	vector<int> Unsche_Seq(m_Jobs);
	GenerateSeedSeqForDPFE(Unsche_Seq);

	Sol.resize(this->m_Factories);
	FacTFT.resize(this->m_Factories);
	vector<int> FacSpan(this->m_Factories);
	for (int f = 0; f < this->m_Factories; f++)
	{
		int pos = rand() % Unsche_Seq.size();
		int Job = Unsche_Seq[pos];
		Sol[f].push_back(Job);
		FacSpan[f] = this->m_TotalPTime[Job];
		FacTFT[f] = GetEij(Sol[f], perEij[f]);
		Unsche_Seq.erase(Unsche_Seq.begin() + pos);
	}

	while (Unsche_Seq.size() > this->m_Factories)
	{
		int SelectFac = min_element(FacSpan.begin(), FacSpan.end()) - FacSpan.begin();
		int BestJobPos = FindBestJobPosForDPFE(Sol[SelectFac], Unsche_Seq, perEij[SelectFac], FacSpan[SelectFac]);
		PushAJobToSeq(Unsche_Seq[BestJobPos], Sol[SelectFac], perEij[SelectFac], FacTFT[SelectFac]);
		Unsche_Seq.erase(Unsche_Seq.begin() + BestJobPos);
	}

	for (int j = 0; j < m_Factories; j++)
	{
		int Job = Unsche_Seq[j];
		int BestPos = -1, BestFac = -1;
		int Min_TFT = FindBestInsertPosFacByMinIncreasedTFT(Sol, FacTFT, perEij, Job, BestFac, BestPos);
		this->InsertAJobToSeq(Job, BestPos, Sol[BestFac], perEij[BestFac], FacTFT[BestFac]);
	}
	int totaltft = accumulate(FacTFT.begin(), FacTFT.end(), 0);
	return totaltft;
}

int ConstructiveHeuristics::DPF(vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij)
{
	vector<int> Unsche_Seq(this->m_Jobs);
	GenerateSeedSeqForDPFE(Unsche_Seq);

	Sol.resize(this->m_Factories);
	FacTFT.resize(this->m_Factories);
	vector<int> FacSpan(this->m_Factories);
	for (int f = 0; f < this->m_Factories; f++)
	{
		int Job = Unsche_Seq[f];
		Sol[f].push_back(Job);
		FacSpan[f] = this->m_TotalPTime[Job];
		FacTFT[f] = GetEij(Sol[f], perEij[f]);
	}
	Unsche_Seq.erase(Unsche_Seq.begin(), Unsche_Seq.begin() + this->m_Factories);

	while (Unsche_Seq.size())
	{
		int SelectFac = min_element(FacSpan.begin(), FacSpan.end()) - FacSpan.begin();
		int BestJobPos = FindBestJobPosForDPFE(Sol[SelectFac], Unsche_Seq, perEij[SelectFac], FacSpan[SelectFac]);
		PushAJobToSeq(Unsche_Seq[BestJobPos], Sol[SelectFac], perEij[SelectFac], FacTFT[SelectFac]);
		Unsche_Seq.erase(Unsche_Seq.begin() + BestJobPos);
	}
	int totaltft = accumulate(FacTFT.begin(), FacTFT.end(), 0);
	return totaltft;
}

int ConstructiveHeuristics::HPF23(float lamda, float mu, vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij)
{
	vector<int> InitSeq;
	this->Rule_HPF23(InitSeq, lamda, mu);//
	Sol.clear();
	Sol.resize(this->m_Factories);
	FacTFT.clear();
	FacTFT.resize(m_Factories);
	vector<int> FacTotalTime(this->m_Factories, 0);

	int aveload = m_TotalLoad / m_Factories;
	int curfac = 0;
	int curpos = 0;
	while (curfac < m_Factories && curpos < m_Jobs)
	{
		if (FacTotalTime[curfac] + this->m_TotalPTime[InitSeq[curpos]] > aveload && curfac != this->m_Factories - 1)
		{
			if (FacTotalTime[curfac] + this->m_TotalPTime[InitSeq[curpos]] - aveload > aveload - FacTotalTime[curfac])
			{
				curfac++;
				continue;
			}
		}
		Sol[curfac].push_back(InitSeq[curpos]);
		FacTotalTime[curfac] += m_TotalPTime[InitSeq[curpos]];
		curpos++;
	}

	for (int f = 0; f < m_Factories; f++)
	{
		vector<int> CopySol(Sol[f]);
		int nJobs = Sol[f].size();
		Sol[f].clear();
		int Job = CopySol[0];
		PushAJobToSeq(Job, Sol[f], perEij[f], FacTFT[f]);
		for (int j = 1; j < nJobs; j++)
		{
			Job = CopySol[j];
			int bestpos = -1;
			int resulting_tft = GetMinTFTAfterInsertSingleFac(Sol[f], Job, bestpos, perEij[f]);
			InsertAJobToSeq(Job, bestpos, Sol[f], perEij[f], FacTFT[f]);
		}
	}
	int totaltft = accumulate(FacTFT.begin(), FacTFT.end(), 0);
	return totaltft; 
}

void ConstructiveHeuristics::GenerateSeedSeq_EHPF2(vector<int>& SeedSeq,float lamda)
{
	Base::Pair<float>* ch = new Base::Pair<float>[this->m_Jobs];
	for (int j = 0; j < this->m_Jobs; j++)
	{
		int Job = j;

		float index = 0;
		for (int m = 0; m < m_Machines; m++)
			index += 1.0f * (m_Machines - m - 1) * this->m_Ptime[Job][m] / (m_Machines - 1);
		index = index * lamda * 2;
		for (int m = 0; m < m_Machines; m++)
			index = index + (1 - lamda) * this->m_Ptime[Job][m];

		ch[j].dim = j;
		ch[j].value = index;
	}

	sort(ch, ch + this->m_Jobs, Base::PairLess<float>());
	for (int j = 0; j < m_Jobs; j++)
		SeedSeq[j] = ch[j].dim;
	delete[] ch;
}

int ConstructiveHeuristics::EHPF2(vector<vector<int>>& Sol, vector<int>& FacTFT)
{
	vector<vector<vector<int>>> perEij(this->m_Factories);
	vector<int> SeedPermutation(this->m_Jobs);
	GenerateSeedSeq_EHPF2(SeedPermutation, 0.55f);
	vector<int> Facload(this->m_Factories, 0);
	FacTFT.resize(this->m_Factories);
	Sol.resize(this->m_Factories);
	for (int f = 0; f < m_Factories; f++)
	{
		this->PushAJobToSeq(SeedPermutation[f], Sol[f], perEij[f], FacTFT[f]);
		Facload[f] = this->m_TotalPTime[SeedPermutation[f]];
	}
	SeedPermutation.erase(SeedPermutation.begin(), SeedPermutation.begin() + this->m_Factories);

	GetTotalLoad();
	float average_load = this->m_TotalLoad / this->m_Factories;
	int cur_fac = 0;
	int curfac_nextpos = 1;
	while (!SeedPermutation.empty())
	{
		int bestchoosedpos = this->FindBestJobPosForEHPF2(perEij[cur_fac], SeedPermutation, curfac_nextpos, 0.7f);
		if (Facload[cur_fac] + this->m_TotalPTime[SeedPermutation[bestchoosedpos]] > average_load && cur_fac != this->m_Factories - 1)
		{
			if (Facload[cur_fac] + this->m_TotalPTime[SeedPermutation[bestchoosedpos]] - average_load > average_load - Facload[cur_fac])
			{
				cur_fac++;
				curfac_nextpos = 1;
				continue;
			}
		}
		this->PushAJobToSeq(SeedPermutation[bestchoosedpos], Sol[cur_fac], perEij[cur_fac], FacTFT[cur_fac]);
		Facload[cur_fac] += this->m_TotalPTime[SeedPermutation[bestchoosedpos]];
		SeedPermutation.erase(SeedPermutation.begin() + bestchoosedpos);
		curfac_nextpos++;
	}

	for (int f = 0; f < m_Factories; f++)
	{
		int SizeofSolf = Sol[f].size();
		vector<int> CopySeq(Sol[f]);
		Sol[f].clear();
		Sol[f].push_back(CopySeq[0]);
		FacTFT[f] = GetEij(Sol[f], perEij[f]);
		for (int j = 1; j < SizeofSolf; j++)
		{
			int bestpos = -1;
			FacTFT[f] = GetMinTFTAfterInsertSingleFac(Sol[f], CopySeq[j], bestpos, perEij[f]);
			InsertAJobToSeq(CopySeq[j], bestpos, Sol[f], perEij[f], FacTFT[f]);

			if (Sol[f].size() > 2)
			{
				vector<int> re_desseq;

				if (bestpos == 0)
				{
					re_desseq.push_back(Sol[f][1]);
					Sol[f].erase(Sol[f].begin() + 1);
				}
				else if (bestpos == Sol[f].size() - 1)
				{
					re_desseq.push_back(Sol[f][Sol[f].size() - 1]);
					Sol[f].pop_back();
				}
				else
				{
					re_desseq.push_back(Sol[f][bestpos - 1]);
					re_desseq.push_back(Sol[f][bestpos + 1]);
					Sol[f].erase(Sol[f].begin() + (bestpos + 1));
					Sol[f].erase(Sol[f].begin() + (bestpos - 1));
				}

				FacTFT[f] = GetEij(Sol[f], perEij[f]);
				for (int rej = 0; rej < re_desseq.size(); rej++)
				{
					int re_job = re_desseq[rej];
					int re_bestpos = -1;
					GetMinTFTAfterInsertSingleFac(Sol[f], re_job, re_bestpos, perEij[f]);
					InsertAJobToSeq(re_job, re_bestpos, Sol[f], perEij[f], FacTFT[f]);
				}
			}
		}
	}

	int totaltft = accumulate(FacTFT.begin(), FacTFT.end(), 0);
	return totaltft;
}

int ConstructiveHeuristics::EHPF2(vector<vector<int>>& Sol, vector<int>& FacTFT, vector<int>& FacSpan)
{
	vector<vector<vector<int>>> perEij(this->m_Factories);
	vector<int> SeedPermutation(this->m_Jobs);
	GenerateSeedSeq_EHPF2(SeedPermutation, 0.55f);
	vector<int> Facload(this->m_Factories, 0);

	for (int f = 0; f < m_Factories; f++)
	{
		this->PushAJobToSeq(SeedPermutation[f], Sol[f], perEij[f], FacTFT[f]);
		Facload[f] = this->m_TotalPTime[SeedPermutation[f]];
	}
	SeedPermutation.erase(SeedPermutation.begin(), SeedPermutation.begin() + this->m_Factories);

	float average_load = this->m_TotalLoad / this->m_Factories;
	int cur_fac = 0;
	int curfac_nextpos = 1;
	while (!SeedPermutation.empty())
	{
		int bestchoosedpos = this->FindBestJobPosForEHPF2(perEij[cur_fac], SeedPermutation, curfac_nextpos, 0.7f);
		if (Facload[cur_fac] + this->m_TotalPTime[SeedPermutation[bestchoosedpos]] > average_load && cur_fac != this->m_Factories - 1)
		{
			if (Facload[cur_fac] + this->m_TotalPTime[SeedPermutation[bestchoosedpos]] - average_load > average_load - Facload[cur_fac])
			{
				cur_fac++;
				curfac_nextpos = 1;
				continue;
			}
		}
		this->PushAJobToSeq(SeedPermutation[bestchoosedpos], Sol[cur_fac], perEij[cur_fac], FacTFT[cur_fac]);
		Facload[cur_fac] += this->m_TotalPTime[SeedPermutation[bestchoosedpos]];
		SeedPermutation.erase(SeedPermutation.begin() + bestchoosedpos);
		curfac_nextpos++;
	}
	for (int f = 0; f < m_Factories; f++)
	{
		int SizeofSolf = Sol[f].size();
		vector<int> CopySeq(Sol[f]);
		Sol[f].clear();
		Sol[f].push_back(CopySeq[0]);
		FacTFT[f] = GetEij(Sol[f], perEij[f]);
		for (int j = 1; j < SizeofSolf; j++)
		{
			int bestpos = -1;
			FacTFT[f] = GetMinTFTAfterInsertSingleFac(Sol[f], CopySeq[j], bestpos, perEij[f]);
			InsertAJobToSeq(CopySeq[j], bestpos, Sol[f], perEij[f], FacTFT[f]);

			if (Sol[f].size() > 2)
			{
				vector<int> re_desseq;

				if (bestpos == 0)
				{
					re_desseq.push_back(Sol[f][1]);
					Sol[f].erase(Sol[f].begin() + 1);
				}
				else if (bestpos == Sol[f].size() - 1)
				{
					re_desseq.push_back(Sol[f][Sol[f].size() - 1]);
					Sol[f].pop_back();
				}
				else
				{
					re_desseq.push_back(Sol[f][bestpos - 1]);
					re_desseq.push_back(Sol[f][bestpos + 1]);
					Sol[f].erase(Sol[f].begin() + (bestpos + 1));
					Sol[f].erase(Sol[f].begin() + (bestpos - 1));
				}

				FacTFT[f] = GetEij(Sol[f], perEij[f]);
				for (int rej = 0; rej < re_desseq.size(); rej++)
				{
					int re_job = re_desseq[rej];
					int re_bestpos = -1;
					GetMinTFTAfterInsertSingleFac(Sol[f], re_job, re_bestpos, perEij[f]);
					InsertAJobToSeq(re_job, re_bestpos, Sol[f], perEij[f], FacTFT[f]);
				}
			}
		}
	}

	for (int f = 0; f < m_Factories; f++)
	{
		FacSpan[f] = perEij[f][Sol[f].size() - 1][m_Machines - 1];
	}
	int totaltft = accumulate(FacTFT.begin(), FacTFT.end(), 0);
	return totaltft;
}

int ConstructiveHeuristics::NEHR2A4(vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij)
{
	Sol.resize(this->m_Factories);
	FacTFT.resize(this->m_Factories, 0);

	vector<int> SeedPermutation(this->m_Jobs);

	DescendingOrder(m_Jobs, m_TotalPTime, SeedPermutation);

	for (int f = 0; f < m_Factories; f++)
		this->PushAJobToSeq(SeedPermutation[f], Sol[f], perEij[f], FacTFT[f]);
	for (int j = this->m_Factories; j < this->m_Jobs; j++)
	{
		int Job = SeedPermutation[j];
		int BestPos = -1, BestFac = -1;
		int Min_TFT = FindBestInsertPosFacByMinIncreasedTFT(Sol, FacTFT, perEij, Job, BestFac, BestPos);
		this->InsertAJobToSeq(Job, BestPos, Sol[BestFac], perEij[BestFac], FacTFT[BestFac]);
	}
	int TotalTFT = accumulate(FacTFT.begin(), FacTFT.end(), 0);
	return TotalTFT;
}

int ConstructiveHeuristics::RandNEHR2A4(vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij)
{
	Sol.resize(this->m_Factories);
	FacTFT.resize(this->m_Factories, 0);
	vector<int> SeedPermutation(m_FullSeq);
	random_shuffle(SeedPermutation.begin(), SeedPermutation.end());
	for (int f = 0; f < m_Factories; f++)
		this->PushAJobToSeq(SeedPermutation[f], Sol[f], perEij[f], FacTFT[f]);
	for (int j = this->m_Factories; j < this->m_Jobs; j++)
	{
		int Job = SeedPermutation[j];
		int BestPos = -1, BestFac = -1;
		int Min_TFT = FindBestInsertPosFacByMinIncreasedTFT(Sol, FacTFT, perEij, Job, BestFac, BestPos);
		this->InsertAJobToSeq(Job, BestPos, Sol[BestFac], perEij[BestFac], FacTFT[BestFac]);
	}
	int TotalTFT = accumulate(FacTFT.begin(), FacTFT.end(), 0);
	return TotalTFT;
}

int ConstructiveHeuristics::PWNEH(vector<vector<int>>& Sol, vector<int>& FacTFT, vector<int>& FacSpan, vector<vector<vector<int>>>& perEij)
{
	Sol.resize(this->m_Factories);
	FacTFT.resize(this->m_Factories);
	FacSpan.resize(this->m_Factories);
	vector<int> SeedPermutation;
	PW(SeedPermutation);

	for (int f = 0; f < m_Factories; f++)
		this->PushAJobToSeq(SeedPermutation[f], Sol[f], perEij[f], FacTFT[f]);
	for (int j = this->m_Factories; j < this->m_Jobs; j++)
	{
		int Job = SeedPermutation[j];
		int BestPos = -1, BestFac = -1;
		int Min_TFT = FindBestInsertPosFacByMinIncreasedTFT(Sol, FacTFT, perEij, Job, BestFac, BestPos);
		this->InsertAJobToSeq(Job, BestPos, Sol[BestFac], perEij[BestFac], FacTFT[BestFac]);
	}
	for (int f = 0; f < m_Factories; f++)
		FacSpan[f] = perEij[f][Sol[f].size() - 1][m_Machines - 1];
	return accumulate(FacTFT.begin(), FacTFT.end(), 0);
}

void ConstructiveHeuristics::GenerateSeedSeq_IF0(vector<int>& SeedSeq)
{
	Base::Pair<float>* ch = new Base::Pair<float>[this->m_Jobs];
	for (int j = 0; j < this->m_Jobs; j++)
	{
		int Job = j;

		vector<int> Ci(this->m_Machines);
		Ci[0] = this->m_Ptime[Job][0];
		for (int m = 1; m < m_Machines; m++)
			Ci[m] = Ci[m - 1] + this->m_Ptime[Job][m];

		float Index_Idletime = 0;
		for (int m = 1; m < m_Machines; m++)
			Index_Idletime += float(this->m_Machines * Ci[m - 1]) / (m + 1);

		float Index_IF0 = ((float)this->m_Jobs / this->m_Factories - 2) * Index_Idletime + Ci[this->m_Machines - 1];

		ch[j].dim = j;
		ch[j].value = Index_IF0;
	}

	sort(ch, ch + this->m_Jobs, Base::PairLess<float>());
	for (int j = 0; j < m_Jobs; j++)
		SeedSeq[j] = ch[j].dim;
	delete[] ch;
}

void ConstructiveHeuristics::GenerateSeedSeqSingle_IF0(vector<int>& SeedSeq)
{
	unsigned int Nsize = SeedSeq.size();
	Base::Pair<float>* ch = new Base::Pair<float>[Nsize];
	for (int j = 0; j < Nsize; j++)
	{
		int Job = SeedSeq[j];

		vector<int> Ci(this->m_Machines);
		Ci[0] = this->m_Ptime[Job][0];
		for (int m = 1; m < m_Machines; m++)
			Ci[m] = Ci[m - 1] + this->m_Ptime[Job][m];

		float Index_Idletime = 0.0f;
		for (int m = 1; m < m_Machines; m++)
			Index_Idletime += float(this->m_Machines * Ci[m - 1]) / (m + 1);

		float Index_IF0 = (Nsize - 2) * Index_Idletime + Ci[this->m_Machines - 1];

		ch[j].dim = Job;
		ch[j].value = Index_IF0;
	}

	sort(ch, ch + Nsize, Base::PairLess<float>());
	for (unsigned int j = 0; j < Nsize; j++)
		SeedSeq[j] = ch[j].dim;
	delete[] ch;
}

void ConstructiveHeuristics::GenerateSeedSeqForDPFE(vector<int>& SeedSeq)
{
	Base::Pair<float>* ch = new Base::Pair<float>[this->m_Jobs];
	for (int j = 0; j < this->m_Jobs; j++)
	{
		int Job = j;

		vector<int> Ci(this->m_Machines);
		Ci[0] = this->m_Ptime[Job][0];
		for (int m = 1; m < m_Machines; m++)
			Ci[m] = Ci[m - 1] + this->m_Ptime[Job][m];

		float Index_Idletime = 0;
		for (int m = 0; m < m_Machines; m++)
			Index_Idletime += 1.0f * m_Machines * (Ci[m] - m_Ptime[Job][m]) / (m + 1);

		float Index_IF0 = (this->m_Jobs/m_Factories) * Index_Idletime / (this->m_Machines*2.0f) + Ci[this->m_Machines - 1];

		ch[j].dim = j;
		ch[j].value = Index_IF0;
	}

	sort(ch, ch + this->m_Jobs, Base::PairLess<float>());
	for (int j = 0; j < m_Jobs; j++)
		SeedSeq[j] = ch[j].dim;
	delete[] ch;
}

int ConstructiveHeuristics::FindBestJobPosForDPFE(const vector<int>& CurSeq,
	const vector<int>& UnscheduledSeq, const vector<vector<int>>& Eij, int& ResultMakeSpan)
{
	int CurSize = CurSeq.size();
	int LPos = CurSize - 1;
	int num_Unsche = UnscheduledSeq.size();

	vector<float> CurWeight(this->m_Machines);

	for (int m = 0; m < this->m_Machines; m++)
	{
		float w1 = (float)this->m_Jobs / this->m_Factories;
		float w2 = (float)CurSize * (this->m_Machines - m - 1) / w1 + m + 1;
		CurWeight[m] = (float)this->m_Machines / w2; 
	}

	int BestJobPos = -1;
	float Min_IBtime = FLT_MAX;
	for (int j = 0; j < num_Unsche; j++)
	{
		int TestJob = UnscheduledSeq[j];

		float IBtime = 0;
		vector<int> Dij(this->m_Machines);

		if (CurSize > 0)
		{
			Dij[0] = max(Eij[LPos][0] + this->m_Ptime[TestJob][0], Eij[LPos][1]);
			for (int m = 1; m < this->m_Machines - 1; m++)
				Dij[m] = max(Dij[m - 1] + this->m_Ptime[TestJob][m], Eij[LPos][m + 1]);
			Dij[this->m_Machines - 1] = Dij[this->m_Machines - 2] + this->m_Ptime[TestJob][this->m_Machines - 1];

			for (int m = 0; m < this->m_Machines; m++)
			{
				IBtime += max(Dij[m] - Eij[LPos][m] - this->m_Ptime[TestJob][m], 0) * CurWeight[m];
			}
			float A = (float)(m_Jobs / m_Factories - CurSize) / (m_Machines * 2.0f);
			IBtime += IBtime * A + Dij[m_Machines - 1];
		}

		if (IBtime < Min_IBtime)
		{
			BestJobPos = j;
			Min_IBtime = IBtime;
			ResultMakeSpan = Dij[this->m_Machines - 1];
		}
	}
	return BestJobPos;
}

int ConstructiveHeuristics::PWIndexFun(vector<vector<int>> Eij, vector<int> Remain_Seq, int k, vector<int>& RemainPtime)
{
	int BestJobpos = 0;
	float minIndex = INT_MAX;
	for (int j = 0; j < Remain_Seq.size(); j++)
	{
		int Job = Remain_Seq[j];
		vector<int> Dij(this->m_Machines);
		if (k)
		{
			Dij[0] = max(Eij[k - 1][0] + this->m_Ptime[Job][0], Eij[k - 1][1]);
			for (int m = 1; m < this->m_Machines - 1; m++)
				Dij[m] = max(Dij[m - 1] + this->m_Ptime[Job][m], Eij[k - 1][m + 1]);
			Dij[this->m_Machines - 1] = Dij[this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];
		}
		else
		{
			Dij[0] = this->m_Ptime[Job][0];
			for (int m = 1; m < this->m_Machines; m++)
			{
				Dij[m] = Dij[m - 1] + this->m_Ptime[Job][m];
			}
		}

		vector<float> Weight(this->m_Machines);
		for (int m = 0; m < this->m_Machines; m++)
		{
			float w = m + 1 + float(k * (this->m_Machines - m - 1)) / (this->m_Jobs - 2);
			Weight[m] = this->m_Machines / w;
		}
		vector<float> NewJobptime(this->m_Machines);
		for (int m = 0; m < this->m_Machines; m++)
		{
			int Totaltime = RemainPtime[m] - this->m_Ptime[Job][m];
			NewJobptime[m] = Totaltime / float(this->m_Jobs - k - 1);
		}

		vector<float> D2ij(this->m_Machines);
		D2ij[0] = max(Dij[0] + NewJobptime[0], Dij[1]);
		for (int m = 1; m < this->m_Machines - 1; m++)
			D2ij[m] = max(D2ij[m - 1] + NewJobptime[m], Dij[m + 1]);
		D2ij[this->m_Machines - 1] = D2ij[this->m_Machines - 2] + NewJobptime[this->m_Machines - 1];

		float Index = 0;
		if (k)
		{
			for (int m = 0; m < this->m_Machines; m++)
			{
				Index += float((Dij[m] - Eij[k - 1][m] - this->m_Ptime[Job][m]) * Weight[m]) * (this->m_Jobs - k - 2) + Weight[m] * (D2ij[m] - Dij[m] - NewJobptime[m]);
			}
		}
		else
		{
			for (int m = 0; m < this->m_Machines; m++)
			{
				Index += float((Dij[m] - this->m_Ptime[Job][m]) * Weight[m]) * (this->m_Jobs - k - 2) + Weight[m] * (D2ij[m] - Dij[m] - NewJobptime[m]);
			}
		}
		if (Index < minIndex)
		{
			BestJobpos = j;
			minIndex = Index;
		}
	}

	for (int m = 0; m < this->m_Machines; m++)
		RemainPtime[m] -= this->m_Ptime[Remain_Seq[BestJobpos]][m];
	return BestJobpos;
}

void ConstructiveHeuristics::PW(vector<int>& Seq)
{
	Seq.clear();
	vector<int> remainseq(this->m_FullSeq);
	vector<int> remain_ptime(this->m_Machines, 0);
	for (int i = 0; i < m_Machines; i++)
	{
		for (int j = 0; j < m_Jobs; j++)
			remain_ptime[i] += this->m_Ptime[j][i];
	}
	vector<vector<int>> Eij;
	int num_remainjob = this->m_Jobs;
	int curpos = 0;
	int total_tft = 0;
	while (num_remainjob > 1)
	{
		int bestjobpos = PWIndexFun(Eij, remainseq, curpos, remain_ptime);
		PushAJobToSeq(remainseq[bestjobpos], Seq, Eij, total_tft);
		remainseq.erase(remainseq.begin() + bestjobpos);
		curpos++;
		num_remainjob--;
	}
	Seq.push_back(remainseq[0]);
}

int ConstructiveHeuristics::FindBestJobPosFor_DPWE(const vector<int>& CurSeq, const vector<int>& UnscheduledSeq, const vector<vector<int>>& Eij, int& ResultMakeSpan)
{
	int CurSize = CurSeq.size();
	int LPos = CurSize - 1;
	int num_Unsche = UnscheduledSeq.size();

	vector<float> CurWeight(this->m_Machines);
	for (int m = 0; m < this->m_Machines; m++)
	{
		float w1 = (float)(this->m_Jobs - 2) / m_Factories;
		float w2 = (float)CurSize * (this->m_Machines - m - 1) / w1 + m + 1;
		CurWeight[m] = (float)this->m_Machines / w2;
	}

	int BestJobPos = -1;
	float Min_IBtime = FLT_MAX;
	for (int j = 0; j < num_Unsche; j++)
	{
		int TestJob = UnscheduledSeq[j];

		float IBtime = 0;
		vector<int> Dij(this->m_Machines);

		if (CurSize > 0)
		{
			Dij[0] = max(Eij[LPos][0] + this->m_Ptime[TestJob][0], Eij[LPos][1]);
			for (int m = 1; m < this->m_Machines - 1; m++)
				Dij[m] = max(Dij[m - 1] + this->m_Ptime[TestJob][m], Eij[LPos][m + 1]);
			Dij[this->m_Machines - 1] = Dij[this->m_Machines - 2] + this->m_Ptime[TestJob][this->m_Machines - 1];

			for (int m = 0; m < this->m_Machines; m++)
			{
				IBtime += max(Dij[m] - Eij[LPos][m] - this->m_Ptime[TestJob][m], 0) * CurWeight[m];
			}
			float A = (float)num_Unsche / this->m_Factories;
			float B = max(A - 2, 0);
			IBtime += IBtime * B + Dij[this->m_Machines - 1];
		}

		if (IBtime < Min_IBtime)
		{
			BestJobPos = j;
			Min_IBtime = IBtime;
			ResultMakeSpan = Dij[this->m_Machines - 1];
		}
	}
	return BestJobPos;
}

void ConstructiveHeuristics::Rule_HPF23(vector<int>& Seq, float lamda, float mu)
{
	vector<int> InitSeq;
	this->Rule_HPF23_FirstJob(InitSeq, lamda);

	Seq.clear();
	Seq.push_back(InitSeq[0]);
	InitSeq.erase(InitSeq.begin());
	vector<vector<int>> Eij;
	this->GetEij(Seq, Eij);
	int Cnt = 1;
	int tft = -1;
	while (InitSeq.size() > 1)
	{
		int BestJobPos = FindBestJobPosForHPF23(Eij, InitSeq, Cnt, mu);
		PushAJobToSeq(InitSeq[BestJobPos], Seq, Eij, tft);
		InitSeq.erase(InitSeq.begin() + BestJobPos);
		Cnt++;
	}
	Seq.push_back(InitSeq[0]);
}

void ConstructiveHeuristics::Rule_HPF23_FirstJob(vector<int>& Permutation, float lamda)
{
	Permutation.clear();
	Permutation.resize(m_Jobs);
	Base::Pair<float>* ch = new Base::Pair<float>[m_Jobs];
	for (int j = 0; j < m_Jobs; j++)
	{
		ch[j].dim = j;
		ch[j].value = 0.0f;
		for (int m = 0; m < m_Machines; m++)
		{
			ch[j].value += float(m_Machines - m - 1) * m_Ptime[j][m] / (m_Machines - 1);
		}
		ch[j].value = ch[j].value * 2 * lamda + (1 - lamda) * this->m_TotalPTime[j];
	}
	sort(ch, ch + m_Jobs, Base::PairLess<float>());
	for (int j = 0; j < m_Jobs; j++)
		Permutation[j] = ch[j].dim;
	delete[]ch;
}

int ConstructiveHeuristics::FindBestJobPosForHPF23(const vector<vector<int>>& Eij, const vector<int>& Remain_Seq, int k,float mu)
{
	int BestJobpos = 0;
	float minIndex = INT_MAX;
	for (int j = 0; j < Remain_Seq.size(); j++)
	{
		int Job = Remain_Seq[j];
		vector<int> Dij(m_Machines);
		if (k)
		{
			Dij[0] = max(Eij[k - 1][0] + m_Ptime[Job][0], Eij[k - 1][1]);
			for (int m = 1; m < m_Machines - 1; m++)
				Dij[m] = max(Dij[m - 1] + m_Ptime[Job][m], Eij[k - 1][m + 1]);
			Dij[this->m_Machines - 1] = Dij[this->m_Machines - 2] + m_Ptime[Job][this->m_Machines - 1];
		}
		else
		{
			Dij[0] = m_Ptime[Job][0];
			for (int m = 1; m < this->m_Machines; m++)
			{
				Dij[m] = Dij[m - 1] + m_Ptime[Job][m];
			}
		}

		float Index = 0;
		for (int m = 0; m < this->m_Machines; m++)
		{
			Index += Dij[m] - Eij[k - 1][m] - m_Ptime[Job][m];
		}
		Index = Index * mu + (Dij[m_Machines - 1] - Eij[k - 1][m_Machines - 1]) * (1.0 - mu);
		if (Index < minIndex)
		{
			BestJobpos = j;
			minIndex = Index;
		}
		Dij.clear();
	}
	return BestJobpos;
}

int ConstructiveHeuristics::FindBestJobPosForEHPF2(const vector<vector<int>>& Eij,
	const vector<int> &RemainSeq, int NextPos, float mu)
{
	int sizeofjobs = RemainSeq.size();
	float bestindex = FLT_MAX;
	int bestchoosedpos = -1;
	for (int p = 0; p < sizeofjobs; p++)
	{
		float Sita = 0.0f;
		int PJob = RemainSeq[p];
		vector<int> Dpj(this->m_Machines);
		Dpj[0] = max(Eij[NextPos - 1][0] + this->m_Ptime[PJob][0], Eij[NextPos - 1][1]);
		for (int m = 1; m < this->m_Machines - 1; m++)
			Dpj[m] = max(Dpj[m - 1] + this->m_Ptime[PJob][m], Eij[NextPos - 1][m + 1]);
		Dpj[this->m_Machines - 1] = Dpj[this->m_Machines - 2] + this->m_Ptime[PJob][this->m_Machines - 1];
		for (int m = 0; m < this->m_Machines; m++)
			Sita += Dpj[m] - Eij[NextPos - 1][m] - this->m_Ptime[PJob][m];
		Sita = Sita * mu + (1 - mu) * (float)(Dpj[this->m_Machines - 1] - Eij[NextPos - 1][this->m_Machines - 1]);
		if (Sita < bestindex)
		{
			bestindex = Sita;
			bestchoosedpos = p;
		}
	}
	return bestchoosedpos;
}

int ConstructiveHeuristics::FindBestJobPosFor_DLR(const vector<int>& CurSeq, 
	const vector<int>& UnscheduledSeq, int& ResultMakeSpan, const vector<vector<int>>& Eij)
{
	int CurSize = CurSeq.size();
	int LPos = CurSize - 1;

	vector<float> CurWeight(this->m_Machines);
	for (int m = 1; m < this->m_Machines; m++)
	{
		float w1 = (float)this->m_Jobs / this->m_Factories - 2;
		float w2 = (float)CurSize * (this->m_Machines - m - 1) / w1 + m + 1;
		CurWeight[m] = (float)this->m_Machines / w2;
	}

	int BestJobPos = -1;
	float Min_INdexIdletime = FLT_MAX;
	float Min_IndexIF = FLT_MAX;
	int num_Unsche = UnscheduledSeq.size();
	for (int j = 0; j < num_Unsche; j++)
	{
		int TestJob = UnscheduledSeq[j];

		float Index_IT = 0;
		float FinIndex = 0;
		vector<int> Dij(this->m_Machines);

		if (CurSize > 0)
		{
			Dij[0] = max(Eij[LPos][0] + this->m_Ptime[TestJob][0], Eij[LPos][1]);
			for (int m = 1; m < this->m_Machines - 1; m++)
				Dij[m] = max(Dij[m - 1] + this->m_Ptime[TestJob][m], Eij[LPos][m + 1]);
			Dij[this->m_Machines - 1] = Dij[this->m_Machines - 2] + this->m_Ptime[TestJob][this->m_Machines - 1];

			for (int m = 1; m < this->m_Machines; m++)
			{
				float IdleTime = (float)max(Dij[m - 1] - Eij[LPos][m], 0);
				Index_IT += IdleTime * CurWeight[m];
			}
		}
		else
		{
			Dij[0] = this->m_Ptime[TestJob][0];
			for (int m = 1; m < this->m_Machines; m++)
			{
				Dij[m] = Dij[m - 1] + this->m_Ptime[TestJob][m];
			}
			for (int m = 1; m < this->m_Machines; m++)
			{
				Index_IT += (float)Dij[m - 1] * CurWeight[m];
			}
		}
		float A = 1.0f * this->m_Jobs / this->m_Factories;
		float B = A - CurSize - 2;///------------------------------此处为原论文公式，缺陷是可能小于0，经实验，B = max(A - CurSize - 2,0) 或 B = (float)(this->m_Jobs - CurSize - 2) / this->m_Factories 较好
		FinIndex = Index_IT * B + Dij[this->m_Machines - 1];

		if (FinIndex < Min_IndexIF)
		{
			BestJobPos = j;
			Min_IndexIF = FinIndex;
			Min_INdexIdletime = Index_IT;
			ResultMakeSpan = Dij[this->m_Machines - 1];
		}
		else if (FinIndex == Min_IndexIF && Index_IT < Min_INdexIdletime)
		{
			BestJobPos = j;
			Min_IndexIF = FinIndex;
			Min_INdexIdletime = Index_IT;
			ResultMakeSpan = Dij[this->m_Machines - 1];
		}
	}
	return BestJobPos;
}

