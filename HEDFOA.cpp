#include "HEDFOA.h"

void HEDFOA::SetParameter_HEDFOA(int SetNP, int SetSN, int Setd, float SetT, int SetC)
{
	this->m_NP = SetNP;
	this->m_SN = SetSN;
	this->m_d = Setd;
	this->m_T = SetT;
	this->m_CPUtime = SetC * m_Machines * m_Jobs;

	this->m_BestSol_Seq.resize(this->m_Factories);
	this->m_BestSol_FacTFT.resize(this->m_Factories);
	this->m_BestSol_FacSpan.resize(this->m_Factories);
	this->m_BestSol_TotalTFT = -1;

	this->m_SwarmSol.resize(this->m_NP);
	this->m_SwarmFacSpan.resize(this->m_NP);
	this->m_SwarmFacTft.resize(this->m_NP);
	this->m_SwarmTotaltft.resize(this->m_NP);
}

int HEDFOA::HEDFOA_Run()
{
	long StartTime = Base::GetElapsedProcessTime();
	this->GetTotalPTime();
	this->GetTotalLoad();

	this->m_BestSol_TotalTFT = EHPF2(m_BestSol_Seq, m_BestSol_FacTFT, m_BestSol_FacSpan);
	for (int n = 0; n < this->m_NP; n++)//每个果蝇中心位置初始化
	{
		this->m_SwarmSol[n] = this->m_BestSol_Seq;
		this->m_SwarmFacTft[n] = this->m_BestSol_FacTFT;
		this->m_SwarmTotaltft[n] = this->m_BestSol_TotalTFT;
		this->m_SwarmFacSpan[n] = this->m_BestSol_FacSpan;
	}
	//更新所有果蝇群的最优
	double Temperature = this->m_T * this->m_TotalLoad / (this->m_Machines * this->m_Jobs * 10);
	while (Base::GetElapsedProcessTime() - StartTime < this->m_CPUtime)
	{
		for (int n = 0; n < this->m_NP; n++)//m_NP 个果蝇群
		{
			vector<vector<vector<int>>> FFSol(this->m_SN);//每个果蝇群有m_SN个果蝇，即每个父代解生成若干个子代解
			vector<vector<int>> FFFacTft(this->m_SN);
			vector<vector<int>> FFFacSpan(this->m_SN);
			vector<int> FFTotalTFT(this->m_SN);

			vector < vector<vector<vector<int>>>> SF_Eij(this->m_SN);
			for (int sf = 0; sf < this->m_SN; sf++)
			{
				FFSol[sf] = this->m_SwarmSol[n];
				FFFacTft[sf] = this->m_SwarmFacTft[n];
				FFTotalTFT[sf] = this->m_SwarmTotaltft[n];
				FFFacSpan[sf] = this->m_SwarmFacSpan[n];

				FFTotalTFT[sf] = this->Smell_based_foraging(FFSol[sf], FFFacSpan[sf], FFFacTft[sf], SF_Eij[sf]);
			}
			int BestFF = min_element(FFTotalTFT.begin(), FFTotalTFT.end()) - FFTotalTFT.begin();//该果蝇群中探索的最好的一个解

			LocalSearch_HEDFOA(FFSol[BestFF], FFFacTft[BestFF], FFFacSpan[BestFF], FFTotalTFT[BestFF], SF_Eij[BestFF]);

			if (FFTotalTFT[BestFF] < this->m_SwarmTotaltft[n])
			{
				this->m_SwarmSol[n] = FFSol[BestFF];
				this->m_SwarmFacTft[n] = FFFacTft[BestFF];
				this->m_SwarmFacSpan[n] = FFFacSpan[BestFF];
				this->m_SwarmTotaltft[n] = FFTotalTFT[BestFF];
			}
			else
			{
				float RandNumber = (float)rand() / RAND_MAX;
				float Index = 1.0f * exp((this->m_SwarmTotaltft[n] - FFTotalTFT[BestFF]) / Temperature);
				if (RandNumber < Index)
				{
					this->m_SwarmSol[n] = FFSol[BestFF];
					this->m_SwarmFacTft[n] = FFFacTft[BestFF];
					this->m_SwarmFacSpan[n] = FFFacSpan[BestFF];
					this->m_SwarmTotaltft[n] = FFTotalTFT[BestFF];
				}
			}
		}

		int bestswarm = min_element(this->m_SwarmTotaltft.begin(), this->m_SwarmTotaltft.end()) - this->m_SwarmTotaltft.begin();
		if (this->m_SwarmTotaltft[bestswarm] < this->m_BestSol_TotalTFT)
		{
			this->m_BestSol_Seq = this->m_SwarmSol[bestswarm];
			this->m_BestSol_FacTFT = this->m_SwarmFacTft[bestswarm];
			this->m_BestSol_TotalTFT = this->m_SwarmTotaltft[bestswarm];
			this->m_BestSol_FacSpan = this->m_SwarmFacSpan[bestswarm];

			//每个果蝇中心位置初始化
			for (int n = 0; n < this->m_NP; n++)
			{
				if (n == bestswarm)
					continue;
				this->m_SwarmSol[n] = this->m_BestSol_Seq;
				this->m_SwarmFacTft[n] = this->m_BestSol_FacTFT;
				this->m_SwarmTotaltft[n] = this->m_BestSol_TotalTFT;
				this->m_SwarmFacSpan[n] = this->m_BestSol_FacSpan;
			}
		}
	}
	this->Check(this->m_BestSol_Seq, this->m_BestSol_FacTFT, this->m_BestSol_TotalTFT);
	return this->m_BestSol_TotalTFT;
}

int HEDFOA::Smell_based_foraging(vector<vector<int>>& Solution, vector<int>& FacSpan, vector<int> &FacTFT, 
	vector<vector<vector<int>>>& perEij)
{
	vector<int> des_seq(this->m_d);
	//先抽离工件
	int critfac = max_element(FacSpan.begin(), FacSpan.end()) - FacSpan.begin();

	int extpos = rand() % Solution[critfac].size();
	des_seq[0] = Solution[critfac][extpos];
	Solution[critfac].erase(Solution[critfac].begin() + extpos);

	for (int j = 1; j < m_d; j++)
	{
		int extfac = -1;
		do
		{
			extfac = rand() % this->m_Factories;
		} while (Solution[extfac].empty());
		int extpos = rand() % Solution[extfac].size();
		des_seq[j] = Solution[extfac][extpos];
		Solution[extfac].erase(Solution[extfac].begin() + extpos);
	}
	perEij.resize(this->m_Factories);
	for (int f = 0; f < m_Factories; f++)
	{
		FacTFT[f] = GetEij(Solution[f], perEij[f]);
	}

	//再把抽取的工件重新插入到部分序列中
	for (int j = 0; j < m_d; j++)
	{
		int Job = des_seq[j];
		int BestPos = -1, BestFac = -1;
		int Min_TFT = FindBestInsertPosFacByMinIncreasedTFT(Solution, FacTFT, perEij, Job, BestFac, BestPos);
		FacTFT[BestFac] = Min_TFT;
		InsertAJobToSeq(Job, BestPos, Solution[BestFac], perEij[BestFac], FacTFT[BestFac]);

		//reinsert the prvious and follow job in original factory
		if (Solution[BestFac].size() > 2)
		{
			vector<int> re_desseq;

			if (BestPos == 0)
			{
				re_desseq.push_back(Solution[BestFac][1]);
				Solution[BestFac].erase(Solution[BestFac].begin() + 1);
			}
			else if (BestPos == Solution[BestFac].size() - 1)
			{
				re_desseq.push_back(Solution[BestFac][Solution[BestFac].size() - 1]);
				Solution[BestFac].pop_back();
			}
			else 
			{
				re_desseq.push_back(Solution[BestFac][BestPos - 1]);
				re_desseq.push_back(Solution[BestFac][BestPos + 1]);
				Solution[BestFac].erase(Solution[BestFac].begin() + (BestPos + 1));
				Solution[BestFac].erase(Solution[BestFac].begin() + (BestPos - 1));
			}

			FacTFT[BestFac] = GetEij(Solution[BestFac], perEij[BestFac]);
			for (int rej = 0; rej < re_desseq.size(); rej++)
			{
				int re_job = re_desseq[rej];
				int re_bestfac = -1, re_bestpos = -1;
				int re_mintft = FindBestInsertPosFacByMinIncreasedTFT(Solution, FacTFT, perEij, re_job, re_bestfac, re_bestpos);
				FacTFT[re_bestfac] = re_mintft;
				InsertAJobToSeq(re_job, re_bestpos, Solution[re_bestfac], perEij[re_bestfac], FacTFT[re_bestfac]);
			}
		}
	}

	for (int f = 0; f < m_Factories; f++)
		FacSpan[f] = perEij[f][Solution[f].size() - 1][m_Machines - 1];
	int totaltft = accumulate(FacTFT.begin(), FacTFT.end(), 0);

	Check(Solution, FacTFT, totaltft);

	return totaltft;
}
