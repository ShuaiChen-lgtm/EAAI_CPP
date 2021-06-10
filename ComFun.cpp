#include "ComFun.h"

ComFun::ComFun()
{
}

ComFun::~ComFun()
{
}

void ComFun::GetOrderSeq(int num_jobs, vector<int>& OrderSeq)
{
	OrderSeq.clear();
	OrderSeq.resize(num_jobs);
	for (int j = 0; j < num_jobs; j++)
		OrderSeq[j] = j;
}

void ComFun::GetRandomSeq(int num_jobs, vector<int>& RefSeq)
{
	GetOrderSeq(num_jobs, RefSeq);
	random_shuffle(RefSeq.begin(), RefSeq.end());
}

bool ComFun::JudgeExist(int Job, const vector<int>& Terget_Seq)
{
	int sizeofseq = Terget_Seq.size();
	for (int j = 0; j < sizeofseq; j++)
		if (Job == Terget_Seq[j])
			return true;
	return false;
}

bool ComFun::JudgeSameForTwoSol(const vector<vector<int>>& _1st_Sol, const vector<vector<int>>& _2nd_Sol)
{
	for (int f = 0; f < m_Factories; f++)
	{
		int Size_of_Sol1st_f = _1st_Sol[f].size();
		int Size_of_Sol2nd_f = _2nd_Sol[f].size();
		if (Size_of_Sol1st_f != Size_of_Sol2nd_f)
			return true;
		else
		{
			for (int j = 0; j < Size_of_Sol1st_f; j++)
			{
				if (_1st_Sol[f][j] != _2nd_Sol[f][j])
					return true;
			}
		}
	}
	return false;
}

bool ComFun::JudgeSameForTwoSol(const vector<vector<int>>& _1st_Sol, const vector<int>& _1st_factft, int _1st_TotalTFT,
	const vector<vector<int>>& _2nd_Sol, const vector<int>& _2nd_factft, int _2nd_TotalTFT)
{
	if (_1st_TotalTFT != _2nd_TotalTFT)
		return true;
	else
	{
		for (int f = 0; f < m_Factories; f++)
		{
			if (_1st_factft[f] != _2nd_factft[f])
				return true;
		}
		bool Flag = JudgeSameForTwoSol(_1st_Sol, _2nd_Sol);
		if (Flag)
			return true;
		else
			return false;
	}
}

float ComFun::JudgeSimilarityForTwoSol(const vector<vector<int>>& _1st_Sol, const vector<vector<int>>& _2nd_Sol)
{
	int simil = 0;
	for (int f = 0; f < m_Factories; f++)
	{
		int sizeofsol = min(_1st_Sol[f].size(), _2nd_Sol[f].size());
		for (int j = 0; j < sizeofsol; j++)
			if (_1st_Sol[f][j] == _2nd_Sol[f][j])
				simil++;
	}
	return (float)simil / m_Jobs;
}

void ComFun::GenerateRowByReduceDimen(const vector<vector<int>>& PreSeq, vector<int>& RefSeq)
{
	int Cnt = 0;
	for (int f = 0; f < m_Factories; f++)
	{
		int Size = PreSeq[f].size();
		for (int j = 0; j < Size; j++)
		{
			RefSeq[Cnt] = PreSeq[f][j];
			Cnt++;
		}
	}
}

bool ComFun::LocateJobInTwoDimArray(int Job, const vector<vector<int>>& TwoDimArray, int& FirstDimPos, int& SecondDimPos)
{
	int SizeofFirstDim = TwoDimArray.size();
	for (int f = 0; f < SizeofFirstDim; f++)
	{
		int SizeofSecondDim = TwoDimArray[f].size();
		for (int s = 0; s < SizeofSecondDim; s++)
		{
			if (TwoDimArray[f][s] == Job)
			{
				FirstDimPos = f;
				SecondDimPos = s;
				return true;
			}
		}
	}
	return false;
}

void ComFun::Check(const vector<vector<int>>& Sol, const vector<int>& FacTFT, int TotalTFT)
{
	int TTFT = 0;
	ofstream ofile;
	vector<bool> JobFlag(this->m_Jobs, false);
	int TotalJobs = 0;
	for (int fac = 0; fac < this->m_Factories; fac++)
	{
		for (unsigned int j = 0; j < Sol[fac].size(); j++)
			JobFlag[Sol[fac][j]] = true;
		TotalJobs += Sol[fac].size();

		int TFT = this->GetTotalFlowtime(Sol[fac]);
		if (TFT != FacTFT[fac])
		{
			ofile.open("C:\\error_log.txt", ios::app);
			ofile << "flowtime is Error" << this->m_Factories << "    " << this->m_Jobs << "     " << this->m_Machines << "     ";
			cout << "flowtime is Error in factory: " << fac << "  ThisTFT is: " << FacTFT[fac] << "  actualTFT is:" << TFT << endl;
			ofile.close();
		}
		TTFT += TFT;
	}
	if (TTFT != TotalTFT)
	{
		ofile.open("C:\\error_log.txt", ios::app);
		cout << "TotalTFT is error!" << endl;
		ofile << "TotalTFT is error!" << this->m_Factories << "    " << this->m_Jobs << "     " << this->m_Machines << "     ";
		ofile.close();
	}

	for (int j = 0; j < this->m_Jobs; j++)
		if (!JobFlag[j])
		{
			ofile.open("C:\\error_log.txt", ios::app);
			cout << "the job" << j << "is missing!" << endl;
			ofile << "the job is missing!" << this->m_Factories << "    " << this->m_Jobs << "     " << this->m_Machines << "     ";
			ofile.close();
		}
	if (TotalJobs > this->m_Jobs)
	{
		ofile.open("C:\\error_log.txt", ios::app);
		cout << "some jobs are repeated!" << endl;
		ofile << "some jobs are repeated!" << this->m_Factories << "    " << this->m_Jobs << "     " << this->m_Machines << "     ";
		ofile.close();
	}
	if (TotalJobs < this->m_Jobs)
	{
		ofile.open("C:\\error_log.txt", ios::app);
		cout << "some jobs are missing!" << endl;
		ofile << "some jobs are missing!" << this->m_Factories << "    " << this->m_Jobs << "     " << this->m_Machines << "     ";
		ofile.close();
	}
}

int ComFun::GetTotalFlowtime(const vector<int>& Seq)
{
	unsigned int Size_Seq = Seq.size();
	int Val_TF = 0;
	if (Size_Seq)
	{
		vector<vector<int>> Eij(Size_Seq, vector<int>(this->m_Machines));
		int Job = Seq[0];
		Eij[0][0] = this->m_Ptime[Job][0];
		for (int m = 1; m < this->m_Machines; m++)
			Eij[0][m] = Eij[0][m - 1] + this->m_Ptime[Job][m];
		for (unsigned int j = 1; j < Size_Seq; j++)
		{
			Job = Seq[j];
			Eij[j][0] = max(Eij[j - 1][0] + this->m_Ptime[Job][0], Eij[j - 1][1]);
			for (int m = 1; m < this->m_Machines - 1; m++)
			{
				Eij[j][m] = max(Eij[j][m - 1] + m_Ptime[Job][m], Eij[j - 1][m + 1]);
			}
			Eij[j][this->m_Machines - 1] = Eij[j][this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];
		}

		for (unsigned int j = 0; j < Size_Seq; j++)
			Val_TF += Eij[j][this->m_Machines - 1];
	}
	return Val_TF;
}

int ComFun::GetTotalFlowtime(const vector<int>& Seq, int& Span)
{
	unsigned int Size_Seq = Seq.size();
	int Val_TF = 0;
	Span = 0;
	if (Size_Seq)
	{
		vector<vector<int>> Eij(Size_Seq, vector<int>(this->m_Machines));
		int Job = Seq[0];
		Eij[0][0] = this->m_Ptime[Job][0];
		for (int m = 1; m < this->m_Machines; m++)
			Eij[0][m] = Eij[0][m - 1] + this->m_Ptime[Job][m];
		for (unsigned int j = 1; j < Size_Seq; j++)
		{
			Job = Seq[j];
			Eij[j][0] = max(Eij[j - 1][0] + this->m_Ptime[Job][0], Eij[j - 1][1]);
			for (int m = 1; m < this->m_Machines - 1; m++)
			{
				Eij[j][m] = max(Eij[j][m - 1] + m_Ptime[Job][m], Eij[j - 1][m + 1]);
			}
			Eij[j][this->m_Machines - 1] = Eij[j][this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];
		}

		for (unsigned int j = 0; j < Size_Seq; j++)
			Val_TF += Eij[j][this->m_Machines - 1];
		Span = Eij[Size_Seq - 1][m_Machines - 1];
	}
	return Val_TF;
}

int ComFun::GetTotalFlowtime(const vector<vector<int>>& Sol, vector<int>& FacTFT, vector<int>& FacSpan)
{
	int totaltft = 0;
	for (int f = 0; f < m_Factories; f++)
	{
		FacTFT[f] = GetTotalFlowtime(Sol[f], FacSpan[f]);
		totaltft += FacTFT[f];
	}
	return totaltft;
}

int ComFun::GetTotalFlowtime(const vector<vector<int>>& Sol, vector<int>& FacTFT, vector<vector<vector<int>>>& perEij)
{
	int totaltft = 0;
	for (int f = 0; f < m_Factories; f++)
	{
		FacTFT[f] = GetEij(Sol[f], perEij[f]);
		totaltft += FacTFT[f];
	}
	return totaltft;
}

int ComFun::GetEij(const vector<int> &Seq, vector<vector<int>>& Eij)
{
	unsigned int Size_Seq = Seq.size();
	Eij.resize(Size_Seq, vector<int>(m_Machines));
	int Val_TF = 0;
	if (Size_Seq)
	{
		int Job = Seq[0];
		Eij[0][0] = this->m_Ptime[Job][0];
		for (int m = 1; m < this->m_Machines; m++)
			Eij[0][m] = Eij[0][m - 1] + this->m_Ptime[Job][m];
		for (unsigned int j = 1; j < Size_Seq; j++)
		{
			Job = Seq[j];
			Eij[j][0] = max(Eij[j - 1][0] + this->m_Ptime[Job][0], Eij[j - 1][1]);
			for (int m = 1; m < this->m_Machines - 1; m++)
			{
				Eij[j][m] = max(Eij[j][m - 1] + m_Ptime[Job][m], Eij[j - 1][m + 1]);
			}
			Eij[j][this->m_Machines - 1] = Eij[j][this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];
		}

		for (unsigned int j = 0; j < Size_Seq; j++)
			Val_TF += Eij[j][this->m_Machines - 1];
	}
	return Val_TF;
}

int ComFun::GetEij(const vector<int> &Seq, vector<vector<int>>& Eij, int& MakeSpan)
{
	unsigned int Size_Seq = Seq.size();
	int Val_TF = 0;
	MakeSpan = 0;
	if (Size_Seq)
	{
		int Job = Seq[0];
		Eij[0][0] = this->m_Ptime[Job][0];
		for (int m = 1; m < this->m_Machines; m++)
			Eij[0][m] = Eij[0][m - 1] + this->m_Ptime[Job][m];
		for (unsigned int j = 1; j < Size_Seq; j++)
		{
			Job = Seq[j];
			Eij[j][0] = max(Eij[j - 1][0] + this->m_Ptime[Job][0], Eij[j - 1][1]);
			for (int m = 1; m < this->m_Machines - 1; m++)
			{
				Eij[j][m] = max(Eij[j][m - 1] + m_Ptime[Job][m], Eij[j - 1][m + 1]);
			}
			Eij[j][this->m_Machines - 1] = Eij[j][this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];
		}

		for (unsigned int j = 0; j < Size_Seq; j++)
			Val_TF += Eij[j][this->m_Machines - 1];
		MakeSpan = Eij[Size_Seq - 1][this->m_Machines - 1];
	}
	return Val_TF;
}

int ComFun::EraseAJob(vector<int>& Seq,int ErasePos, vector<vector<int>>& Eij)
{
	Seq.erase(Seq.begin() + ErasePos);
	int size_seq = Seq.size();
	Eij.resize(size_seq);
	if (ErasePos ==0)
	{
		GetEij(Seq, Eij);
	}
	else if (ErasePos < size_seq)
	{
		for (int j = ErasePos; j < size_seq; j++)
		{
			int Job = Seq[j];
			Eij[j][0] = max(Eij[j - 1][0] + this->m_Ptime[Job][0], Eij[j - 1][1]);
			for (int m = 1; m < this->m_Machines - 1; m++)
			{
				Eij[j][m] = max(Eij[j][m - 1] + m_Ptime[Job][m], Eij[j - 1][m + 1]);
			}
			Eij[j][this->m_Machines - 1] = Eij[j][this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];
		}
	}
	int totaltft = 0;
	for (int j = 0; j < size_seq; j++)
		totaltft += Eij[j][m_Machines - 1];
	return totaltft;
}

int ComFun::EraseAJob_2(const vector<vector<int>>& PreEij, vector<int>& Seq, int ErasePos, vector<vector<int>>& Eij)
{
	Seq.erase(Seq.begin() + ErasePos);
	int size_seq = Seq.size();
	if (ErasePos == 0)
	{
		GetEij(Seq, Eij);
	}
	else if (ErasePos <= size_seq)
	{
		for (int j = 0; j < ErasePos; j++)
			for (int m = 0; m < m_Machines; m++)
				Eij[j][m] = PreEij[j][m];
		for (int j = ErasePos; j < size_seq; j++)
		{
			int Job = Seq[j];
			Eij[j][0] = max(Eij[j - 1][0] + this->m_Ptime[Job][0], Eij[j - 1][1]);
			for (int m = 1; m < this->m_Machines - 1; m++)
			{
				Eij[j][m] = max(Eij[j][m - 1] + m_Ptime[Job][m], Eij[j - 1][m + 1]);
			}
			Eij[j][this->m_Machines - 1] = Eij[j][this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];
		}
	}
	int totaltft = 0;
	for (int j = 0; j < size_seq; j++)
		totaltft += Eij[j][m_Machines - 1];
	return totaltft;
}

int ComFun::GetTFTAfterPushAJob(const vector<int>& Seq, int PushJob, const vector<vector<int>>& Eij)
{
	unsigned int Size_Seq = Seq.size();
	if (Size_Seq == 0)
		return m_TotalPTime[PushJob];
	else
	{
		vector<int> Dij(m_Machines);
		Dij[0] = max(Eij[Size_Seq - 1][0] + this->m_Ptime[PushJob][0], Eij[Size_Seq - 1][1]);
		for (int m = 1; m < this->m_Machines - 1; m++)
		{
			Dij[m] = max(Dij[m - 1] + m_Ptime[PushJob][m], Eij[Size_Seq - 1][m + 1]);
		}
		Dij[this->m_Machines - 1] = Dij[this->m_Machines - 2] + this->m_Ptime[PushJob][this->m_Machines - 1];
		return Dij[m_Machines - 1];
	}
}

void ComFun::GetEijAfterInsertAJob(const vector<int>& Seq, int InsertJob, int InsertPos, vector<vector<int>>& Eij, int& Tft)
{
	unsigned int Size_Seq = Seq.size();
	Eij.resize(Size_Seq + 1, vector<int>(this->m_Machines));

	int Job = -1;
	if (InsertPos == 0)
	{
		Job = InsertJob;
		Eij[0][0] = this->m_Ptime[Job][0];
		for (int m = 1; m < this->m_Machines; m++)
			Eij[0][m] = Eij[0][m - 1] + this->m_Ptime[Job][m];

		for (unsigned int j = 1; j <= Size_Seq; j++)
		{
			Job = Seq[j - 1];
			Eij[j][0] = max(Eij[j - 1][0] + this->m_Ptime[Job][0], Eij[j - 1][1]);
			for (int m = 1; m < this->m_Machines - 1; m++)
			{
				Eij[j][m] = max(Eij[j][m - 1] + m_Ptime[Job][m], Eij[j - 1][m + 1]);
			}
			Eij[j][this->m_Machines - 1] = Eij[j][this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];
		}
	}
	else
	{
		Job = InsertJob;
		Eij[InsertPos][0] = max(Eij[InsertPos - 1][0] + this->m_Ptime[Job][0], Eij[InsertPos - 1][1]);
		for (int m = 1; m < this->m_Machines - 1; m++)
		{
			Eij[InsertPos][m] = max(Eij[InsertPos][m - 1] + m_Ptime[Job][m], Eij[InsertPos - 1][m + 1]);
		}
		Eij[InsertPos][this->m_Machines - 1] = Eij[InsertPos][this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];

		for (unsigned int j = InsertPos + 1; j <= Size_Seq; j++)
		{
			Job = Seq[j - 1];
			Eij[j][0] = max(Eij[j - 1][0] + this->m_Ptime[Job][0], Eij[j - 1][1]);
			for (int m = 1; m < this->m_Machines - 1; m++)
			{
				Eij[j][m] = max(Eij[j][m - 1] + m_Ptime[Job][m], Eij[j - 1][m + 1]);
			}
			Eij[j][this->m_Machines - 1] = Eij[j][this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];
		}
	}

	Tft = 0;
	for (unsigned int j = 0; j <= Size_Seq; j++)
		Tft += Eij[j][this->m_Machines - 1];
}

int ComFun::GetTFTAfterInsertAJob(const vector<int>& Seq, int InsertJob, int InsertPos, const vector<vector<int>>& PreEij)
{
	unsigned int Size_Seq = Seq.size();
	vector<int> Dij(m_Machines);
	int Val_TF = 0;
	if (InsertPos == 0)
	{
		int Job = InsertJob;
		Dij[0] = this->m_Ptime[Job][0];
		for (int m = 1; m < this->m_Machines; m++)
			Dij[m] = Dij[m - 1] + this->m_Ptime[Job][m];
		Val_TF += Dij[m_Machines - 1];

		for (unsigned int j = 1; j <= Size_Seq; j++)
		{
			Job = Seq[j - 1];
			Dij[0] = max(Dij[0] + this->m_Ptime[Job][0], Dij[1]);
			for (int m = 1; m < this->m_Machines - 1; m++)
			{
				Dij[m] = max(Dij[m - 1] + m_Ptime[Job][m], Dij[m + 1]);
			}
			Dij[this->m_Machines - 1] = Dij[this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];
			Val_TF += Dij[this->m_Machines - 1];
		}
	}
	else
	{
		for (int p = 0; p < InsertPos; p++)
			Val_TF += PreEij[p][m_Machines - 1];
		int Job = InsertJob;
		Dij[0] = max(PreEij[InsertPos - 1][0] + this->m_Ptime[Job][0], PreEij[InsertPos - 1][1]);
		for (int m = 1; m < this->m_Machines - 1; m++)
		{
			Dij[m] = max(Dij[m - 1] + m_Ptime[Job][m], PreEij[InsertPos - 1][m + 1]);
		}
		Dij[this->m_Machines - 1] = Dij[this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];
		Val_TF += Dij[this->m_Machines - 1];
		for (int j = 1; j < Size_Seq + 1 - InsertPos; j++)
		{
			Job = Seq[InsertPos + j - 1];
			Dij[0] = max(Dij[0] + this->m_Ptime[Job][0], Dij[1]);
			for (int m = 1; m < m_Machines - 1; m++)
			{
				Dij[m] = max(Dij[m - 1] + m_Ptime[Job][m], Dij[m + 1]);
			}
			Dij[m_Machines - 1] = Dij[m_Machines - 2] + m_Ptime[Job][m_Machines - 1];
			Val_TF += Dij[this->m_Machines - 1];
		}
	}
	return Val_TF;
}

int ComFun::GetTFTAfterSwapAJob(const vector<int>& Seq, int SwapPos, int SwapJob, const vector<vector<int>>& Eij)
{
	int Size_Seq = Seq.size();
	vector<int> Dj(m_Machines);
	int factft = 0;
	if (SwapPos == 0)
	{
		int Job = SwapJob;
		Dj[0] = this->m_Ptime[Job][0];
		for (int m = 1; m < this->m_Machines; m++)
			Dj[m] = Dj[m - 1] + this->m_Ptime[Job][m];
		factft += Dj[m_Machines - 1];
		for (unsigned int j = 1; j < Size_Seq; j++)
		{
			Job = Seq[j];
			Dj[0] = max(Dj[0] + this->m_Ptime[Job][0], Dj[1]);
			for (int m = 1; m < this->m_Machines - 1; m++)
			{
				Dj[m] = max(Dj[m - 1] + m_Ptime[Job][m], Dj[m + 1]);
			}
			Dj[this->m_Machines - 1] = Dj[this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];
			factft += Dj[m_Machines - 1];
		}
	}
	else
	{
		for (int j = 0; j < SwapPos; j++)
			factft += Eij[j][m_Machines - 1];
		int Job = SwapJob;
		Dj[0] = max(Eij[SwapPos - 1][0] + this->m_Ptime[Job][0], Eij[SwapPos - 1][1]);
		for (int m = 1; m < this->m_Machines - 1; m++)
		{
			Dj[m] = max(Dj[m - 1] + m_Ptime[Job][m], Eij[SwapPos - 1][m + 1]);
		}
		Dj[this->m_Machines - 1] = Dj[this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];
		factft += Dj[m_Machines - 1];
		for (int j = SwapPos + 1; j < Size_Seq; j++)
		{
			Job = Seq[j];
			Dj[0] = max(Dj[0] + this->m_Ptime[Job][0], Dj[1]);
			for (int m = 1; m < this->m_Machines - 1; m++)
			{
				Dj[m] = max(Dj[m - 1] + m_Ptime[Job][m], Dj[m + 1]);
			}
			Dj[this->m_Machines - 1] = Dj[this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];
			factft += Dj[m_Machines - 1];
		}
	}
	return factft;
}

void ComFun::GetEijAfterSwapAJob(const vector<int>& Seq, int SwapPos, int SwapJob, vector<vector<int>>& Eij, int& Tft)
{
	int Size_Seq = Seq.size();
	if (SwapPos == 0)
	{
		int Job = SwapJob;
		Eij[0][0] = this->m_Ptime[Job][0];
		for (int m = 1; m < this->m_Machines; m++)
			Eij[0][m] = Eij[0][m - 1] + this->m_Ptime[Job][m];
		for (unsigned int j = 1; j < Size_Seq; j++)
		{
			Job = Seq[j];
			Eij[j][0] = max(Eij[j - 1][0] + this->m_Ptime[Job][0], Eij[j - 1][1]);
			for (int m = 1; m < this->m_Machines - 1; m++)
			{
				Eij[j][m] = max(Eij[j][m - 1] + m_Ptime[Job][m], Eij[j - 1][m + 1]);
			}
			Eij[j][this->m_Machines - 1] = Eij[j][this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];
		}
	}
	else
	{
		int Job = SwapJob;
		Eij[SwapPos][0] = max(Eij[SwapPos - 1][0] + this->m_Ptime[Job][0], Eij[SwapPos - 1][1]);
		for (int m = 1; m < this->m_Machines - 1; m++)
		{
			Eij[SwapPos][m] = max(Eij[SwapPos][m - 1] + m_Ptime[Job][m], Eij[SwapPos - 1][m + 1]);
		}
		Eij[SwapPos][this->m_Machines - 1] = Eij[SwapPos][this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];

		for (int j = SwapPos + 1; j < Size_Seq; j++)
		{
			Job = Seq[j];
			Eij[j][0] = max(Eij[j - 1][0] + this->m_Ptime[Job][0], Eij[j - 1][1]);
			for (int m = 1; m < this->m_Machines - 1; m++)
			{
				Eij[j][m] = max(Eij[j][m - 1] + m_Ptime[Job][m], Eij[j - 1][m + 1]);
			}
			Eij[j][this->m_Machines - 1] = Eij[j][this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];
		}
	}

	Tft = 0;
	for (int j = 0; j < Size_Seq; j++) {
		Tft += Eij[j][this->m_Machines - 1];
	}
}

int ComFun::GetTFTAfterSwapJobsInFac(vector<int>& Seq, int ActivePos, int PassivePos, const vector<vector<int>>& Eij)
{
	int Size_Seq = Seq.size();
	swap(Seq[ActivePos], Seq[PassivePos]);
	int FrontPos = min(ActivePos, PassivePos);

	int factft = 0;
	vector<int> Dj(m_Machines);
	if (FrontPos == 0)
	{
		int Job = Seq[0];
		Dj[0] = this->m_Ptime[Job][0];
		for (int m = 1; m < this->m_Machines; m++)
			Dj[m] = Dj[m - 1] + this->m_Ptime[Job][m];
		factft += Dj[m_Machines - 1];
		for (unsigned int j = 1; j < Size_Seq; j++)
		{
			Job = Seq[j];
			Dj[0] = max(Dj[0] + this->m_Ptime[Job][0], Dj[1]);
			for (int m = 1; m < this->m_Machines - 1; m++)
			{
				Dj[m] = max(Dj[m - 1] + m_Ptime[Job][m], Dj[m + 1]);
			}
			Dj[this->m_Machines - 1] = Dj[this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];
			factft += Dj[m_Machines - 1];
		}
	}
	else
	{
		for (int j = 0; j < FrontPos; j++)
			factft += Eij[j][m_Machines - 1];
		int Job = Seq[FrontPos];
		Dj[0] = max(Eij[FrontPos - 1][0] + this->m_Ptime[Job][0], Eij[FrontPos - 1][1]);
		for (int m = 1; m < this->m_Machines - 1; m++)
		{
			Dj[m] = max(Dj[m - 1] + m_Ptime[Job][m], Eij[FrontPos - 1][m + 1]);
		}
		Dj[this->m_Machines - 1] = Dj[this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];
		factft += Dj[m_Machines - 1];
		for (int j = FrontPos + 1; j < Size_Seq; j++)
		{
			Job = Seq[j];
			Dj[0] = max(Dj[0] + this->m_Ptime[Job][0], Dj[1]);
			for (int m = 1; m < this->m_Machines - 1; m++)
			{
				Dj[m] = max(Dj[m - 1] + m_Ptime[Job][m], Dj[m + 1]);
			}
			Dj[this->m_Machines - 1] = Dj[this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];
			factft += Dj[m_Machines - 1];
		}
	}

	swap(Seq[ActivePos], Seq[PassivePos]);
	return factft;
}

int ComFun::GetMinTFTAfterInsertSingleFac(const vector<int>& PreSeq, int InJob, int& BestPos, const vector<vector<int>>& Eij)
{
	int Size_Seq = PreSeq.size();
	vector<int> TFT_Array(Size_Seq + 1);
	for (int j = 0; j <= Size_Seq; j++)
	{
		int Val_TF = GetTFTAfterInsertAJob(PreSeq, InJob, j, Eij);
		TFT_Array[j] = Val_TF;
	}
	BestPos = min_element(TFT_Array.begin(), TFT_Array.end()) - TFT_Array.begin();
	return TFT_Array[BestPos];
}

int ComFun::GetTFAfterInsertAJobInsideFactory(const vector<int>& Seq, int MovePos, int InsertPos, const vector<vector<int>>& Eij)
{
	unsigned int Size_Seq = Seq.size();
	int factft = 0;
	if (InsertPos < MovePos)
	{
		vector<int> Dj(m_Machines);
		if (InsertPos == 0)
		{
			int Job = Seq[MovePos];
			Dj[0] = this->m_Ptime[Job][0];
			for (int m = 1; m < this->m_Machines; m++)
				Dj[m] = Dj[m - 1] + this->m_Ptime[Job][m];
			factft += Dj[m_Machines - 1];
		}
		else
		{
			for (int j = 0; j < InsertPos; j++)
				factft += Eij[j][m_Machines - 1];

			int Job = Seq[MovePos];
			Dj[0] = max(Eij[InsertPos - 1][0] + this->m_Ptime[Job][0], Eij[InsertPos - 1][1]);
			for (int m = 1; m < this->m_Machines - 1; m++)
			{
				Dj[m] = max(Dj[m - 1] + m_Ptime[Job][m], Eij[InsertPos - 1][m + 1]);
			}
			Dj[this->m_Machines - 1] = Dj[this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];
			factft += Dj[m_Machines - 1];
		}

		for (unsigned int j = InsertPos + 1; j < Size_Seq; j++)
		{
			int Job = -1;
			if (j <= MovePos)
				Job = Seq[j - 1];
			else
				Job = Seq[j];
			Dj[0] = max(Dj[0] + this->m_Ptime[Job][0], Dj[1]);
			for (int m = 1; m < this->m_Machines - 1; m++)
			{
				Dj[m] = max(Dj[m - 1] + m_Ptime[Job][m], Dj[m + 1]);
			}
			Dj[this->m_Machines - 1] = Dj[this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];
			factft += Dj[this->m_Machines - 1];
		}
	}
	else
	{
		vector<int> Dj(m_Machines);
		if (MovePos == 0)
		{
			int Job = Seq[1];
			Dj[0] = this->m_Ptime[Job][0];
			for (int m = 1; m < this->m_Machines; m++)
				Dj[m] = Dj[m - 1] + this->m_Ptime[Job][m];
			factft += Dj[m_Machines - 1];
		}
		else
		{
			for (int j = 0; j < MovePos; j++)
				factft += Eij[j][m_Machines - 1];
			int Job = Seq[MovePos + 1];
			Dj[0] = max(Eij[MovePos - 1][0] + this->m_Ptime[Job][0], Eij[MovePos - 1][1]);
			for (int m = 1; m < this->m_Machines - 1; m++)
			{
				Dj[m] = max(Dj[m - 1] + m_Ptime[Job][m], Eij[MovePos - 1][m + 1]);
			}
			Dj[this->m_Machines - 1] = Dj[this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];
			factft += Dj[m_Machines - 1];
		}

		for (int j = MovePos + 1; j < Size_Seq; j++)
		{
			int Job = -1;
			if (j < InsertPos - 1)
				Job = Seq[j + 1];
			else if (j == InsertPos - 1)
				Job = Seq[MovePos];
			else
				Job = Seq[j];
			Dj[0] = max(Dj[0] + this->m_Ptime[Job][0], Dj[1]);
			for (int m = 1; m < this->m_Machines - 1; m++)
			{
				Dj[m] = max(Dj[m - 1] + m_Ptime[Job][m], Dj[m + 1]);
			}
			Dj[this->m_Machines - 1] = Dj[this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];
			factft += Dj[m_Machines - 1];
		}
	}
	return factft;
}

int ComFun::FindBestInsertPosInsideFactory(vector<int>& Seq, int MovePos, const vector<vector<int>>& Eij, int& BestPos)
{
	int Size_Seq = Seq.size();
	vector<int> TFT_Array(Size_Seq + 1, INT_MAX);
	for (int j = 0; j <= Size_Seq; j++)
	{
		if (j == MovePos || j == MovePos + 1)
			continue;
		TFT_Array[j] = GetTFAfterInsertAJobInsideFactory(Seq, MovePos, j, Eij);
	}
	BestPos = min_element(TFT_Array.begin(), TFT_Array.end()) - TFT_Array.begin();
	return TFT_Array[BestPos];
}

int ComFun::FindBestInsertPosFacByMinIncreasedTFT(const vector<vector<int>>& Sol, const vector<int>& FacTFT,
	const vector<vector<vector<int>>>& perEij, int InJob, int& BestFac, int& BestPos)
{
	vector<int> IncreasedTFT_Array(this->m_Factories);
	vector<int> BestPos_Array(this->m_Factories, -1);
	for (int f = 0; f < this->m_Factories; f++)
	{
		int resulting_tft = GetMinTFTAfterInsertSingleFac(Sol[f], InJob, BestPos_Array[f], perEij[f]);
		IncreasedTFT_Array[f] = resulting_tft - FacTFT[f];
	}
	BestFac = min_element(IncreasedTFT_Array.begin(), IncreasedTFT_Array.end()) - IncreasedTFT_Array.begin();
	BestPos = BestPos_Array[BestFac];
	return (IncreasedTFT_Array[BestFac] + FacTFT[BestFac]);
}

int ComFun::GetBestSwapPosInOriginalFactory(vector<int>& PreSeq, int SwapPos, 
	vector<vector<int>> Eij, int PreTFT, int& BestSwapPos)
{
	int Size_PreSeq = PreSeq.size();
	int MinTFT = INT_MAX;
	for (int j = 0; j < SwapPos; j++)
	{
		int TestTFT = this->GetTFTAfterSwapJobsInFac(PreSeq, SwapPos, j, Eij);

		if (TestTFT < MinTFT)
		{
			MinTFT = TestTFT;
			BestSwapPos = j;
		}
	}
	return PreTFT - MinTFT;
}

int ComFun::GetBestSwapPosAfterSwapAmongTwoFactory(const vector<int>& ActiveSeq, const vector<vector<int>>& ActiveEij, int PreTFT_ActiveSeq, int SwapPos,
	const vector<int>& PassiveSeq, const vector<vector<int>>& PassiveEij, int PreTFT_PassiveSeq, int& BestSwapPos)
{
	int Size_ActiveSeq = ActiveSeq.size();
	int Size_PassiveSeq = PassiveSeq.size();

	int MINTFT_TwoFac = INT_MAX;
	for (int j = 0; j < Size_PassiveSeq; j++)
	{
		int TestTFT_FirstFac = this->GetTFTAfterSwapAJob(ActiveSeq, SwapPos, PassiveSeq[j], ActiveEij);
		int TestTFT_SecondFac = this->GetTFTAfterSwapAJob(PassiveSeq, j, ActiveSeq[SwapPos], PassiveEij);
		int TestTFT_TwoFac = TestTFT_FirstFac + TestTFT_SecondFac;
		if (TestTFT_TwoFac < MINTFT_TwoFac)
		{
			BestSwapPos = j;
			MINTFT_TwoFac = TestTFT_TwoFac;
		}
	}
	return (PreTFT_ActiveSeq + PreTFT_PassiveSeq - MINTFT_TwoFac);
}

