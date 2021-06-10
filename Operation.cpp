#include "Operation.h"

Operation::Operation()
{

}

Operation::~Operation()
{

}

//构造解
void Operation::PushAJobToSeq(int Job, vector<int>& Seq, vector<vector<int>>& Eij, int& totaltft)
{
	unsigned int Size_Seq = Seq.size();
	Eij.resize(Size_Seq + 1, vector<int>(this->m_Machines));//-----+1,二维的vector扩容，原数据不变
	if (Size_Seq == 0)
	{
		Eij[0][0] = this->m_Ptime[Job][0];
		for (int m = 1; m < this->m_Machines; m++)
			Eij[0][m] = Eij[0][m - 1] + this->m_Ptime[Job][m];
	}
	else
	{
		Eij[Size_Seq][0] = max(Eij[Size_Seq - 1][0] + this->m_Ptime[Job][0], Eij[Size_Seq - 1][1]);
		for (int m = 1; m < this->m_Machines - 1; m++)
		{
			Eij[Size_Seq][m] = max(Eij[Size_Seq][m - 1] + m_Ptime[Job][m], Eij[Size_Seq - 1][m + 1]);
		}
		Eij[Size_Seq][this->m_Machines - 1] = Eij[Size_Seq][this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];
	}
	totaltft += Eij[Size_Seq][m_Machines - 1];
	Seq.push_back(Job);
}

void Operation::InsertAJobToSeq(int InsertJob, int InsertPos, vector<int>& Seq, vector<vector<int>>& Eij, int& tft)
{
	unsigned int Size_Seq = Seq.size();
	Eij.resize(Size_Seq + 1, vector<int>(this->m_Machines));//-----+1,二维的vector扩容，原数据不变

	int Job = -1;
	if (InsertPos == 0)
	{
		//插入位置为0，全部都重新算
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
		//插入位置不是0时，则从位置0到（插入位置-1）间的 Eij 不用再计算

		Job = InsertJob;//-------插入位置处
		Eij[InsertPos][0] = max(Eij[InsertPos - 1][0] + this->m_Ptime[Job][0], Eij[InsertPos - 1][1]);
		for (int m = 1; m < this->m_Machines - 1; m++)
		{
			Eij[InsertPos][m] = max(Eij[InsertPos][m - 1] + m_Ptime[Job][m], Eij[InsertPos - 1][m + 1]);
		}
		Eij[InsertPos][this->m_Machines - 1] = Eij[InsertPos][this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];

		//--------------------插入位置之后
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

	tft = 0;
	for (unsigned int j = 0; j <= Size_Seq; j++)
		tft += Eij[j][this->m_Machines - 1];

	Seq.insert(Seq.begin() + InsertPos, InsertJob);
}

void Operation::InsertAJobToItsSeq(int MovePos, int InsertPos, vector<int>& Seq, vector<vector<int>>& Eij, int& tft)
{
	//InsertPos!=MovePos, MovePos+1
	unsigned int Size_Seq = Seq.size();
	if (MovePos > InsertPos)//其一，插入位置在原位置其之前
	{
		//-------插入位置处
		if (InsertPos == 0)
		{
			int Job = Seq[MovePos];
			Eij[0][0] = this->m_Ptime[Job][0];
			for (int m = 1; m < this->m_Machines; m++)
				Eij[0][m] = Eij[0][m - 1] + this->m_Ptime[Job][m];
		}
		else
		{
			int Job = Seq[MovePos];
			Eij[InsertPos][0] = max(Eij[InsertPos - 1][0] + this->m_Ptime[Job][0], Eij[InsertPos - 1][1]);
			for (int m = 1; m < this->m_Machines - 1; m++)
			{
				Eij[InsertPos][m] = max(Eij[InsertPos][m - 1] + m_Ptime[Job][m], Eij[InsertPos - 1][m + 1]);
			}
			Eij[InsertPos][this->m_Machines - 1] = Eij[InsertPos][this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];
		}
		for (unsigned int j = InsertPos + 1; j < Size_Seq; j++)
		{
			int Job = -1;
			if (j <= MovePos)
				Job = Seq[j - 1];//两点之间
			else
				Job = Seq[j];//移动点之后
			Eij[j][0] = max(Eij[j - 1][0] + this->m_Ptime[Job][0], Eij[j - 1][1]);
			for (int m = 1; m < this->m_Machines - 1; m++)
			{
				Eij[j][m] = max(Eij[j][m - 1] + m_Ptime[Job][m], Eij[j - 1][m + 1]);
			}
			Eij[j][this->m_Machines - 1] = Eij[j][this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];
		}

		Seq.insert(Seq.begin() + InsertPos, Seq[MovePos]);
		Seq.erase(Seq.begin() + MovePos + 1);
	}
	else//其二，插入位置在其之后
	{
		if (MovePos == 0)
		{
			int Job = Seq[MovePos + 1];
			Eij[MovePos][0] = this->m_Ptime[Job][0];
			for (int m = 1; m < this->m_Machines; m++)
				Eij[MovePos][m] = Eij[MovePos][m - 1] + this->m_Ptime[Job][m];
		}
		else
		{
			int Job = Seq[MovePos + 1];
			Eij[MovePos][0] = max(Eij[MovePos - 1][0] + this->m_Ptime[Job][0], Eij[MovePos - 1][1]);
			for (int m = 1; m < this->m_Machines - 1; m++)
			{
				Eij[MovePos][m] = max(Eij[MovePos][m - 1] + m_Ptime[Job][m], Eij[MovePos - 1][m + 1]);
			}
			Eij[MovePos][this->m_Machines - 1] = Eij[MovePos][this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];
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
			Eij[j][0] = max(Eij[j - 1][0] + this->m_Ptime[Job][0], Eij[j - 1][1]);
			for (int m = 1; m < this->m_Machines - 1; m++)
			{
				Eij[j][m] = max(Eij[j][m - 1] + m_Ptime[Job][m], Eij[j - 1][m + 1]);
			}
			Eij[j][this->m_Machines - 1] = Eij[j][this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];
		}

		Seq.insert(Seq.begin() + InsertPos, Seq[MovePos]);
		Seq.erase(Seq.begin() + MovePos);
	}
	tft = 0;
	for (unsigned int j = 0; j < Size_Seq; j++)
		tft += Eij[j][this->m_Machines - 1];
}


//邻域操作
void Operation::SwapTwoJobsInsideFactory(int pos1, int pos2, vector<int>& Seq, vector<vector<int>>& Eij, int& Tft, 
	int& TotalTFT)
{
	int Size_Seq = Seq.size();
	swap(Seq[pos1], Seq[pos2]);

	int FrontPos = min(pos1, pos2);

	if (FrontPos == 0)
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
	}
	else
	{
		int Job = Seq[FrontPos];
		Eij[FrontPos][0] = max(Eij[FrontPos - 1][0] + this->m_Ptime[Job][0], Eij[FrontPos - 1][1]);
		for (int m = 1; m < this->m_Machines - 1; m++)
		{
			Eij[FrontPos][m] = max(Eij[FrontPos][m - 1] + m_Ptime[Job][m], Eij[FrontPos - 1][m + 1]);
		}
		Eij[FrontPos][this->m_Machines - 1] = Eij[FrontPos][this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];

		for (int j = FrontPos + 1; j < Size_Seq; j++)
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

	TotalTFT -= Tft;
	Tft = 0;
	for (int j = 0; j < Size_Seq; j++) {
		Tft += Eij[j][this->m_Machines - 1];
	}
	TotalTFT += Tft;
}

void Operation::SwapTwoJobsInsideFactory(int pos1, int pos2, vector<int>& Seq, vector<vector<int>>& Eij)
{
	int Size_Seq = Seq.size();
	swap(Seq[pos1], Seq[pos2]);

	int FrontPos = min(pos1, pos2);

	if (FrontPos == 0)
	{
		int Job = Seq[0];
		Eij[0][0] = this->m_Ptime[Job][0];
		for (int m = 1; m < this->m_Machines; m++)
			Eij[0][m] = Eij[0][m - 1] + this->m_Ptime[Job][m];
		for (int j = 1; j < Size_Seq; j++)
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
		int Job = Seq[FrontPos];
		Eij[FrontPos][0] = max(Eij[FrontPos - 1][0] + this->m_Ptime[Job][0], Eij[FrontPos - 1][1]);
		for (int m = 1; m < this->m_Machines - 1; m++)
		{
			Eij[FrontPos][m] = max(Eij[FrontPos][m - 1] + m_Ptime[Job][m], Eij[FrontPos - 1][m + 1]);
		}
		Eij[FrontPos][this->m_Machines - 1] = Eij[FrontPos][this->m_Machines - 2] + this->m_Ptime[Job][this->m_Machines - 1];

		for (int j = FrontPos + 1; j < Size_Seq; j++)
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
}

void Operation::SwapTwoJobsBetweenFactories(int fac1, int fac2, int pos1, int pos2,
	vector<vector<int>>& Sol, vector<vector<vector<int>>>& perEij, vector<int>& FacTFT, int& Totaltft)
{
	GetEijAfterSwapAJob(Sol[fac1], pos1, Sol[fac2][pos2], perEij[fac1], FacTFT[fac1]);
	GetEijAfterSwapAJob(Sol[fac2], pos2, Sol[fac1][pos1], perEij[fac2], FacTFT[fac2]);
	swap(Sol[fac1][pos1], Sol[fac2][pos2]);
	Totaltft = accumulate(FacTFT.begin(), FacTFT.end(), 0);
}


bool Operation::CarryBestSwapByMinGlobalTFT(int MoveFac, int MovePos,
	vector<vector<int>>& Sol, vector<int>& FacTFT, int& TotalTFT, vector<vector<vector<int>>>& perEij)
{
	bool ImpFlag = false;
	vector<int> BestSwapPos(m_Factories);
	vector<int> MinDescreasedTFT(m_Factories);//原值-移动后的值，即正值是有优化的

	MinDescreasedTFT[MoveFac] = GetBestSwapPosInOriginalFactory(Sol[MoveFac], MovePos, perEij[MoveFac], FacTFT[MoveFac], BestSwapPos[MoveFac]);

	for (int f = 0; f < MoveFac; f++)
	{
		MinDescreasedTFT[f] = GetBestSwapPosAfterSwapAmongTwoFactory(Sol[MoveFac], perEij[MoveFac], FacTFT[MoveFac], MovePos,
			Sol[f], perEij[f], FacTFT[f], BestSwapPos[f]);
	}

	int BestFac = max_element(MinDescreasedTFT.begin(), MinDescreasedTFT.end()) - MinDescreasedTFT.begin();
	if (MinDescreasedTFT[BestFac] > 0)
	{
		if (MoveFac == BestFac)
			SwapTwoJobsInsideFactory(MovePos, BestSwapPos[BestFac],
				Sol[MoveFac], perEij[MoveFac], FacTFT[MoveFac], TotalTFT);
		else
			SwapTwoJobsBetweenFactories(MoveFac, BestFac, MovePos, BestSwapPos[BestFac],
				Sol, perEij, FacTFT, TotalTFT);
		ImpFlag = true;
	}
	return ImpFlag;
}

bool Operation::CarryBestInsertByMinGlobalTFT(int MoveFac, int MovePos,
	vector<vector<int>>& Sol, vector<int>& FacTFT, int& TotalTFT, vector<vector<vector<int>>>& perEij)
{
	bool ImpFlag = false;
	int ExtJob = Sol[MoveFac][MovePos];
	//先求原工厂
	vector<vector<int>> NewEij(perEij[MoveFac]);
	//static vector<vector<int>> NewEij;
	//NewEij = perEij[MoveFac];
	int NewTFT_ExtFac = EraseAJob(Sol[MoveFac], MovePos, NewEij);

	vector<int> BestInsertPos(this->m_Factories, -1);
	vector<int> ResultMinTotalTFT(this->m_Factories, -1);
	vector<int> TestMinTFT(this->m_Factories, -1);

	TestMinTFT[MoveFac] = this->GetMinTFTAfterInsertSingleFac(Sol[MoveFac], ExtJob, BestInsertPos[MoveFac], NewEij);
	ResultMinTotalTFT[MoveFac] = TotalTFT - FacTFT[MoveFac] + TestMinTFT[MoveFac];

	for (int f = 0; f < m_Factories; f++)
	{
		if (f == MoveFac)
			continue;
		TestMinTFT[f] = this->GetMinTFTAfterInsertSingleFac(Sol[f], ExtJob, BestInsertPos[f], perEij[f]);
		ResultMinTotalTFT[f] = TotalTFT - FacTFT[MoveFac] + NewTFT_ExtFac - FacTFT[f] + TestMinTFT[f];
	}

	int BestFac = min_element(ResultMinTotalTFT.begin(), ResultMinTotalTFT.end()) - ResultMinTotalTFT.begin();
	if (ResultMinTotalTFT[BestFac] < TotalTFT)
	{
		if (BestFac == MoveFac)
		{
			TotalTFT -= FacTFT[BestFac];
			FacTFT[BestFac] = TestMinTFT[BestFac];
			GetEijAfterInsertAJob(Sol[BestFac], ExtJob, BestInsertPos[BestFac], NewEij, TestMinTFT[BestFac]);
			perEij[BestFac] = NewEij;
			TotalTFT += FacTFT[BestFac];
		}
		else
		{
			TotalTFT = TotalTFT - FacTFT[MoveFac] - FacTFT[BestFac];
			FacTFT[MoveFac] = NewTFT_ExtFac;
			perEij[MoveFac] = NewEij;
			FacTFT[BestFac] = TestMinTFT[BestFac];
			GetEijAfterInsertAJob(Sol[BestFac], ExtJob, BestInsertPos[BestFac], perEij[BestFac], TestMinTFT[BestFac]);
			TotalTFT = TotalTFT + FacTFT[MoveFac] + FacTFT[BestFac];
		}
		Sol[BestFac].insert(Sol[BestFac].begin() + BestInsertPos[BestFac], ExtJob);
		return true;
	}
	else
	{
		Sol[MoveFac].insert(Sol[MoveFac].begin() + MovePos, ExtJob);
	}
	return false;
}

bool Operation::CarryBestInsertByMinGlobalTFT(int MoveFac, int MovePos,
	vector<vector<int>>& Sol, vector<int>& FacTFT, vector<int>& FacSpan, int& TotalTFT, vector<vector<vector<int>>>& perEij)
{
	bool ImpFlag = false;
	int ExtJob = Sol[MoveFac][MovePos];
	//先求原工厂
	vector<vector<int>> NewEij(perEij[MoveFac]);
	int NewTFT_ExtFac = EraseAJob(Sol[MoveFac], MovePos, NewEij);

	vector<int> BestInsertPos(this->m_Factories, -1);
	vector<int> ResultMinTotalTFT(this->m_Factories, -1);
	vector<int> TestMinTFT(this->m_Factories, -1);

	TestMinTFT[MoveFac] = this->GetMinTFTAfterInsertSingleFac(Sol[MoveFac], ExtJob, BestInsertPos[MoveFac], NewEij);
	ResultMinTotalTFT[MoveFac] = TotalTFT - FacTFT[MoveFac] + TestMinTFT[MoveFac];

	for (int f = 0; f < m_Factories; f++)
	{
		if (f == MoveFac)
			continue;
		TestMinTFT[f] = this->GetMinTFTAfterInsertSingleFac(Sol[f], ExtJob, BestInsertPos[f], perEij[f]);
		ResultMinTotalTFT[f] = TotalTFT - FacTFT[MoveFac] + NewTFT_ExtFac - FacTFT[f] + TestMinTFT[f];
	}

	int BestFac = min_element(ResultMinTotalTFT.begin(), ResultMinTotalTFT.end()) - ResultMinTotalTFT.begin();
	if (ResultMinTotalTFT[BestFac] < TotalTFT)
	{
		if (BestFac == MoveFac)
		{
			TotalTFT -= FacTFT[BestFac];
			FacTFT[BestFac] = TestMinTFT[BestFac];
			GetEijAfterInsertAJob(Sol[BestFac], ExtJob, BestInsertPos[BestFac], NewEij, TestMinTFT[BestFac]);
			perEij[BestFac] = NewEij;
			TotalTFT += FacTFT[BestFac];
		}
		else
		{
			TotalTFT = TotalTFT - FacTFT[MoveFac] - FacTFT[BestFac];
			FacTFT[MoveFac] = NewTFT_ExtFac;
			perEij[MoveFac] = NewEij;
			FacTFT[BestFac] = TestMinTFT[BestFac];
			GetEijAfterInsertAJob(Sol[BestFac], ExtJob, BestInsertPos[BestFac], perEij[BestFac], TestMinTFT[BestFac]);
			TotalTFT = TotalTFT + FacTFT[MoveFac] + FacTFT[BestFac];
		}
		Sol[BestFac].insert(Sol[BestFac].begin() + BestInsertPos[BestFac], ExtJob);
		FacSpan[BestFac] = perEij[BestFac][Sol[BestFac].size() - 1][m_Machines - 1];
		FacSpan[MoveFac] = perEij[MoveFac][Sol[MoveFac].size() - 1][m_Machines - 1];
		return true;
	}
	else
	{
		Sol[MoveFac].insert(Sol[MoveFac].begin() + MovePos, ExtJob);
	}
	return false;
}


bool Operation::CarryBestInsertByMinGlobalTFT_randfactory(int MoveFac, int MovePos,
	vector<vector<int>>& Sol, vector<int>& FacTFT, int& TotalTFT, vector<vector<vector<int>>>& perEij)
{
	vector<int> ref_factory;
	GetOrderSeq(m_Factories, ref_factory);
	random_shuffle(ref_factory.begin(), ref_factory.end());

	bool ImpFlag = false;
	int ExtJob = Sol[MoveFac][MovePos];

	vector<vector<int>> NewEij(perEij[MoveFac]);
	int NewTFT_ExtFac = EraseAJob(Sol[MoveFac], MovePos, NewEij);

	//不包含原工厂
	vector<int> BestInsertPos(this->m_Factories, -1);
	vector<int> ResultMinTotalTFT(this->m_Factories, INT_MAX);
	vector<int> TestMinTFT(this->m_Factories, INT_MAX);

	for (int tf = 0; tf < m_Factories; tf++)
	{
		int f = ref_factory[tf];
		if (f == MoveFac)
			continue;
		else
		{
			TestMinTFT[f] = this->GetMinTFTAfterInsertSingleFac(Sol[f], ExtJob, BestInsertPos[f], perEij[f]);
			ResultMinTotalTFT[f] = TotalTFT - FacTFT[MoveFac] + NewTFT_ExtFac - FacTFT[f] + TestMinTFT[f];
		}
	}

	int BestFac = min_element(ResultMinTotalTFT.begin(), ResultMinTotalTFT.end()) - ResultMinTotalTFT.begin();
	if (ResultMinTotalTFT[BestFac] < TotalTFT)
	{
		if (BestFac == MoveFac)
		{
			TotalTFT -= FacTFT[BestFac];
			FacTFT[BestFac] = TestMinTFT[BestFac];
			GetEijAfterInsertAJob(Sol[BestFac], ExtJob, BestInsertPos[BestFac], NewEij, TestMinTFT[BestFac]);
			perEij[BestFac] = NewEij;
			TotalTFT += FacTFT[BestFac];
		}
		else
		{
			TotalTFT = TotalTFT - FacTFT[MoveFac] - FacTFT[BestFac];
			FacTFT[MoveFac] = NewTFT_ExtFac;
			perEij[MoveFac] = NewEij;
			FacTFT[BestFac] = TestMinTFT[BestFac];
			GetEijAfterInsertAJob(Sol[BestFac], ExtJob, BestInsertPos[BestFac], perEij[BestFac], TestMinTFT[BestFac]);
			TotalTFT = TotalTFT + FacTFT[MoveFac] + FacTFT[BestFac];
		}
		Sol[BestFac].insert(Sol[BestFac].begin() + BestInsertPos[BestFac], ExtJob);
		return true;
	}
	else
	{
		Sol[MoveFac].insert(Sol[MoveFac].begin() + MovePos, ExtJob);
	}
	return false;
}

bool Operation::CarryBestPermutationByMinGlobalTFT_randfactory(int plimit, int MoveFac, int MovePos,
	vector<vector<int>>& Sol, vector<int>& FacTFT, int& TotalTFT, vector<vector<vector<int>>>& perEij)
{
	vector<int> ref_factory;
	GetOrderSeq(m_Factories, ref_factory);
	random_shuffle(ref_factory.begin(), ref_factory.end());

	int MoveJob = Sol[MoveFac][MovePos];
	vector<vector<int>> PreEij_movefac(perEij[MoveFac]);
	int NewTFT_MoveFac = EraseAJob(Sol[MoveFac], MovePos, perEij[MoveFac]);

	bool flag = false;
	for (int tf = 0; tf < m_Factories; tf++)
	{
		int f = ref_factory[tf];
		if (f == MoveFac)
			continue;
		int Nsize = Sol[f].size();

		vector<int> ref_seq;
		GetOrderSeq(Nsize, ref_seq);
		random_shuffle(ref_seq.begin(), ref_seq.end());
		int limitjobs = min(Nsize, plimit);
		for (int j = 0; j < limitjobs; j++)
		{
			int pos = ref_seq[j];
			int job = Sol[f][pos];
			vector<vector<int>> PreEij_Seq(perEij[f]);
			EraseAJob(Sol[f], pos, perEij[f]);

			int bpos_inmovefac = -1, bpos_inseq = -1;
			int MinTFT_MoveFac = GetMinTFTAfterInsertSingleFac(Sol[MoveFac], job, bpos_inmovefac, perEij[MoveFac]);
			int MinTFT_Seq = GetMinTFTAfterInsertSingleFac(Sol[f], MoveJob, bpos_inseq, perEij[f]);

			if (MinTFT_MoveFac + MinTFT_Seq < FacTFT[MoveFac] + FacTFT[f])
			{
				TotalTFT = TotalTFT - FacTFT[MoveFac] - FacTFT[f];
				InsertAJobToSeq(MoveJob, bpos_inseq, Sol[f], perEij[f], FacTFT[f]);
				InsertAJobToSeq(job, bpos_inmovefac, Sol[MoveFac], perEij[MoveFac], FacTFT[MoveFac]);
				TotalTFT = TotalTFT + FacTFT[MoveFac] + FacTFT[f];

				flag = true;
				break;
			}
			else
			{
				Sol[f].insert(Sol[f].begin() + pos, job);
				perEij[f] = PreEij_Seq;
			}

		}
		if (flag)
			break;
	}
	if (!flag)
	{
		Sol[MoveFac].insert(Sol[MoveFac].begin() + MovePos, MoveJob);
		perEij[MoveFac] = PreEij_movefac;
	}
	return flag;
}

