#include "Base.h"
#include "stdafx.h"

void Base::ReadInstanceFileNameList(string Dir, vector<string>& m_InstanceFileNameList)
{
	m_InstanceFileNameList.clear();
	vector<int> Known_Cmax;
	ifstream ifile;
	ifile.open(Dir + "\\" + "InstanceFileNameList.txt");
	while (true)
	{
		string FName;
		int CMax;
		ifile >> FName >> CMax;
		if (ifile.peek() != EOF)
		{
			m_InstanceFileNameList.push_back(Dir + "\\" + FName);
			Known_Cmax.push_back(CMax);
		}
		else
			break;
	}
	ifile.close();
}

void Base::ReadInstance(int InsNo, vector<string> m_InstanceFileNameList, int& Jobs, int& Machines, int& Factories, vector<vector<int>>& pTime)
{
	ifstream ifile;
	ifile.open(m_InstanceFileNameList[InsNo]);
	ifile >> Jobs >> Machines >> Factories;
	pTime.clear();
	pTime.resize(Jobs, vector<int>(Machines));
	int f;
	for (int j = 0; j < Jobs; j++)
		for (int m = 0; m < Machines; m++)
		{
			ifile >> f;
			ifile >> pTime[j][m];
		}
	ifile.close();
}

long Base::GetElapsedProcessTime()
{
	FILETIME createTime;
	FILETIME exitTime;
	FILETIME kernelTime;
	FILETIME userTime;

	if (GetProcessTimes(GetCurrentProcess(), &createTime, &exitTime, &kernelTime, &userTime) != 0)
	{
		//  Returns total user time.
		SYSTEMTIME userSystemTime;
		if (FileTimeToSystemTime(&userTime, &userSystemTime) != -1)
			return (userSystemTime.wDay - 1) * 24 * 3600 * 1000
			+ userSystemTime.wHour * 3600 * 1000 +
			userSystemTime.wMinute * 60 * 1000 +
			userSystemTime.wSecond * 1000 +
			userSystemTime.wMilliseconds;
		else
			return 0;

	}
	else
		return 0;
}
