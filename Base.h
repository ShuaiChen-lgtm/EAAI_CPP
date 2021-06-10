//Base.h
#pragma once
#include <vector>
#include <string>

using namespace std;

namespace Base
{
	template <typename T>
	void Free(vector<T>& p);

	template <typename T>
	void Free(vector<vector<T>>& p);

	template <typename T>
	void Free(vector<vector<vector<T>>>& p);

	template <typename T>
	struct Pair
	{
		int dim;
		T value;
	};

	template <typename T>
	struct PairGreater {
		bool operator () (Pair <T> a, Pair <T> b)
		{
			return a.value > b.value;
		}
	};

	template <typename T>
	struct PairLess {
		bool operator () (Pair <T> a, Pair <T> b)
		{
			return a.value < b.value;
		}
	};

	void ReadInstanceFileNameList(string Dir, vector<string>& m_InstanceFileNameList);
	void ReadInstance(int InsNo, vector<string> m_InstanceFileNameList, int& Jobs, int& Machines, int& Factories, vector<vector<int>>& pTime);
	long GetElapsedProcessTime();
}

template <typename T>
void Base::Free(vector<T>& p)
{
	vector<T> pTemp;
	p.swap(pTemp);
}

template <typename T>
void Base::Free(vector<vector<T>>& p)
{
	vector<vector<T>> pTemp;
	p.swap(pTemp);
}

template <typename T>
void Base::Free(vector<vector<vector<T>>>& p)
{
	vector<vector<vector<T>>> pTemp;
	p.swap(pTemp);
}

