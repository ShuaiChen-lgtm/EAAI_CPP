#include "Problem.h"

Problem::Problem()
{
	m_BestSol_Seq.resize(0);
	m_BestSol_FacTFT.resize(0);
	m_BestSol_TotalTFT = INT_MAX;
}

Problem::~Problem()
{
}

void Problem::SetInstance(int factories, int jobs, int machines, const vector<vector<int>>& pTime)
{
	this->m_Jobs = jobs;
	this->m_Factories = factories;
	this->m_Machines = machines;

	this->m_Ptime.resize(this->m_Jobs, vector<int>(this->m_Machines));
	for (int j = 0; j < this->m_Jobs; j++)
		for (int m = 0; m < this->m_Machines; m++)
			this->m_Ptime[j][m] = pTime[j][m];
}

void Problem::GetTotalPTime()
{
	this->m_TotalPTime.clear();
	this->m_TotalPTime.resize(this->m_Jobs, 0);
	for (int j = 0; j < this->m_Jobs; j++)
	{
		this->m_TotalPTime[j] = this->m_Ptime[j][0];
		for (int m = 1; m < this->m_Machines; m++)
			this->m_TotalPTime[j] += this->m_Ptime[j][m];
	}
}

void Problem::GetTotalLoad()
{
	this->m_TotalLoad = 0;
	for (int j = 0; j < m_Jobs; j++)
	{
		this->m_TotalLoad += this->m_TotalPTime[j];
	}
}

