#include "ConstructiveHeuristics.h"
#include "EA.h"
#include "IGAndILS.h"
#include "MyPBIG.h"
#include "HIG.h"
#include "Base.h"
#include "MIG.h"
#include "IGR.h"
#include "IG_2S.h"
#include "HEDFOA.h"


using namespace std;

void Con_Heu_Runmain(int SubDir, int Raps,int id_algorithm)
{
	vector<string> name_algorithm = { "0","NEHR2A4","DNEH","DLR","DLRE","DPFE","EHPF2" };
	vector<int> ChoosedIns = { 0,1,2,3,4,5,6,7,8,9 };
	cout << "THE TESTED ALGORITHMs ARE: Constructive heuristics!" << endl;
	cout << "factories:" << SubDir << endl;

	string dir = "C:\\DBFSP_Ins\\";
	int Instances = 12;
	vector<long> CPUTime;
	ofstream ofile;
	ofile.open("C:\\" + name_algorithm[id_algorithm] + ".csv");

	for (int subd = 0; subd < SubDir; subd++)
	{
		ostringstream fname;
		fname << dir;
		fname << (subd + 2);//工厂数
		vector<string> m_InstanceFileNameList;
		Base::ReadInstanceFileNameList(fname.str(), m_InstanceFileNameList);
		cout << "the numbers of factory:" << fname.str() << endl;

		long Starttime = Base::GetElapsedProcessTime();
		for (int Ins = 0; Ins < Instances; Ins++)
		{
			for (int k = 0; k < ChoosedIns.size(); k++)
			{
				vector<vector<int>> pTime;
				int factories = -1, jobs = -1, machines = -1;
				Base::ReadInstance(Ins * 10 + ChoosedIns[k], m_InstanceFileNameList, jobs, machines, factories, pTime);
				//cout << factories << ends << jobs << ends << machines << endl;
				vector<int> TFTArray(Raps);
				for (int rap = 0; rap < Raps; rap++)
				{
					srand(2021 + Instances + rap);
					ConstructiveHeuristics conheu;
					conheu.SetInstance(factories, jobs, machines, pTime);
					TFTArray[rap] = conheu.Con_Heu_Run(id_algorithm);
					//cout << TFTArray[rap] << ends;
				}
				//cout << endl;
				ofile << factories << "," << jobs << "," << machines << ",";
				for (int rap = 0; rap < Raps; rap++)
					ofile << TFTArray[rap] << ",";
				ofile << endl;
			}
		}
		long Endtime = Base::GetElapsedProcessTime();
		CPUTime.push_back(Endtime - Starttime);
	}
	ofile.close();

	ofile.open("C:\\ConHeu_CPUtime.csv", ios::app);
	ofile << name_algorithm[id_algorithm] << ",";
	long Total = 0;
	for (unsigned int i = 0; i < CPUTime.size(); i++)
	{
		ofile << CPUTime[i] << ",";
	}
	ofile << endl;
}

void EA_Runmain(int SubDir, int Raps,long time_factor)
{
	string A;
	if (time_factor == 10)
		A = "10";
	if (time_factor == 20)
		A = "20";
	if (time_factor == 40)
		A = "40";
	if (time_factor == 80)
		A = "80";

	int num_pop = 4;
	//vector<int> ChoosedIns = { 2,6 };
	vector<int> ChoosedIns = { 0,1,2,3,4,5,6,7,8,9 };
	cout << "EA" << endl;
	cout << "The number of factories:" << SubDir << endl;

	string dir = "C:\\DBFSP_Ins\\";
	int Instances = 12;
	vector<long> CPUTime;
	ofstream ofile;
	ofile.open("C:\\EA" + A + "new_5.csv");
	for (int subd = 0; subd < SubDir; subd++)
	{
		ostringstream fname;
		fname << dir;
		fname << (subd + 5);//工厂数
		vector<string> m_InstanceFileNameList;
		Base::ReadInstanceFileNameList(fname.str(), m_InstanceFileNameList);
		cout << "the numbers of factory:" << fname.str() << endl;
		for (int Ins = 0; Ins < Instances; Ins++)//Instances
		{
			for (int k = 0; k < ChoosedIns.size(); k++)//ChoosedIns.size()
			{
				srand(Raps * 2020 + Instances + 1);
				vector<vector<int>> pTime;
				int factories = -1, jobs = -1, machines = -1;
				Base::ReadInstance(Ins * 10 + ChoosedIns[k], m_InstanceFileNameList, jobs, machines, factories, pTime);
				cout << factories << ends << jobs << ends << machines << ends << k << endl;
				vector<int> TFTArray(Raps);
				for (int rap = 0; rap < Raps; rap++)
				{
					EA ea;
					ea.SetInstance(factories, jobs, machines, pTime);
					ea.SetParameters_EA(num_pop, time_factor);
					TFTArray[rap] = ea.EA_Run();
					cout << TFTArray[rap] << "," << ends;
				}
				cout << endl;
				ofile << factories << "," << jobs << "," << machines << ",";
				for (int rap = 0; rap < Raps; rap++)
					ofile << TFTArray[rap] << ",";
				ofile << endl;
			}
		}
	}
	ofile.close();
}

void IG_Runmain(int SubDir, int Raps, long time_factor)
{
	string A;
	if (time_factor == 10)
		A = "10";
	else if (time_factor == 20)
		A = "20";
	else if (time_factor == 40)
		A = "40";
	else if (time_factor == 80)
		A = "80";
	//d=9，T=0.7,2020,9,19调参
	int num_extjob = 9;
	float T = 0.7;
	int Type_initial = 0;//DLR
	vector<int> ChoosedIns = { 0,1,2,3,4,5,6,7,8,9 };
	//vector<int> ChoosedIns = { 2,6 };
	cout << "IGP10" << endl;
	cout << "The number of factories:" << SubDir << endl;

	string dir = "C:\\DBFSP_Ins\\";
	int Instances = 12;
	vector<long> CPUTime;
	ofstream ofile;
	ofile.open("C:\\IGPnew" + A + "new_5.csv");
	for (int subd = 0; subd < SubDir; subd++)
	{
		ostringstream fname;
		fname << dir;
		fname << (subd + 5);//工厂数
		vector<string> m_InstanceFileNameList;
		Base::ReadInstanceFileNameList(fname.str(), m_InstanceFileNameList);
		cout << "the numbers of factory:" << fname.str() << endl;
		for (int Ins = 0; Ins < Instances; Ins++)//Instances
		{
			for (int k = 0; k < ChoosedIns.size(); k++)//ChoosedIns.size()
			{
				srand(Raps * 2020 + Instances + 1);
				vector<vector<int>> pTime;
				int factories = -1, jobs = -1, machines = -1;
				Base::ReadInstance(Ins * 10 + ChoosedIns[k], m_InstanceFileNameList, jobs, machines, factories, pTime);
				cout << factories << ends << jobs << ends << machines << ends << k << ends;
				vector<int> TFTArray(Raps);
				for (int rap = 0; rap < Raps; rap++)
				{
					IGAndILS ig;
					ig.SetInstance(factories, jobs, machines, pTime);
					ig.SetParameters_IG(time_factor, num_extjob, T, Type_initial);
					TFTArray[rap] = ig.IG_Run();
					cout << TFTArray[rap] << "," << endl;
				}
				cout << endl;
				ofile << factories << "," << jobs << "," << machines << ",";
				for (int rap = 0; rap < Raps; rap++)
					ofile << TFTArray[rap] << ",";
				ofile << endl;
			}
		}
	}
	ofile.close();
}

void ILS_Runmain(int SubDir, int Raps,long time_factor)
{
	string A;
	if (time_factor == 10)
		A = "10";
	else if (time_factor == 20)
		A = "20";
	else if (time_factor == 40)
		A = "40";
	else if (time_factor == 80)
		A = "80";
	int num_pi = 10;
	int num_jobs = 3;
	float T = 1.1;
	int Type_Initial = 2;

	vector<int> ChoosedIns = { 0,1,2,3,4,5,6,7,8,9 };
	//vector<int> ChoosedIns = { 2,6 };
	cout << "ILS" << endl;
	cout << "The number of factories:" << SubDir << endl;

	string dir = "C:\\DBFSP_Ins\\";
	int Instances = 12;
	vector<long> CPUTime;
	ofstream ofile;
	ofile.open("C:\\ILSnew" + A + "_7.csv");
	for (int subd = 0; subd < SubDir; subd++)
	{
		ostringstream fname;
		fname << dir;
		fname << (subd + 7);//工厂数
		vector<string> m_InstanceFileNameList;
		Base::ReadInstanceFileNameList(fname.str(), m_InstanceFileNameList);
		cout << "the numbers of factory:" << fname.str() << endl;
		for (int Ins = 0; Ins < Instances; Ins++)//Instances
		{
			for (int k = 0; k < ChoosedIns.size(); k++)//ChoosedIns.size()
			{
				srand(Raps * 2020 + Instances + 1);
				vector<vector<int>> pTime;
				int factories = -1, jobs = -1, machines = -1;
				Base::ReadInstance(Ins * 10 + ChoosedIns[k], m_InstanceFileNameList, jobs, machines, factories, pTime);
				cout << factories << ends << jobs << ends << machines << ends << k << ends;
				vector<int> TFTArray(Raps);
				for (int rap = 0; rap < Raps; rap++)
				{
					IGAndILS ils;
					ils.SetInstance(factories, jobs, machines, pTime);
					ils.SetParametersSecond_ILS(time_factor, num_pi, num_jobs, T, Type_Initial);
					TFTArray[rap] = ils.ILS_Run();
					cout << TFTArray[rap] << "," << endl;
				}
				cout << endl;
				ofile << factories << "," << jobs << "," << machines << ",";
				for (int rap = 0; rap < Raps; rap++)
					ofile << TFTArray[rap] << ",";
				ofile << endl;
			}
		}
	}
	ofile.close();
}

void MyPBIG_Runmain(int SubDir, int Raps,long time_factor)
{
	string A;
	if (time_factor == 10)
		A = "10";
	else if (time_factor == 20)
		A = "20";
	else if (time_factor == 40)
		A = "40";
	else if (time_factor == 80)
		A = "80";

	int num_Pop = 10;
	int num_extjobs = 14;
	float lamda = 0.5;

	vector<int> ChoosedIns = { 0,1,2,3,4,5,6,7,8,9 };
	cout << "MyPBIG" << endl;
	cout << "The number of factories:" << SubDir << endl;

	string dir = "C:\\DBFSP_Ins\\";
	int Instances = 12;
	vector<long> CPUTime;
	ofstream ofile;
	ofile.open("C:\\PBIG" + A + "2.csv");
	for (int subd = 0; subd < SubDir; subd++)
	{
		ostringstream fname;
		fname << dir;
		fname << (subd + 2);
		vector<string> m_InstanceFileNameList;
		Base::ReadInstanceFileNameList(fname.str(), m_InstanceFileNameList);
		cout << "the numbers of factory:" << fname.str() << endl;
		for (int Ins = 0; Ins < Instances; Ins++)
		{
			for (int k = 0; k < ChoosedIns.size(); k++)
			{
				vector<vector<int>> pTime;
				int factories = -1, jobs = -1, machines = -1;
				Base::ReadInstance(Ins * 10 + ChoosedIns[k], m_InstanceFileNameList, jobs, machines, factories, pTime);
				cout << factories << ends << jobs << ends << machines << ends << k << ends;
				vector<int> TFTArray(Raps);
				for (int rap = 0; rap < Raps; rap++)
				{
					srand(Instances + 2021 + rap);
					MyPBIG pbig;
					pbig.SetInstance(factories, jobs, machines, pTime);
					pbig.SetParameters2_PBIG(num_Pop, time_factor, num_extjobs, lamda);
					TFTArray[rap] = pbig.MyPBIG_Run();
				}
				cout << endl;
				ofile << factories << "," << jobs << "," << machines << ",";
				for (int rap = 0; rap < Raps; rap++)
					ofile << TFTArray[rap] << ",";
				ofile << endl;
			}
		}
	}
	ofile.close();
}

void HIG_Runmain(int SubDir, int Raps,int time_factor)
{
	string A;
	if (time_factor == 10)
		A = "10";
	else if (time_factor == 20)
		A = "20";
	else if (time_factor == 40)
		A = "40";
	else if (time_factor == 80)
		A = "80";

	float ValueT = 4;
	int Iiter = 5000;
	float cooling = 0.9;
	float TLTlow = 0.05;
	float TLThigh = 0.1;
	int Desmax = 3;
	int Desmin = 7;

	vector<int> ChoosedIns = { 0,1,2,3,4,5,6,7,8,9 };
	//vector<int> ChoosedIns = { 2,6 };
	cout << "HIG" << endl;
	cout << "The number of factories:" << SubDir << endl;

	string dir = "C:\\DBFSP_Ins\\";
	int Instances = 12;
	vector<long> CPUTime;
	ofstream ofile;
	ofile.open("C:\\HIG" + A + "_6.csv");
	if (!ofile)
	{
		cout << "打开写入文件失败！" << endl;
		exit(666);
	}
	for (int subd = 0; subd < SubDir; subd++)
	{
		ostringstream fname;
		fname << dir;
		fname << (subd + 6);//工厂数
		vector<string> m_InstanceFileNameList;
		Base::ReadInstanceFileNameList(fname.str(), m_InstanceFileNameList);
		cout << "the numbers of factory:" << fname.str() << endl;
		for (int Ins = 0; Ins < Instances; Ins++)//Instances
		{
			for (int k = 0; k < ChoosedIns.size(); k++)//ChoosedIns.size()
			{
				srand(Raps * 2020 + Instances + 1);
				vector<vector<int>> pTime;
				int factories = -1, jobs = -1, machines = -1;
				Base::ReadInstance(Ins * 10 + ChoosedIns[k], m_InstanceFileNameList, jobs, machines, factories, pTime);
				cout << factories << ends << jobs << ends << machines << ends;
				vector<int> TFTArray(Raps);
				for (int rap = 0; rap < Raps; rap++)
				{
					HIG hig;
					hig.SetInstance(factories, jobs, machines, pTime);
					hig.SetParameters(ValueT, Iiter, cooling, TLTlow, TLThigh, Desmin, Desmax, time_factor);
					TFTArray[rap] = hig.HIG_Run();
					cout << TFTArray[rap] << "," << endl;
				}
				cout << endl;
				ofile << factories << "," << jobs << "," << machines << ",";
				for (int rap = 0; rap < Raps; rap++)
					ofile << TFTArray[rap] << ",";
				ofile << endl;
			}
		}
	}
	ofile.close();
}

void MIG_Runmain(int SubDir, int Raps, long limittime_factor)
{
	string A;
	if (limittime_factor == 10)
		A = "10";
	else if (limittime_factor == 20)
		A = "20";
	else if (limittime_factor == 40)
		A = "40";
	else if (limittime_factor == 80)
		A = "80";

	//int dmax = 7;
	//int  dmin = 2;
	//int Iter = 1000;
	//float lamda = 0.9;
	//float T0 = 5;

	/*int dmax = 13;
	int  dmin = 9;
	int Iter = 5000;
	float lamda = 0.98;
	float T0 = 18;*/

	int dmax = 13;
	int  dmin = 9;
	int Iter = 5000;
	float lamda = 0.98;
	float T0 = 20;

	vector<int> ChoosedIns = { 0,1,2,3,4,5,6,7,8,9 };
	//vector<int> ChoosedIns = { 2,6 };
	cout << "MIG" << endl;
	cout << "The number of factories:" << SubDir << endl;

	string dir = "C:\\DBFSP_Ins\\";
	int Instances = 12;
	ofstream ofile;
	ofile.open("C:\\MIG" + A + "_6.csv");
	for (int subd = 0; subd < SubDir; subd++)
	{
		ostringstream fname;
		fname << dir;
		fname << (subd + 6);//工厂数
		vector<string> m_InstanceFileNameList;
		Base::ReadInstanceFileNameList(fname.str(), m_InstanceFileNameList);
		cout << "the numbers of factory:" << fname.str() << endl;
		for (int Ins = 0; Ins < Instances; Ins++)//Instances
		{
			for (int k = 0; k < ChoosedIns.size(); k++)//ChoosedIns.size()
			{
				srand(Raps * 2020 + Instances + 1);
				vector<vector<int>> pTime;
				int factories = -1, jobs = -1, machines = -1;
				Base::ReadInstance(Ins * 10 + ChoosedIns[k], m_InstanceFileNameList, jobs, machines, factories, pTime);
				cout << factories << ends << jobs << ends << machines << ends << k << ends;
				vector<int> TFTArray(Raps);
				for (int rap = 0; rap < Raps; rap++)
				{
					MIG mig;
					mig.SetInstance(factories, jobs, machines, pTime);
					mig.SetParameters(T0, lamda, Iter, dmin, dmax, limittime_factor);
					TFTArray[rap] = mig.MIG_Run();
					cout << TFTArray[rap] << "," << endl;
				}
				cout << endl;
				ofile << factories << "," << jobs << "," << machines << ",";
				for (int rap = 0; rap < Raps; rap++)
					ofile << TFTArray[rap] << ",";
				ofile << endl;
			}
		}
	}
	ofile.close();
}

void IGR_Runmain(int SubDir, int Raps)
{
	//the last adopted value
	int d = 13;
	int nlimp = 10;
	float T = 0.9f;

	float lamda = 0.35f;
	float mu = 0.9f;

	int time_factor = 10;

	vector<int> ChoosedIns = { 0,1,2,3,4,5,6,7,8,9 };
	cout << "IGR" << endl;
	cout << "The number of factories:" << SubDir << endl;

	string dir = "C:\\DBFSP_Ins\\";
	int Instances = 12;
	ofstream ofile;
	ofile.open("C:\\IGRnew10_2.csv");
	for (int subd = 0; subd < SubDir; subd++)
	{
		ostringstream fname;
		fname << dir;
		fname << (subd + 2);//工厂数
		vector<string> m_InstanceFileNameList;
		Base::ReadInstanceFileNameList(fname.str(), m_InstanceFileNameList);
		cout << "the numbers of factory:" << fname.str() << endl;
		for (int Ins = 0; Ins < Instances; Ins++)//Instances
		{
			for (int k = 0; k < ChoosedIns.size(); k++)//ChoosedIns.size()
			{
				srand(Raps * 2020 + Instances + 1);
				vector<vector<int>> pTime;
				int factories = -1, jobs = -1, machines = -1;
				Base::ReadInstance(Ins * 10 + ChoosedIns[k], m_InstanceFileNameList, jobs, machines, factories, pTime);
				cout << factories << ends << jobs << ends << machines << ends;
				vector<int> TFTArray(Raps);
				for (int rap = 0; rap < Raps; rap++)
				{
					IGR igr;
					igr.SetInstance(factories, jobs, machines, pTime);
					igr.SetParameter(d, nlimp, T, time_factor);
					TFTArray[rap] = igr.IGR_Run(lamda, mu);
					cout << TFTArray[rap] << "," << endl;
				}
				cout << endl;
				ofile << factories << "," << jobs << "," << machines << ",";
				for (int rap = 0; rap < Raps; rap++)
					ofile << TFTArray[rap] << ",";
				ofile << endl;
			}
		}
	}
	ofile.close();
}

void IGRTune_Runmain(int SubDir, int Raps) 
{
	//the first calibration
	/*vector<int> njobs = { 5,10,15 };
	vector<int> nlimp = { 20,40,60 };
	vector<float> T_factor = { 0.3,0.6,0.9 };*/

	//the second calibration
	/*vector<int> njobs = { 11,13,15 };
	vector<int> nlimp = { 10,20,30 };
	vector<float> T_factor = { 0.7,0.9,1.1 };*/

	//supply the second calibration
	/*vector<int> njobs = { 11,13,15 };
	vector<int> nlimp = { 5};
	vector<float> T_factor = { 0.7,0.9,1.1 };*/

	//supply the last instance 
	vector<int> njobs = { 11,13,15 };
	vector<int> nlimp = { 5,10,20 };
	vector<float> T_factor = { 0.7,0.9,1.1 };

	int time_factor = 10;

	cout << "IGRTune" << endl;
	cout << "The number of factories:" << SubDir << endl;
	string dir = "C:\\DBFSP_CaliIns\\";
	int Instances = 2;
	vector<long> CPUTime;
	ofstream ofile;

	ofile.open("C:\\IGRTuninglastsupply_7.csv");
	if (ofile)
		cout << "打开写入文件成功";
	else
		cout << "打开写入文件失败" << endl;
	for (int subd = 0; subd < SubDir; subd++)
	{
		ostringstream fname;
		fname << dir;
		fname << (subd + 7);//工厂数
		vector<string> m_InstanceFileNameList;
		Base::ReadInstanceFileNameList(fname.str(), m_InstanceFileNameList);
		cout << "the numbers of factory:" << fname.str() << endl;
		for (int Ins = 22; Ins < 22 + Instances; Ins++)
		{
			srand(2020 + Raps + SubDir);
			vector<vector<int>> pTime;
			int factories = -1, jobs = -1, machines = -1;
			Base::ReadInstance(Ins, m_InstanceFileNameList, jobs, machines, factories, pTime);
			cout << factories << ends << jobs << ends << machines << endl;
			for (int nj = 0; nj < njobs.size(); nj++)
			{
				for (int nl = 0; nl < nlimp.size(); nl++)
				{
					for (int tf = 0; tf < T_factor.size(); tf++)
					{
						cout << njobs[nj] << ends << nlimp[nl] << ends << T_factor[tf] << endl;
						vector<int> TFTArray(Raps);
						for (int rap = 0; rap < Raps; rap++)
						{
							IGR igr;
							igr.SetInstance(factories, jobs, machines, pTime);
							igr.SetParameter(njobs[nj], nlimp[nl], T_factor[tf], time_factor);
							TFTArray[rap] = igr.IGR_Run(0.35f, 0.9f);//lamda=0.35,mu=0.9
							cout << TFTArray[rap] << "," << ends;
						}
						cout << endl;
						ofile << factories << "," << jobs << "," << machines << ",";
						ofile << njobs[nj] << "," << nlimp[nl] << "," << T_factor[tf] << ",";
						for (int rap = 0; rap < Raps; rap++)
							ofile << TFTArray[rap] << ",";
						ofile << endl;
					}
				}
			}
		}
	}
	ofile.close();
}

void IG2S_Runmain(int SubDir, int Raps, long time_factor)
{
	string A;
	if (time_factor == 10)
		A = "10";
	else if (time_factor == 20)
		A = "20";
	else if (time_factor == 40)
		A = "40";
	else if (time_factor == 80)
		A = "80";

	int sd1 = 12;
	int sd2 = 9;
	double sT = 1.0;
	double sratio = 0.97;

	vector<int> ChoosedIns = { 0,1,2,3,4,5,6,7,8,9 };
	//vector<int> ChoosedIns = { 2,6 };
	cout << "IG2S" << endl;
	cout << "The number of factories:" << SubDir << endl;

	string dir = "C:\\DBFSP_Ins\\";
	int Instances = 12;
	vector<long> CPUTime;
	ofstream ofile;
	ofile.open("C:\\IG2S" + A + "_6.csv");
	for (int subd = 0; subd < SubDir; subd++)
	{
		ostringstream fname;
		fname << dir;
		fname << (subd + 6);//工厂数
		vector<string> m_InstanceFileNameList;
		Base::ReadInstanceFileNameList(fname.str(), m_InstanceFileNameList);
		cout << "the numbers of factory:" << fname.str() << endl;
		for (int Ins = 0; Ins < Instances; Ins++)//Instances
		{
			for (int k = 0; k < ChoosedIns.size(); k++)//ChoosedIns.size()
			{
				vector<vector<int>> pTime;
				int factories = -1, jobs = -1, machines = -1;
				Base::ReadInstance(Ins * 10 + ChoosedIns[k], m_InstanceFileNameList, jobs, machines, factories, pTime);
				cout << factories << ends << jobs << ends << machines << ends << k << ends;
				vector<int> TFTArray(Raps);
				for (int rap = 0; rap < Raps; rap++)
				{
					srand(Instances + 2021 + rap);
					IG_2S ig2s;
					ig2s.SetInstance(factories, jobs, machines, pTime);
					ig2s.SetParameter_IG2S(sd1, sd2, sT, sratio, time_factor);
					TFTArray[rap] = ig2s.IG2S_Run();
					cout << TFTArray[rap] << "," << endl;
				}
				cout << endl;
				ofile << factories << "," << jobs << "," << machines << ",";
				for (int rap = 0; rap < Raps; rap++)
					ofile << TFTArray[rap] << ",";
				ofile << endl;
			}
		}
	}
	ofile.close();
}

void HEDFOA_Runmain(int SubDir, int Raps, long time_factor)
{
	string A;
	if (time_factor == 10)
		A = "10";
	else if (time_factor == 20)
		A = "20";
	else if (time_factor == 40)
		A = "40";
	else if (time_factor == 80)
		A = "80";

	int SNP = 1;
	int SSN = 3;
	int Sd = 10;
	float ST = 1.0f;

	vector<int> ChoosedIns = { 0,1,2,3,4,5,6,7,8,9 };
	//vector<int> ChoosedIns = { 2,6 };
	cout << "HEDFOA" << endl;
	cout << "The number of factories:" << SubDir << endl;

	string dir = "C:\\DBFSP_Ins\\";
	int Instances = 12;
	vector<long> CPUTime;
	ofstream ofile;
	ofile.open("C:\\HEDFOA" + A + "_2.csv");
	for (int subd = 0; subd < SubDir; subd++)
	{
		ostringstream fname;
		fname << dir;
		fname << (subd + 2);//工厂数
		vector<string> m_InstanceFileNameList;
		Base::ReadInstanceFileNameList(fname.str(), m_InstanceFileNameList);
		cout << "the numbers of factory:" << fname.str() << endl;
		for (int Ins = 0; Ins < Instances; Ins++)//Instances
		{
			for (int k = 0; k < ChoosedIns.size(); k++)//ChoosedIns.size()
			{
				vector<vector<int>> pTime;
				int factories = -1, jobs = -1, machines = -1;
				Base::ReadInstance(Ins * 10 + ChoosedIns[k], m_InstanceFileNameList, jobs, machines, factories, pTime);
				cout << factories << ends << jobs << ends << machines << ends << k << endl;
				vector<int> TFTArray(Raps);
				for (int rap = 0; rap < Raps; rap++)
				{
					srand(2021 + Instances + rap);
					HEDFOA hedfoa;
					hedfoa.SetInstance(factories, jobs, machines, pTime);
					hedfoa.SetParameter_HEDFOA(SNP, SSN, Sd, ST, time_factor);
					TFTArray[rap] = hedfoa.HEDFOA_Run();
				}
				ofile << factories << "," << jobs << "," << machines << ",";
				for (int rap = 0; rap < Raps; rap++)
					ofile << TFTArray[rap] << ",";
				ofile << endl;
			}
		}
	}
	ofile.close();
}

void IG2S_TuningRunmain(int SubDir, int Raps)
{
	long time_factor = 10;
	//the first calibration 
	vector<float> S_ValueT = { 0.4f,0.8f,1.2f,1.6f };
	vector<int> S_d1 = { 4,8,12,16 };
	vector<int> S_d2 = { 3,6,9,12 };
	vector<float> S_Ratio = { 0.91f,0.94f,0.97f,1.0f };

	vector<vector<int>> L16 = {
	{0,0,0,0},{0,1,1,1},{0,2,2,2},{0,3,3,3},
	{1,0,1,2},{1,1,0,3},{1,2,3,0},{1,3,2,1},
	{2,0,2,3},{2,1,3,2},{2,2,0,1},{2,3,1,0},
	{3,0,3,1},{3,1,2,0},{3,2,1,3},{3,3,0,2}
	};

	cout << "THE TESTED ALGORITHM IS: IG2S" << endl;
	cout << "factories:" << SubDir << endl;
	string dir = "C:\\DBFSP_CaliIns\\";
	int Instances = 24;
	vector<long> CPUTime;
	ofstream ofile;

	ofile.open("C:\\IG2STuning_7.csv");
	if (ofile)
		cout << "打开写入文件成功";
	else
		cout << "打开写入文件失败" << endl;
	for (int subd = 0; subd < SubDir; subd++)
	{
		ostringstream fname;
		fname << dir;
		fname << (subd + 7);//工厂数
		vector<string> m_InstanceFileNameList;
		Base::ReadInstanceFileNameList(fname.str(), m_InstanceFileNameList);
		cout << "the numbers of factory:" << fname.str() << endl;
		for (int Ins = 0; Ins < Instances; Ins++)
		{
			srand(2020 + Raps + SubDir);
			vector<vector<int>> pTime;
			int factories = -1, jobs = -1, machines = -1;
			Base::ReadInstance(Ins, m_InstanceFileNameList, jobs, machines, factories, pTime);
			cout << factories << ends << jobs << ends << machines << endl;
			for (int i = 0; i < L16.size(); i++)
			{
				float _valueT = S_ValueT[L16[i][0]];
				int _d1 = S_d1[L16[i][1]];
				int _d2 = S_d2[L16[i][2]];
				float _Ratio = S_Ratio[L16[i][3]];
				cout << _valueT << ends << _d1 << ends << _d2 << ends << _Ratio << endl;
				vector<int> TFTArray(Raps);
				for (int rap = 0; rap < Raps; rap++)
				{
					IG_2S ig2s;
					ig2s.SetInstance(factories, jobs, machines, pTime);
					ig2s.SetParameter_IG2S(_d1, _d2, _valueT, _Ratio, time_factor);
					TFTArray[rap] = ig2s.IG2S_Run();
					cout << TFTArray[rap] << ends << "," << ends;
				}
				cout << endl;
				ofile << factories << "," << jobs << "," << machines << ",";
				ofile << _valueT << "," << _d1 << "," << _d2 << "," << _Ratio << ",";
				for (int rap = 0; rap < Raps; rap++)
					ofile << TFTArray[rap] << ",";
				ofile << endl;
			}
		}
	}
	ofile.close();
}

void IG2S_TuningRunmain2(int SubDir, int Raps)
{
	long time_factor = 10;
	//the results of the first calibration
	int _d1 = 12;
	int _d2 = 9;
	vector<float> S_ValueT = { 0.8f,1.0f,1.2f,1.4f };
	vector<float> S_Ratio = { 0.94f,0.97f,1.0f };

	cout << "THE TESTED ALGORITHM IS: IG2S" << endl;
	cout << "factories:" << SubDir << endl;
	string dir = "C:\\DBFSP_CaliIns\\";
	int Instances = 24;
	vector<long> CPUTime;
	ofstream ofile;

	ofile.open("C:\\IG2STuning2_7.csv");
	if (ofile)
		cout << "打开写入文件成功";
	else
		cout << "打开写入文件失败" << endl;
	for (int subd = 0; subd < SubDir; subd++)
	{
		ostringstream fname;
		fname << dir;
		fname << (subd + 7);//工厂数
		vector<string> m_InstanceFileNameList;
		Base::ReadInstanceFileNameList(fname.str(), m_InstanceFileNameList);
		cout << "the numbers of factory:" << fname.str() << endl;
		for (int Ins = 0; Ins < Instances; Ins++)
		{
			srand(2020 + Raps + SubDir);
			vector<vector<int>> pTime;
			int factories = -1, jobs = -1, machines = -1;
			Base::ReadInstance(Ins, m_InstanceFileNameList, jobs, machines, factories, pTime);
			cout << factories << ends << jobs << ends << machines << endl;
			for (int tt = 0; tt < S_ValueT.size(); tt++)
			{
				for (int rr = 0; rr < S_Ratio.size(); rr++)
				{
					float _valueT = S_ValueT[tt];
					float _Ratio = S_Ratio[rr];
					cout << _valueT << ends << _Ratio << endl;
					vector<int> TFTArray(Raps);
					for (int rap = 0; rap < Raps; rap++)
					{
						IG_2S ig2s;
						ig2s.SetInstance(factories, jobs, machines, pTime);
						ig2s.SetParameter_IG2S(_d1, _d2, _valueT, _Ratio, time_factor);
						TFTArray[rap] = ig2s.IG2S_Run();
						cout << TFTArray[rap] << ends << "," << ends;
					}
					cout << endl;
					ofile << factories << "," << jobs << "," << machines << ",";
					ofile << _valueT << "," << _d1 << "," << _d2 << "," << _Ratio << ",";
					for (int rap = 0; rap < Raps; rap++)
						ofile << TFTArray[rap] << ",";
					ofile << endl;
				}
			}
		}
	}
	ofile.close();
}

void HEDFOA_TuningRunmain(int SubDir, int Raps)
{
	long time_factor = 10;
	//the first calibration 

	vector<int> S_NP = { 1,2,3,4 };
	vector<int> S_SN = { 1,2,3,4 };
	vector<int> S_d = { 6, 8,10,12 };
	vector<float> S_T = { 0.6f,0.8f,1.0f,1.2f };

	/*vector<vector<int>> L16;
	for (int np = 0; np < S_NP.size(); np++)
	{
		for (int sn = 0; sn < S_SN.size(); sn++)
		{
			for (int d = 0; d < S_d.size(); d++)
			{
				for (int t = 0; t < S_T.size(); t++)
				{
					vector<int> Test(4);
					Test[0] = np; Test[1] = sn; Test[2] = d; Test[3] = t;
					L16.push_back(Test);
				}
			}
		}
	}*/

	vector<vector<int>> L16 = {
	{0,0,0,0},{0,1,1,1},{0,2,2,2},{0,3,3,3},
	{1,0,1,2},{1,1,0,3},{1,2,3,0},{1,3,2,1},
	{2,0,2,3},{2,1,3,2},{2,2,0,1},{2,3,1,0},
	{3,0,3,1},{3,1,2,0},{3,2,1,3},{3,3,0,2}
	};

	cout << "THE TESTED ALGORITHM IS: HEDFOA" << endl;
	cout << "factories:" << SubDir << endl;
	string dir = "C:\\DBFSP_CaliIns\\";
	int Instances = 24;
	vector<long> CPUTime;
	ofstream ofile;

	ofile.open("C:\\HEDFOATuning6_2.csv");
	if (ofile)
		cout << "打开写入文件成功";
	else
		cout << "打开写入文件失败" << endl;
	for (int subd = 0; subd < SubDir; subd++)
	{
		ostringstream fname;
		fname << dir;
		fname << (subd + 2);//工厂数
		vector<string> m_InstanceFileNameList;
		Base::ReadInstanceFileNameList(fname.str(), m_InstanceFileNameList);
		cout << "the numbers of factory:" << fname.str() << endl;
		for (int Ins = 20; Ins < Instances; Ins++)
		{
			srand(2020 + Raps + SubDir);
			vector<vector<int>> pTime;
			int factories = -1, jobs = -1, machines = -1;
			Base::ReadInstance(Ins, m_InstanceFileNameList, jobs, machines, factories, pTime);
			cout << factories << ends << jobs << ends << machines << endl;
			for (int i = 0; i < L16.size(); i++)
			{
				int NP = S_NP[L16[i][0]];
				int SN = S_SN[L16[i][1]];
				int d = S_d[L16[i][2]];
				float T = S_T[L16[i][3]];

				cout << NP<< ends << SN << ends << d << ends << T << endl;
				vector<int> TFTArray(Raps);
				for (int rap = 0; rap < Raps; rap++)
				{
					HEDFOA hedfoa;
					hedfoa.SetInstance(factories, jobs, machines, pTime);
					hedfoa.SetParameter_HEDFOA(NP, SN, d, T, time_factor);
					TFTArray[rap] = hedfoa.HEDFOA_Run();
					cout << TFTArray[rap] << ends << "," << ends;
				}
				cout << endl;
				ofile << factories << "," << jobs << "," << machines << ",";
				ofile << NP << "," << SN << "," << d << "," << T << ",";
				for (int rap = 0; rap < Raps; rap++)
					ofile << TFTArray[rap] << ",";
				ofile << endl;
			}
		}
	}
	ofile.close();
}

int main()
{
	//MyPBIG_Runmain(3, 5, 20);
}

