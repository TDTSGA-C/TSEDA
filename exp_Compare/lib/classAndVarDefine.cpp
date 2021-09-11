#include "classAndVarDefine.h"


vector<Task> Tasks;
vector<vector<int>> TskLstInLvl;
vector<vector<double>> ParChildTranFileSizeSum;
vector<int> LevelIdOfTask;
vector<Resource> Rscs;
vector<chromosome> population;
vector<chromosome> NewPopulation;
vector<vector<chromosome>> populations;
Paramet_CGA Parameter_CGA;
Paramet_HGA Parameter_HGA;
Paramet_LWSGA Parameter_LWSGA;
Paramet_TSEDA Parameter_TSEDA;
ComConst comConst;
double ModelScale;
vector<set<int> > Descendants;
vector<set<int> > Ancestors;
