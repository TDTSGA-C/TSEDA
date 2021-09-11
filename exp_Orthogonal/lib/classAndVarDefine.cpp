#include "classAndVarDefine.h"


vector<Task> Tasks;
vector<vector<int>> TskLstInLvl;
vector<vector<double>> ParChildTranFileSizeSum;
vector<int> LevelIdOfTask;
vector<Resource> Rscs;
vector<chromosome> population;
vector<chromosome> NewPopulation;
vector<vector<chromosome>> populations;
Paramet_TSEDA Parameter_TSEDA;
ComConst comConst;
double ModelScale;
vector<Orthogonal> orthogonal;
vector<set<int> > Descendants;
vector<set<int> > Ancestors;