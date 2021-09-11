#ifndef CSTCHANGE_CLASSDEFINE_H

#include "math.h"
#include <vector>
#include <string>
#include <set>
#include <map>

#define CSTCHANGE_CLASSDEFINE_H
using namespace std;
//{file}
class vfile {
public:
    string FileName;     //file name
    int source;          //the source of file, -1:from the shared server; i: from task i
    double size;         //the size of file
};

//{task}
class Task {
public:
    double length;          //the length of task
    vector<int> ElgRsc;     //the set of resources which are eligible to  perform this task -xy4
    vector<int> parents;    //the set of parent tasks
    vector<int> children;   //the set of child tasks
    vector<vfile> IFile;    //the set of input files
    vector<vfile> OFile;    //the set of output files
};

//{resources}
class Resource { //-w Resource
public:
    vector<int> ElgTsk;     //the set of tasks which it is eligible to perform -xy4
    double pc, bw;          //processing capacity, bandwidth

    Resource(int id, double pc, double bw) {
        this->bw = bw;
        this->pc = pc;
    }
};

class chromosome {
public:
    vector<int> RscAlcLst;       //resources allocation, task-to-resource mapping, (Match String)
    vector<int> TskSchLst;       //task scheduling order (Scheduling String)
    vector<double> EndTime;      //the finish time of task
    double FitnessValue;         //Fitness (makespan)
    //{sort chromosome and remove the same chromosome according to fitness value}
//    bool operator<(const chromosome &otherChromosome)const {
//        return this->FitnessValue + 1e-6 < otherChromosome.FitnessValue;
//    }
//    //{sort chromosome according to fitness value and remove the same chromosome according to code }
//    bool operator<(const chromosome &otherChromosome)const {
//        int flag = -1;
//        for (int i = 0; i < this->TskSchLst.size(); i++) {
//            if (this->RscAlcLst[i] != otherChromosome.RscAlcLst[i]) {
//                flag = i;
//                break;
//            }
//            if (this->TskSchLst[i] != otherChromosome.TskSchLst[i]) {
//                flag = i;
//                break;
//            }
//        }
//        if (flag == -1) {
//            return false;
//        } else {
//            if(fabs(this->FitnessValue-otherChromosome.FitnessValue)<1e-6) {
//                if(this->RscAlcLst[flag] != otherChromosome.RscAlcLst[flag]){
//                    return this->RscAlcLst[flag]<otherChromosome.RscAlcLst[flag];
//                }else {
//                    return this->TskSchLst[flag]<otherChromosome.TskSchLst[flag];
//                }
//            }else {
//                return this->FitnessValue < otherChromosome.FitnessValue;
//            }
//        }
//    }
};

class Paramet_TSEDA {
public:
    int NumOfChormPerPop;//种群规模
    int NumOfEliteOfPop;//精英个数
    double RunTimeRatioOfStg1;//Stage1终止条件
    double alpha;//概率矩阵更新系数
    double beta;//概率矩阵更新系数
    int NumOfImproveOfPop;//改进个数
};

class ComConst {
public:
    int NumOfTsk;             //the number of Tasks
    int NumOfRsc;             //the number of resources
};

class Orthogonal{
public:
    double NumOfChormPerPopRate;//种群规模系数
    double NumOfEliteOfPopRate;//精英率
    double alpha;//概率矩阵更新系数
    double beta;//概率矩阵更新系数
    double NumOfImproveOfPopRate;//改进率
    double RunTimeRatioOfStg1Rate;//Stage1终止条件系数
};

extern vector<Task> Tasks;
extern vector<vector<int>> TskLstInLvl; //task list (set) in each level
extern vector<vector<double>> ParChildTranFileSizeSum;
extern vector<int> LevelIdOfTask;      //the level of task
extern vector<Resource> Rscs;
extern vector<chromosome> population;
extern vector<chromosome> NewPopulation;
extern vector<vector<chromosome>> populations;
extern Paramet_TSEDA Parameter_TSEDA;
extern ComConst comConst;
extern double ModelScale;
extern vector<set<int> > Descendants;
extern vector<set<int> > Ancestors;
extern vector<Orthogonal> orthogonal;
#endif //CSTCHANGE_CLASSDEFINE_H