#include <cstdlib>
#include "GenerateAChrom.h"
#include "GenOperator.h"
#include "tools.hpp"

//{calculate the average execution time of tasks}
void W_Cal_Average(vector<double>& w) {
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        double sum = 0;
        int RscSize = Tasks[i].ElgRsc.size();
        for (int j = 0; j < RscSize; ++j)
            sum += 1.0 / Rscs[Tasks[i].ElgRsc[j]].pc;
        w[i] = Tasks[i].length * sum / RscSize;
    }
}

//{calculate the average transfer time among tasks}
void C_Cal_Average(vector<vector<double>>& c) {
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        if(Tasks[i].parents.size() == 0){
            continue;
        }
        for (int j = 0; j < Tasks[i].parents.size(); ++j) {
            int parent = Tasks[i].parents[j];
            double sum1 = 0;
            double sum = 0;
            sum = ParChildTranFileSizeSum[parent][i] / VALUE;
            for (int k = 0; k < Tasks[i].ElgRsc.size(); ++k) {
                for (int y = 0; y < Tasks[parent].ElgRsc.size(); ++y) {
                    if (Tasks[i].ElgRsc[k] == Tasks[parent].ElgRsc[y]) {
                        continue;
                    } else {
                        sum1 += sum * 8 / XY_MIN(Rscs[Tasks[i].ElgRsc[k]].bw, Rscs[Tasks[parent].ElgRsc[y]].bw);
                    }
                }
            }
            c[parent][i] = sum1 / (double) (Tasks[i].ElgRsc.size() * Tasks[parent].ElgRsc.size());
        }
    }
}

//calculate the rank of tasks based on independent IO using transfer time C[i][j]
void Calculate_Rank_b(vector<double>& RankList, vector<vector<double>>& c, vector<double>& w){
    for (int i = 0; i < TskLstInLvl[TskLstInLvl.size()-1].size(); ++i) {
        int TaskId=TskLstInLvl[TskLstInLvl.size()-1][i];
        RankList[TaskId] = w[TaskId];
    }
    for(int i =TskLstInLvl.size()-2 ;i >=0 ;--i){
        for (int j = 0; j < TskLstInLvl[i].size(); ++j) {
            int TaskId=TskLstInLvl[i][j];
            double ChildMaxRankc = 0;
            for (int k = 0; k < Tasks[TaskId].children.size(); ++k) {
                int tem = Tasks[TaskId].children[k];
                double CompareObject = RankList[tem] + c[TaskId][tem];
                if(ChildMaxRankc  < CompareObject ){
                    ChildMaxRankc = CompareObject;
                }
            }
            RankList[TaskId] =w[TaskId] + ChildMaxRankc;
        }
    }
}

//{initialize chromosome to allocate spaces}
void IntChr(chromosome& chrom) {
    chrom.TskSchLst.resize(comConst.NumOfTsk);
    chrom.RscAlcLst.resize(comConst.NumOfTsk);
    chrom.EndTime.resize(comConst.NumOfTsk);
}

double HrsDcd_EFT(chromosome& ch) {
    for (int i = 0; i < comConst.NumOfTsk; ++i)
        ch.RscAlcLst[i] = -1;
    vector<set<double> > ITL;                           //the idle time-slot lists  for all resources
    double makespan = 0;
    for (int j = 0; j < comConst.NumOfRsc; ++j) {
        set<double> a;
        a.insert(0.0);
        a.insert(99999999 * 1.0);
        ITL.push_back(a);
    }
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        int RscIndex = -1;
        int TaskIndex = ch.TskSchLst[i];
        double FinalEndTime = 100000000000;
        double FinalStartTime = 0;
        for (int j = 0; j < Tasks[TaskIndex].ElgRsc.size(); ++j) {
            double ReadyTime = 0;
            int v = Tasks[TaskIndex].ElgRsc[j];
            if(Tasks[TaskIndex].parents.size() != 0){
                for (int n = 0; n < Tasks[TaskIndex].parents.size(); ++n) {
                    int ParentIndex = Tasks[TaskIndex].parents[n];
                    int ParentRscIndex = ch.RscAlcLst[ParentIndex];
                    double max = ch.EndTime[ParentIndex];
                    if(v != ParentRscIndex){
                        double TransferData = ParChildTranFileSizeSum[ParentIndex][TaskIndex];
                        max += TransferData / VALUE * 8 / (XY_MIN(Rscs[v].bw,Rscs[ParentRscIndex].bw));
                    }
                    if (ReadyTime < max){
                        ReadyTime = max;
                    }
                }
            }
            double ExeTime = Tasks[TaskIndex].length / Rscs[v].pc;
            double StartTime = 0;
            double EndTime = 0;
            //{Find an idle time-slot as early as possible from ITL}
            set<double>::iterator pre  = ITL[v].begin();
            set<double>::iterator post = ITL[v].begin();
            ++post;
            while(post != ITL[v].end()) {
                if((*post - *pre) >= ExeTime && ReadyTime <= (*post)-ExeTime) {
                    StartTime = XY_MAX(*pre, ReadyTime);
                    break;
                } else {
                    ++pre;
                    ++pre;
                    ++post;
                    ++post;
                }
            }
            EndTime = StartTime + ExeTime;
            //{find/record the earliest finish time}
            if (EndTime < FinalEndTime) {
                FinalStartTime = StartTime;
                FinalEndTime = EndTime;
                RscIndex = v;
            }
        }
        ch.EndTime[TaskIndex] = FinalEndTime;
        ch.RscAlcLst[TaskIndex] = RscIndex;
        //{update ITL}
        if(ITL[RscIndex].find(FinalStartTime) != ITL[RscIndex].end()) {
            ITL[RscIndex].erase(FinalStartTime);
        } else {
            ITL[RscIndex].insert(FinalStartTime);
        }
        if(ITL[RscIndex].find(ch.EndTime[TaskIndex]) != ITL[RscIndex].end()) {
            ITL[RscIndex].erase(ch.EndTime[TaskIndex]);
        } else {
            ITL[RscIndex].insert(ch.EndTime[TaskIndex]);
        }
        makespan = XY_MAX(makespan, FinalEndTime);
    }
    ch.FitnessValue = makespan;
    return makespan;
}

chromosome GnrChr_HEFT(vector<double> Rank_b) {
    vector<int> ind(comConst.NumOfTsk);
    chromosome TemChrom;
    IntChr(TemChrom);
    IndexSort(ind, Rank_b);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        TemChrom.TskSchLst[i] = ind[comConst.NumOfTsk - i - 1];
    }
    HrsDcd_EFT(TemChrom);
    return TemChrom;
}

void InitProModelOfResAlc(vector<vector<double> >& PMR) {
    for(int i = 0; i < comConst.NumOfTsk; ++i) {
        int numOfAvailVMs = (int)Tasks[i].ElgRsc.size();
        double prob = 1.0 / numOfAvailVMs;
        for(int j = 0; j < numOfAvailVMs; ++j) {
            PMR[i][Tasks[i].ElgRsc[j]] =  prob;
        }
    }
}

void initProModelTask(vector<vector<double> >& PMS,
                      vector<int>& NumOfAncestors,
                      vector<int>& NumOfNonAncestors,
                      vector<int>& STS) {
    for(int i = 0; i < comConst.NumOfTsk; ++i) {
        int left  = NumOfAncestors[i];
        int right = NumOfNonAncestors[i];
        for(int j = left; j <= right; ++j) {
            PMS[i][j] = 1;
            ++STS[j];
        }
    }

    for(int j = 0; j < comConst.NumOfTsk; ++j) {
        double val = 1.0 / STS[j];
        for(int i = 0; i < comConst.NumOfTsk; ++i) {
            if(PMS[i][j]) {
                PMS[i][j] = val;
            }
        }
    }
}

int ChooseTask(vector<vector<double> >& PMS, int order, vector<int>& RTI) {
    double sum = 0;
    for (int i = 0; i < RTI.size(); ++i) {
        sum += PMS[RTI[i]][order];
    }

    if(fabs(sum - 0) < PrecisionValue) {
        return (*RTI.begin());
    }

    vector<double> SltProb(RTI.size());
    for (int i = 0; i < RTI.size(); ++i) {
        SltProb[i] = PMS[RTI[i]][order] / sum;
    }

    double rd_num = double(rand()%100) / 100;
    double ProbSum = 0;
    for (int i = 0; i < RTI.size(); ++i) {
        ProbSum += SltProb[i];
        if(ProbSum > rd_num) {
            return RTI[i];
        }
    }

    return -1;
}

int ChooseRes(vector<vector<double> >& PMVM, int taskId) {
    double rd_num = double(rand()%100) / 100;
    double sum = 0;
    for(int j = 0; j < comConst.NumOfRsc; ++j) {
        sum += PMVM[taskId][j];
        if(sum > rd_num) {
            return j;
        }
    }
    return -1;
}

chromosome GnrTskLstOfChr(vector<vector<double> >& PMS) {
    chromosome chrom;
    IntChr(chrom);
    for (int i = 0; i < comConst.NumOfTsk; i++) {
        chrom.TskSchLst[i] = -1;
    }


    vector<int > upr(comConst.NumOfTsk,0.0);
    vector<int> RTI;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i]=Tasks[i].parents.size();
        if (upr[i]==0){
            RTI.push_back(i);
        }
    }

    //startDecode
    for (int i = 0; i < comConst.NumOfTsk; i++) {
        int taskIndex = ChooseTask(PMS, i,RTI);
        if (taskIndex == -1) {
            cout << "wrong at func ChooseTask of taskIndex == -1" << endl;
            exit(0);
        }
        chrom.TskSchLst[i] = taskIndex;
        RTI.erase(find(RTI.begin(), RTI.end(), taskIndex));
        for (int k = 0; k < Tasks[taskIndex].children.size(); ++k) {
            upr[Tasks[taskIndex].children[k]]--;
            if (upr[Tasks[taskIndex].children[k]] == 0){
                RTI.push_back(Tasks[taskIndex].children[k]);
            }
        }
    }
    return chrom;
}

void GnrRscLstOfChr(chromosome& chrom,
                       vector<vector<double> >& PMR) {
    for (int i = 0; i < comConst.NumOfTsk; i++) {
        chrom.RscAlcLst[i] = -1;
    }

    //startDecode
    for (int i = 0; i < comConst.NumOfTsk; i++) {
        int taskIndex = chrom.TskSchLst[i];
        //选择虚拟机
        int RscIndex = ChooseRes(PMR,taskIndex);
        if(RscIndex == -1) {
            cout<<"wrong at func ChooseVM of VMIndex == -1"<<endl;
            exit(0);
        }
        chrom.RscAlcLst[taskIndex] = RscIndex;
    }
}