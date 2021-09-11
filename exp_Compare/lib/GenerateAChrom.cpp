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
    chrom.Code_RK.resize(comConst.NumOfTsk);
    chrom.Code_TD.resize(comConst.NumOfRsc);
    chrom.RscAlcLst.resize(comConst.NumOfTsk);
    chrom.EndTime.resize(comConst.NumOfTsk);
}

//{generate a topological sort randomly}
vector<int> GnrSS_TS() {
    vector<int> SS;
    vector<int> upr(comConst.NumOfTsk); //the variables for recording the numbers of unscheduled parent tasks
    vector<int> RTI;                    //the set for recording ready tasks whose parent tasks have been scheduled or not exist
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i] = Tasks[i].parents.size();
        if (upr[i]==0){
            RTI.push_back(i);
        }
    }

    while (!RTI.empty()) {
        int RandVec = rand() % RTI.size();
        int v = RTI[RandVec];
        vector<int>::iterator iter = RTI.begin() + RandVec;
        RTI.erase(iter);
        for (int i = 0; i < Tasks[v].children.size(); ++i) {
            --upr[Tasks[v].children[i]];
            if (upr[Tasks[v].children[i]] == 0) RTI.push_back(Tasks[v].children[i]);
        }
        SS.push_back(v);
    }
    return SS;
}

//{generate a task scheduling order by the levels of tasks from small to large}
//{Those haveing the same level are ranked arbitrarily among them}
vector<int> GnrSS_Lvl() {
    vector<int> ch;
    vector<vector<int>> tem = TskLstInLvl;
    for (int i = 0; i < TskLstInLvl.size(); ++i) {
        random_shuffle(tem[i].begin(), tem[i].end());   //arrange the tasks in each level
        for (int j = 0; j < tem[i].size(); ++j) {
            ch.push_back(tem[i][j]);
        }
    }
    return ch;
}

double HrsDcd_EFT(chromosome& ch) {
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

// {in the SS, tasks are arranged according to the level form small to large, and the tasks in the same level are arranged in descend of the number of child tasks}
chromosome GnrChr_HEFT_Baseline() {
    chromosome TemChrom;
    IntChr(TemChrom);
    int ScheduleOrder = 0;
    for (int j = 0; j < TskLstInLvl.size(); ++j) {
        if (TskLstInLvl[j].size() < 2) {
            TemChrom.TskSchLst[ScheduleOrder++]=TskLstInLvl[j][0];
            continue;
        }
        vector<int> SonTaskNum;
        for (int i = 0; i < TskLstInLvl[j].size(); ++i)
            SonTaskNum.push_back(Tasks[TskLstInLvl[j][i]].children.size());

        vector<int> ind(TskLstInLvl[j].size());
        IndexSort(ind, SonTaskNum);
        for (int i = TskLstInLvl[j].size() - 1; i >= 0; i--) {
            TemChrom.TskSchLst[ScheduleOrder++] = TskLstInLvl[j][ind[i]];
        }
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

void InitProModelOfTskSch(vector<vector<double> >& PMS,
                      vector<int>& NumOfAncestors,
                      vector<int>& NumOfNonDescendants,
                      vector<int>& STS) {
    for(int i = 0; i < comConst.NumOfTsk; ++i) {
        int left  = NumOfAncestors[i];
        int right = NumOfNonDescendants[i];
        for(int j = left; j < right; ++j) {
            PMS[i][j] = 1;
            ++STS[j];
        }
    }

    for(int j = 0; j < comConst.NumOfTsk; ++j) {
        for(int i = 0; i < comConst.NumOfTsk; ++i) {
            PMS[i][j] = PMS[i][j] / STS[j];
        }
    }
}

void GnrRscLstOfChr(chromosome& chrom, vector<vector<double> >& PMR) {
    for (int i = 0; i < comConst.NumOfTsk; i++) {
        double rd_num = double(rand()%100) / 100;
        double sum = 0;
        for(int j = 0; j < Tasks[i].ElgRsc.size(); ++j){
            sum += PMR[i][Tasks[i].ElgRsc[j]];
            if(sum > rd_num) {
                chrom.RscAlcLst[i] = Tasks[i].ElgRsc[j];
                break;
            }
        }
    }
}

chromosome GnrTskLstOfChr(vector<vector<double> >& PMS) {
    chromosome chrom;
    IntChr(chrom);
    vector<int > upr(comConst.NumOfTsk,0);
    list<int> RTI;
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        upr[i]=Tasks[i].parents.size();
        if (upr[i]==0){
            RTI.push_back(i);
        }
    }

    for (int i = 0; i < comConst.NumOfTsk; i++) {
        double sum = 0;
        for(int j : RTI){
            sum += PMS[j][i];
        }
        vector<double> SltProb(RTI.size());
        int n = 0;
        for (int j : RTI) {
            SltProb[n] = PMS[j][i] / sum;
            n++;
        }
        double rd_num = double(rand()%100) / 100;
        double ProbSum = 0;
        int taskIndex;
        n = 0;
        for (int j : RTI) {
            ProbSum += SltProb[n];
            if(ProbSum > rd_num) {
                taskIndex = j;
                break;
            }
            n++;
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