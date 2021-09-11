#include "common.h"
#include "tools.hpp"

using namespace std;

//{calculate the level of tasks}
void CalculateLevelList() {
    LevelIdOfTask.resize(comConst.NumOfTsk);
    vector<int> InDegree;   //variables for recording the number of parent tasks whose level have not been calculated;
    vector<int> stk;        //a set for recording the index of tasks whose inDegree is equal to 0;
    InDegree.assign(comConst.NumOfTsk, 0);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        InDegree[i] = Tasks[i].parents.size();
    }
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        if (InDegree[i] == 0) stk.push_back(i);
    }
    int MaxLevel = 0;
    while (!stk.empty()) {
        int v = stk[0];
        LevelIdOfTask[v] = 0;
        for (int i = 0; i < Tasks[v].parents.size(); ++i) {
            if (LevelIdOfTask[Tasks[v].parents[i]] >= LevelIdOfTask[v]) {
                LevelIdOfTask[v] = LevelIdOfTask[Tasks[v].parents[i]] + 1;
            }
        }
        if(LevelIdOfTask[v] + 1> MaxLevel) {
            MaxLevel = LevelIdOfTask[v] + 1;
            TskLstInLvl.resize(MaxLevel);
        }
        TskLstInLvl[LevelIdOfTask[v]].push_back(v);
        stk.erase(stk.begin());
        for (int i = 0; i < Tasks[v].children.size(); ++i) {
            InDegree[Tasks[v].children[i]]--;
            if (InDegree[Tasks[v].children[i]] == 0) {
                stk.push_back(Tasks[v].children[i]);
            }
        }
    }
}

bool SortPopOnFitValueByAscend(chromosome& a, chromosome& b) {
    return a.FitnessValue < b.FitnessValue;
}

bool SortValueByDescend(pair<int,double>& a, pair<int,double>& b) {
    return a.second > b.second;
}

void IndexSort(vector<int>& ind, vector<double>& fitness) {
    vector<double> result;
    for (int i = 0; i < ind.size(); ++i) {
        result.push_back(fitness[i]);
        ind[i] = i;
    }
    sort(ind.begin(), ind.end(), [&result](int v1, int v2) { return result[v1] < result[v2]; });
}

//计算每个任务的子孙任务
void CalculateDescendants() {
    Descendants.resize(comConst.NumOfTsk);
    for(int i = TskLstInLvl.size()-2; i >= 0; --i) {
        for(int taskId : TskLstInLvl[i]) {
            for(int childId : Tasks[taskId].children) {
                Descendants[taskId].insert(childId);
                Descendants[taskId].insert(Descendants[childId].begin(),Descendants[childId].end());
            }
        }
    }
}

//计算每个任务的祖先任务
void CalculateAncestors() {
    Ancestors.resize(comConst.NumOfTsk);
    for(int i = 1; i < TskLstInLvl.size(); ++i) {
        for(int taskId : TskLstInLvl[i]) {
            for(int parentId : Tasks[taskId].parents) {
                Ancestors[taskId].insert(parentId);
                Ancestors[taskId].insert(Ancestors[parentId].begin(),Ancestors[parentId].end());
            }
        }
    }
}

