//
// Created by Administrator on 2021/5/9.
//

#include "TSEDA.h"
#include "tools.hpp"
#include "config.h"
#include "GenerateAChrom.h"
#include "GenOperator.h"

double runTSEDA(string XmlFile, string RscAlcFile, Orthogonal& orthogonal,double& SchTime, int& iteration) {
    clock_t start = clock();
    ReadFile(XmlFile, RscAlcFile);
    ConfigParameter_TSEDA(orthogonal);
    CalculateLevelList();
    CalculateDescendants();
    CalculateAncestors();

    vector<double> Rank_b(comConst.NumOfTsk, 0);
    vector<double> ww(comConst.NumOfTsk, 0);
    vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk, 0));
    W_Cal_Average(ww);
    C_Cal_Average(cc);
    Calculate_Rank_b(Rank_b, cc, ww);

    //初始化任务调度顺序概率模型
    vector<int> NumOfAncestors(comConst.NumOfTsk);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        NumOfAncestors[i] = Ancestors[i].size();
    }
    //递归计算，子孙任务的数量
    vector<int> NumOfDescendants(comConst.NumOfTsk);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        NumOfDescendants[i] = Descendants[i].size();
    }
    //总任务数-子孙任务数量-1
    vector<int> NumOfNonDescendants(comConst.NumOfTsk);
    for (int i = 0; i < comConst.NumOfTsk; ++i) {
        NumOfNonDescendants[i] = comConst.NumOfTsk -  NumOfDescendants[i] - 1;
    }

    vector<vector<double> > PMR(comConst.NumOfTsk, vector<double>(comConst.NumOfRsc));
    InitProModelOfResAlc(PMR);
    vector<int> STS(comConst.NumOfTsk ,0);
    vector<vector<double> > PMS(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk ,0));
    initProModelTask(PMS, NumOfAncestors, NumOfNonDescendants, STS);

    chromosome chrom_HEFT_b = GnrChr_HEFT(Rank_b);         //generate a chromosome according to HEFT
    NewPopulation.push_back(chrom_HEFT_b);
    double BestFitness = 999999999999;
    while ((double) (clock() - start) / CLOCKS_PER_SEC < Parameter_TSEDA.RunTimeRatioOfStg1 * SchTime) {
        while (NewPopulation.size() < Parameter_TSEDA.NumOfChormPerPop) {
            chromosome chrom = GnrTskLstOfChr(PMS);
            HrsDcd_EFT(chrom);
            NewPopulation.push_back(chrom);
        }
        population.insert(population.end(), NewPopulation.begin(), NewPopulation.end());
        sort(population.begin(), population.end(), SortPopOnFitValueByAscend);
        vector<chromosome>::iterator Iter = population.begin();
        advance(Iter, Parameter_TSEDA.NumOfChormPerPop);
        population.assign(population.begin(), Iter);
        if (population[0].FitnessValue + PrecisionValue < BestFitness) {
            BestFitness = population[0].FitnessValue;
        }
        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            for (int j = 0; j < comConst.NumOfRsc; ++j) {
                int count = 0;
                for (int n = 0; n < Parameter_TSEDA.NumOfEliteOfPop; ++n) {
                    if (population[n].RscAlcLst[i] == j) {
                        ++count;
                    }
                }
                PMR[i][j] = (1 - Parameter_TSEDA.alpha) * PMR[i][j] +
                            Parameter_TSEDA.alpha * (1.0 * count / Parameter_TSEDA.NumOfEliteOfPop);
            }
        }

        for (int i = 0; i < comConst.NumOfTsk; ++i) {
            for (int j = 0; j < comConst.NumOfTsk; ++j) {
                int count = 0;
                for (int n = 0; n < Parameter_TSEDA.NumOfEliteOfPop; ++n) {
                    if (population[n].TskSchLst[i] == j) {
                        ++count;
                    }
                }
                PMS[j][i] = (1 - Parameter_TSEDA.beta) * PMS[j][i] +
                            Parameter_TSEDA.beta * (1.0 * count / Parameter_TSEDA.NumOfEliteOfPop);
            }
        }
        ++iteration;
        NewPopulation.clear();
    }
    while ((double) (clock() - start) / CLOCKS_PER_SEC < SchTime){
        while(NewPopulation.size() < Parameter_TSEDA.NumOfChormPerPop){
            chromosome chrom = GnrTskLstOfChr(PMS);
            GnrRscLstOfChr(chrom,PMR);
            DcdEvl(chrom,true);
            NewPopulation.push_back(chrom);
        }
        sort(NewPopulation.begin(), NewPopulation.end(), SortPopOnFitValueByAscend);
        #pragma omp parallel for
        for (int n = 0; n < Parameter_TSEDA.NumOfImproveOfPop; ++n) {
            IFBDI(NewPopulation[n]);
            LBCAI(NewPopulation[n]);
        }
        population.insert(population.end(), NewPopulation.begin(), NewPopulation.end());
        sort(population.begin(), population.end(), SortPopOnFitValueByAscend);
        vector<chromosome>::iterator Iter = population.begin();
        advance(Iter, Parameter_TSEDA.NumOfChormPerPop);
        population.assign(population.begin(), Iter);
        if ( population[0].FitnessValue + PrecisionValue <  BestFitness ) {
            BestFitness = population[0].FitnessValue;
        }
        for(int i = 0; i < comConst.NumOfTsk; ++i) {
            for(int j = 0; j < comConst.NumOfRsc; ++j) {
                int count = 0;
                for(int n = 0; n < Parameter_TSEDA.NumOfEliteOfPop; ++n) {
                    if(population[n].RscAlcLst[i] == j) {
                        ++count;
                    }
                }
                PMR[i][j] = (1-Parameter_TSEDA.alpha)*PMR[i][j] + Parameter_TSEDA.alpha*(1.0*count/Parameter_TSEDA.NumOfEliteOfPop);
            }
        }

        for(int i = 0; i < comConst.NumOfTsk; ++i) {
            for(int j = 0; j < comConst.NumOfTsk; ++j) {
                int count = 0;
                for(int n = 0; n < Parameter_TSEDA.NumOfEliteOfPop; ++n) {
                    if(population[n].TskSchLst[i] == j) {
                        ++count;
                    }
                }
                PMS[j][i] = (1-Parameter_TSEDA.beta)*PMS[j][i] + Parameter_TSEDA.beta*(1.0*count/Parameter_TSEDA.NumOfEliteOfPop);
            }
        }
        ++iteration;
        NewPopulation.clear();
    }
    ClearALL();
    SchTime = (double) (clock() - start) / CLOCKS_PER_SEC;
    return BestFitness;
}