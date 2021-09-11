//
// Created by Administrator on 2020/9/3.
//
#include "common.h"
#include "tools.hpp"
#include "config.h"
#include "GenerateAChrom.h"

double runHEFT(string XmlFile, string RscAlcFile, double& MaxTime) {
    ReadFile(XmlFile, RscAlcFile);
    CalculateLevelList();
    clock_t start = clock();
    vector<double> Rank_b(comConst.NumOfTsk, 0);
    vector<double> ww(comConst.NumOfTsk, 0);
    vector<vector<double>> cc(comConst.NumOfTsk, vector<double>(comConst.NumOfTsk,0));
    W_Cal_Average(ww);
    C_Cal_Average(cc);
    Calculate_Rank_b(Rank_b,cc,ww);
    chromosome Chrom_HEFT_b = GnrChr_HEFT(Rank_b);
    MaxTime = (double) (clock() - start) / CLOCKS_PER_SEC;
    ClearALL();
    return Chrom_HEFT_b.FitnessValue;
}