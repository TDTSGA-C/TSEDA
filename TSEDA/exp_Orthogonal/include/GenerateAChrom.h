

#ifndef CSTCHANGE_GENERATEACHROM_H
#define CSTCHANGE_GENERATEACHROM_H

#include "common.h"
void W_Cal_Average(vector<double>& w);
void C_Cal_Average(vector<vector<double>>& c);
void Calculate_Rank_b(vector<double>& rankList, vector<vector<double>>& c, vector<double>& w);
void IntChr(chromosome& chrom);
double HrsDcd_EFT(chromosome& chrom);
chromosome GnrChr_HEFT(vector<double> rnk);
void InitProModelOfResAlc(vector<vector<double> >& PMR);
void initProModelTask(vector<vector<double> >& PTS,
                      vector<int>& ParNum,
                      vector<int>& son_num,
                      vector<int>& STS);
int ChooseTask(vector<vector<double> >& PTS, int order, vector<int>& RT);
int ChooseRes(vector<vector<double> >& PMVM, int taskId);
chromosome GnrTskLstOfChr(vector<vector<double> >& PMS);
void GnrRscLstOfChr(chromosome& chrom,
                    vector<vector<double> >& PMR) ;
#endif //CSTCHANGE_GENERATEACHROM_H
