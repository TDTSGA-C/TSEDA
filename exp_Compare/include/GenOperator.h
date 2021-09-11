//
// Created by smallfish on 18-4-4.
//

#ifndef CSTCHANGE_CROSSOVER_H
#include "common.h"
#define CSTCHANGE_CROSSOVER_H
double DcdEvl(chromosome& chrom, bool isForward);
double IFBDI(chromosome& chrom);
void LBCAI(chromosome& chrom);
void CalSlctProb_Rank(double RtOfSltPrb, vector<double>& A,int& NumOfChormPerPop);
int SltChr(vector<double>& A);
void SelectionTournament(int& parent_1, int& parent_2 ,int& NumOfChormPerPop);
void CrsMS_MP(chromosome& chrom1, chromosome& chrom2);
void MtnSS_TS(chromosome& a);
void MtnMS_MP(chromosome& ch);
void MtnMS_SP(chromosome& a);
void CrsSS_TS_R(chromosome& ch1, chromosome& ch2, int CrossPoint);
void CrsMS_SP(chromosome& ch1, chromosome& ch2);
void Crossover_CGA(chromosome& pop1, chromosome& pop2);
void Mutation_CGA(chromosome& a);
void CrsMS_DP(chromosome& ch1, chromosome& ch2);
void GnrTskSchLst_HGA(chromosome& chrom);
void Crossover_HGA(chromosome& chrom1, chromosome& chrom2);
void Mutation_HGA(chromosome& chrom);
void Crs_Lvl(chromosome& chrom1, chromosome& chrom2);
void CrsSS_ExcTskInLvl(chromosome& chrom1, chromosome& chrom2);
void Crossover_LWSGA(chromosome& chrom1, chromosome& chrom2);
void MtnSS_ExcTskInLvl(chromosome& chrom);
void Mtn_rebuild_level(chromosome& ch);
void Mutation_LWSGA(chromosome& chrom);
void RscLoadAdjust_HGA(vector<chromosome>& Pop);
#endif //CSTCHANGE_CROSSOVER_H
