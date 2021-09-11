#ifndef CSTCHANGE_TOOLS_HPP
#include "common.h"
#define CSTCHANGE_TOOLS_HPP

void CalculateLevelList();
void CalculateDescendants();
void CalculateAncestors();
bool SortPopOnFitValueByAscend(chromosome& a, chromosome& b) ;
bool SortValueByDescend(pair<int,double>& a, pair<int,double>& b);
void IndexSort(vector<int>& ind, vector<double>& fitness);
#endif //CSTCHANGE_TOOLS_HPP
