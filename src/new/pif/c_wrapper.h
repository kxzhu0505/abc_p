
#pragma once
#include <iostream>
#include <cstdio>
#include <vector>
#include <memory>

// #include "/home/kxzhu/partition/metis/include/metis.h" comment metis by zli
#include "partNtk.h"

using namespace std;

namespace ymc {

int hello();
// int try_metis(); comment metis by zli
int test_yaig();
Abc_Ntk_t* pif(Abc_Ntk_t* pNtk, uint32_t nParts, uint32_t sCluster, char* libFileName, char* dirName);
















} //for namespace



