
#pragma once
#include <iostream>
#include <cstdio>
#include <vector>
#include <memory>

#include <metis.h>
#include "partNtk.h"
#include "omp.h"

using namespace std;

namespace ymc {

int hello();
int try_metis();
int test_yaig();
Abc_Ntk_t* pif(Abc_Ntk_t* pNtk, uint32_t nParts, char* libFileName);
















} //for namespace



