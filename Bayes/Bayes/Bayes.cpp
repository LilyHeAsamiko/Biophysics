#include "stdafx.h"
#include "stdafx.cpp"
#include <math.h>
#include <vector>
#include <string>
#include <iostream>

using namespace System;


int main(array<System::String ^> ^args)
{
	//std::vector<std::string> seq;
	std::vector<std::string> seq{ "T", "S", "A", "P", "Q" };
	int n = 6;
	struct tree Tr;
	Tr = draw_tree(seq, n, Tr);  
}
