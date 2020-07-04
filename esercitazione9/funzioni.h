#include "random.h"
#include "popolazione.h"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>
using namespace std;
void setrandom(Random &rnd);
double distanza (double x1, double y1, double x2, double y2);
void swap (double &x, double &y);
int pbc (int x, int tot);
bool cerca(int x, int *y, int c);
void sort (double *x, int n);
