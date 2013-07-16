#include <iostream>
#include "Array.h"

using namespace std;

typedef unsigned int INT;

int main()
{
    int **i = new int *[3];
    for(int j=0;j<3;j++) i[j] = new int;

    *i[0] = 42;
    cout << *i[0] << endl;

    for(int j=0;j<3;j++) delete i[j];
    delete [] i;

    return 0;
}
