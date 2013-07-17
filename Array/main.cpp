#include <iostream>
#include <string.h>

#include "Array.h"

using namespace std;
using namespace ij;

typedef int * T;

int main()
{
    INT n = 2;
    Array<T> a1(Vector<INT> (2,n));

    Array<T> a2;
    a2 = a1;

    Vector<INT> x (2); x[0] = 1; x[1] = 0;
    Vector<INT> y (2); y[0] = 0; y[1] = 1;

    for(INT i=0;i<n;i++)
    {
        for(INT j=0;j<n;j++)
        {
            a1.at(i*x+j*y) = new int;
            *a1.at(i*x+j*y) = i+j;
        }
    }

    cout << *a1.at(x) << endl;
    return 0;
}
