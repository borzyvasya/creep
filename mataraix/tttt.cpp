#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <string>
#include <math.h>
#include <float.h>
#include <algorithm>

using namespace std;
const  int n = 3, p = 10, MAXITER = 10000;
const  double eps = 0.0001, one = 1.0, zero = 0.0;
void ZeidelIt(double[], double[][n], double[], int);
void ShowVector(double[], int, string);
void ShowMatrix(double[][n], int, int, string);
double MAXVAL(double[], double[], int);
double norm(double[][n], int);
bool check(double[][n], int); // check for diagonal domi

int main()
{
    int msz = n, nsz = n, i, j, k = 0;
    double B[n][n], c[n], x[n], xtmp[n];
    double A[n][n] = { 3.,1.,-1.,2.,6.,3.,-1.,1.,4. };
    double vb[n] = { 7.,-2.,4. };
    double eps2, criteria = one, normb, w = 1.04;
    for (i = 0; i < n; i++)x[i] = zero; for (i = 0; i <
        n; i++)xtmp[i] = zero;
    for (i = 0; i < n; i++)c[i] = zero; for (i = 0; i <
        n; i++)for (j = 0; j < n; j++)B[i][i] = zero;
    cout << " Matrix & Vector data...." << endl;
    ShowMatrix(A, msz, nsz, " Matrix A data : ");
    cout << endl;
    ShowVector(vb, nsz, " Vector b data : ");
    for (i = 0; i < msz; i++) 
        for (j = 0; j < nsz; j++)
            if (i != j)B[i][j] = -A[i][j] / A[i][i];

    for (i = 0; i < nsz; i++) c[i] = vb[i] / A[i][i];
    ShowMatrix(B, msz, nsz, " Matrix B data : ");
    ShowVector(c, nsz, " Vector c data : ");

    normb = norm(B, msz);
    if (!check(A, nsz))
        cout << endl << "Matrix[A]: a full diagonal predomi nance was not found!" << endl;
    else
        cout << endl << "Matrix[A]: is a full diagonal predominance" << endl;
    cout << endl << "MATRIX Norm for | B | = " << normb <<endl;
    if (normb >= one)
        eps2 = eps;

    else
        eps2 = eps * (one - normb) / normb;
    cout << endl << endl << " SOLUTION OF THE SLAE by the SOR method using, with eps = " 
        << eps2 << ", " << "w = " << w << endl;
    cout << "___________________________________" << endl << endl;
    do
    {
        ZeidelIt(x, B, c, nsz);
        for (i = 0; i < nsz; i++) xtmp[i] = w * x[i] + (one - w) * xtmp[i];
        criteria = fabs(MAXVAL(x, xtmp, nsz));
        for (i = 0; i < nsz; i++)xtmp[i] = x[i];
        k += 1;
        cout << setw(p) << scientific << setprecision(2) <<
            "SOLUTION(" << k << ")" << ", CRITERIA = " << criteria <<
            endl;
    } while ((criteria > eps2 || criteria == zero) && k <=
        MAXITER);
    if (k > MAXITER)
    {
        cerr << "FAILED TO CONVERGENCE : maximum iterations was exceeded!" << endl;
    }
    else
    {
        cout << endl << " Finished!" << endl;
        cout << " Iterations was made = " << k << endl;
        cout << " Exit criteria = " << criteria << endl;
        cout << " Solution is : " << endl << endl;
        for (i = 0; i < nsz; i++)cout << " x(" << setw(2) <<
            i + 1 << ") = " << fixed << setprecision(3) << setw(p) << x[i]
            << endl;
    }
    cin.get();
    return EXIT_SUCCESS;
}
void ShowMatrix(double Matrix[][n], int msz, int nsz,
    string str)
{
    cout << str << endl;
    int i, j;
        for (i = 0; i < msz; ++i)
        {
            cout << endl;
            for (j = 0; j < nsz; ++j) cout << " " << setw(p) <<
                fixed << setprecision(4) << Matrix[i][j];
        }
    cout << endl;
}
void ShowVector(double Vector[], int nsz, string str)
{
    cout << str << endl;
    for (int j = 0; j < nsz; j++)cout << " " << setw(p) <<
        fixed << setprecision(4) << Vector[j];
    cout << endl;
}
void ZeidelIt(double x[], double B[][n], double c[], int
    nsz)
{
    int i, j;
    double sum;
    for (i = 0; i < nsz; i++)
    {
        sum = c[i];
        for (j = 0; j < nsz; j++)if (i != j)sum = sum + B[i]
            [j] * x[j];
            x[i] = sum;
    }
}
double MAXVAL(double v1[], double v2[], int sz)
{
    int i;
    double maxval = 0.0, * tmp;
    tmp = new double[sz];
    for (i = 0; i < sz; i++)tmp[i] = v1[i] - v2[i];
    maxval = *max_element(tmp, tmp + sz);
    delete[] tmp;
    return maxval;
}
double norm(double A[][n], int m)
{
    int i, j;
    double fnorm = zero, * tmp, sum;
    tmp = new double[m];
    for (i = 0; i < m; ++i)
    {
        sum = zero;
        for (j = 0; j < m; ++j)sum = sum + fabs(A[i][j]);
        tmp[i] = sum;
    }
    fnorm = *max_element(tmp, tmp + m);
    delete[] tmp;
    return fnorm;
}
bool check(double A[][n], int nsz)
{
    int i, j, isdiag = true;
    double sum;
    for (i = 0; i < nsz; i++)
    {
        sum = zero;
        for (j = 0; j < nsz; j++) sum += fabs(A[i][j]);
        sum = sum - fabs(A[i][i]);
        if (sum >= A[i][i])
        {
            isdiag = false;
            break;
        }
    }
    return isdiag;
}