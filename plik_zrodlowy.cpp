//przedzial <150;245> funkcja p0*x^7+p1*x^6+p2*x^5+p3*x^4+p4*x^3+p5*x^2+p6*x+p7*c
//zmienne: 
//EXT PARAMETER                APPROXIMATE        STEP         FIRST
// NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
// 1  p0           3.48899e-14   5.16012e-16   8.84243e-20  -9.83238e+11
// 2  p1          -5.19831e-12   1.34129e-13   2.02847e-17   6.35468e+09
// 3  p2          -1.93046e-09   3.18631e-11   1.30109e-15   1.46084e+07
// 4  p3          -2.13845e-07   7.20047e-09   2.90615e-13  -7.63812e+04
// 5  p4           4.05831e-05   1.57636e-06   6.36610e-11  -5.45408e+02
// 6  p5           2.25273e-02   3.33389e-04   1.36153e-08   1.79462e-01
// 7  p6           3.38742e+00   6.56963e-02   2.83251e-06   1.66506e-02
// 8  p7          -9.93351e+02   1.08143e+01   5.72160e-04  -4.11678e-05

//przedzial <245;350> gauss
//zmienne:
//EXT PARAMETER                                   STEP         FIRST
//  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
// 1  Constant     1.67624e+02   3.92727e+00   8.49572e-03   1.02554e-05
// 2  Mean         2.97723e+02   6.32400e-01   1.83293e-03   4.11436e-04
// 3  Sigma        2.97756e+01   6.89059e-01   1.85642e-05   3.43117e-02

//przedzial <350;500> pol5 (wielomian piatego stopnia)
//zmienne:
//Chi2                      =      23.1953
//NDf                       =           34
//p0                        =     -706.373   +/-   70257.8
//p1                        =      8.62646   +/-   825.052
//p2                        =   -0.0454071   +/-   3.86228
//p3                        =  0.000132775   +/-   0.00900983
//p4                        = -2.06272e-07   +/-   1.04744e-05
//p5                        =  1.29316e-10   +/-   4.85515e-09
#include <map>
#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <random>
#include <iomanip>
#include <array>
using namespace std;
# define M_PI           3.14159265358979323846
piecewise_linear_distribution<double> triangular_distribution(double min, double peak, double max)
{
    array<double, 3> i{ min, peak, max };
    array<double, 3> w{ 0, 1, 0 };
    return piecewise_linear_distribution<double>{i.begin(), i.end(), w.begin()};
}
int main()
{
    ofstream wynik;
    wynik.open("wynik.txt");
    srand(time(NULL));
    for (int i = 0; i < 5000; i++)
    {


        double czas = 0, ludzie_x, ludzie_p, d, x, s = 0, lam1 = 4, y, a;
        int x1 = -1;
        double c = 20;
        double b = 0.2;

        random_device rd;
        mt19937 gen(rd());
        //Rotacja ludzi: - rozklad poissona
        while (s <= lam1)
        {
            d = (double)rand() / RAND_MAX;
            y = -log(1 - d);
            s += y;
            x1++;
        }
        //rozklad trojkatny przy uzyciu rozkladu liniowego
        if (x1 <= 2)
        {
            auto dist = triangular_distribution(1, 3, 6);
            x = dist(gen);
            czas += x;
        }
        else if (x1 <= 4)
        {
            auto dist = triangular_distribution(3, 7, 11);
            x = dist(gen);
            czas += x;
        }
        else if (x1 <= 6)
        {
            auto dist = triangular_distribution(10, 13, 16);
            x = dist(gen);
            czas += x;
        }
        else if (x1 <= 8)
        {
            auto dist = triangular_distribution(15, 19, 23);
            x = dist(gen);
            czas += x;
        }
        else
        {
            auto dist = triangular_distribution(20, 25, 30);
            x = dist(gen);
            czas += x;
        }
        //wyjazd z przystanku:
        //funkcja f(x)=cos(1/12x-pi/6) -przedzial <0;6pi>
        d = (double)rand() / RAND_MAX;
        x = 12 * (asin(d / (sqrt(3) - 1) - 1. / 2) + M_PI / 6.);
        czas += x;
        //dojazd do swiatel:
        d = (double)rand() / RAND_MAX;
        if (d == 1)
        {
            while (d == 1)
                d = (double)rand() / RAND_MAX;
        }
        //wykladnicza
        x = c - (1 / b) * log(1 - d);
        czas += x;
        //swiatla
        d = (double)rand() / RAND_MAX;
        //zielone swiatlo - rozklad normalny
        if (d <= 0.2142)
        {
            normal_distribution<> nd{ 15,2 };
            czas += nd(gen);
        }
        //zolte ; zielone
        else if (d <= 0.2357)
        {
            normal_distribution<> nd{ 18,2 };
            czas += nd(gen);
        }
        //zolte ; czerwone
        else if (d <= 0.2571)
        {
            normal_distribution<> nd{ 120,15 };
            czas += nd(gen);
        }
        //czerwone 
        else if (d <= 0.9)
        {
            normal_distribution<> nd{ 100,10 };
            czas += nd(gen);
        }
        //podwojne czerwone ; korki
        else
        {
            normal_distribution<> nd{ 200,20 };
            czas += nd(gen);
        }

        //droga prosta:
        //funkcja cos(1/6x-pi/6) przedzial <12pi;16pi>
        d = (double)rand() / RAND_MAX;
        x = 10 * M_PI + 6 * acos(1. / 2 - (3 * d) / 2.);
        czas += x;
        //zakret z pierwszenstwem:


        //spokojny przejazd - box miller
        double mi = 20, sig = 3, u1, u2, z, r, m, suma = 0;
        u1 = (double)rand() / RAND_MAX;
        u2 = (double)rand() / RAND_MAX;
        r = sqrt(-2 * (log(u1)));
        z = r * sin(2 * M_PI * u2);
        x = z * sig + mi;
        //wymuszenie przez kogos - metoda w oparciu o centralne twierdzenie graniczne
        d = (double)rand() / RAND_MAX;
        if (d > 0.85)
        {
            m = 6;
            mi = 0.1;
            sig = sqrt(1. / 24.);
            for (int i = 0; i < m; i++)
            {
                x = (double)rand() / RAND_MAX;
                suma += x;

            }
            y = (((1 / m) * suma) - mi) / ((sig) / (sqrt(m)));
            if (y < 0) {
                y = -y;
            }
            czas += y;
        }
        //rondo:
        //liczba samochodow przed nami - bernoulliego
        int k = 0, lprob = 5;
        double prawd = 0.4;
        for (int i = 0; i < lprob; i++)
        {
            r = (double)rand() / RAND_MAX;
            if (r <= prawd)
            {

                k++;
            }
        }

        //wjezdzamy natychmiast - rozklad liniowy
        if (k == 0)
        {
            a = 4.5;
            b = 10;
            d = (double)rand() / RAND_MAX;
            x = d * (b - a) + a;
        }
        //po 1 samochodzie
        else if (k == 1)
        {
            a = 7;
            b = 14;
            d = (double)rand() / RAND_MAX;
            x = d * (b - a) + a;
        }
        //po 2
        else if (k == 2)
        {
            a = 12;
            b = 22;
            d = (double)rand() / RAND_MAX;
            x = d * (b - a) + a;
        }
        //po 3
        else if (k == 3)
        {
            a = 20;
            b = 33;
            d = (double)rand() / RAND_MAX;
            x = d * (b - a) + a;
        }
        //po 4
        else if (k == 4)
        {
            a = 30;
            b = 46;
            d = (double)rand() / RAND_MAX;
            x = d * (b - a) + a;
        }
        //po 5
        else if (k == 5)
        {
            a = 43;
            b = 60;
            d = (double)rand() / RAND_MAX;
            x = d * (b - a) + a;
        }
        czas += x;
        //kolejne swiatla - rozklad trojkatny:
        d = (double)rand() / RAND_MAX;
        //zielone swiatlo
        if (d <= 0.625)
        {
            auto dist = triangular_distribution(20, 29, 36);
            x = dist(gen);
        }
        //zolte ; zielone
        else if (d <= 0.65625)
        {
            auto dist = triangular_distribution(25, 33, 40);
            x = dist(gen);
        }
        // czerwone
        else if (d <= 0.96875)
        {
            auto dist = triangular_distribution(55, 64, 70);
            x = dist(gen);
        }
        //zolte ; czerwone 
        else
        {
            auto dist = triangular_distribution(63, 70, 78);
            x = dist(gen);
        }
        czas += x;
        //przyjazd - wykladnicza:
        c = 45;
        b = 1. / 10;
        d = (double)rand() / RAND_MAX;
        x = c - (1 / b) * log(1 - d);
        czas += x;
        wynik << czas<<endl;
    }
}
