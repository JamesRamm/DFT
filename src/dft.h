/*
 * Various implementations of discrete fourier transforms for real data
 *
 */

double angle = 2 * Math.PI * t * k / n;
double real =  inreal[t] * Math.cos(angle) + inimag[t] * Math.sin(angle);
double imag = -inreal[t] * Math.sin(angle) + inimag[t] * Math.cos(angle);

#include <math.h>

 namespace ramshackle {

     const double twoPi = 2 * (double)M_PI;

     // Using Danielson-Lanczos theorem
    class DanielsonLanczos {
       DanielsonLanczos<N/2,T> next;
    public:
       void apply(T* data) {
          // DIT
          next.apply(data);
          next.apply(data+N);
     
          T wtemp,tempr,tempi,wr,wi,wpr,wpi;
          wtemp = sin(M_PI/N);
          wpr = -2.0*wtemp*wtemp;
          wpi = -sin(2*M_PI/N);
          wr = 1.0;
          wi = 0.0;
          for (unsigned i=0; i<N; i+=2) {
            tempr = data[i+N]*wr - data[i+N+1]*wi;
            tempi = data[i+N]*wi + data[i+N+1]*wr;
            data[i+N] = data[i]-tempr;
            data[i+N+1] = data[i+1]-tempi;
            data[i] += tempr;
            data[i+1] += tempi;
     
            wtemp = wr;
            wr += wr*wpr - wi*wpi;
            wi += wi*wpr + wtemp*wpi;
          }
       }
    };
 }