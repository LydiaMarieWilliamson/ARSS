#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#ifdef WIN32
#   include "<windows.h>" // For Now().
#else
#   include <sys/time.h> // For Now().
#endif
#include "EnBMP.h"
#include "DeBMP.h"
#include "EnWAV.h"
#include "DeWAV.h"

typedef long Integer;

// The Utility Routines.
void CheckOut(bool Fatal, char *Format, ...) {
   va_list AP; va_start(AP, Format); vfprintf(stderr, Format, AP); va_end(AP);
   fputc('\n', stderr);
   if (Fatal) fprintf(stderr, "Exiting with error.\n"), exit(EXIT_FAILURE);
}

Integer Now(void) { // In milliseconds.
#ifdef WIN32
   return (Integer)GetTickCount();
#else
   struct timeval T; gettimeofday(&T, NULL);
   return (Integer)1000*T.tv_sec + T.tv_usec/1000;
#endif
}

inline Integer RoundDn(double X) { double Xr = fmod(X, 1.0); if (Xr >= 0.0) Xr -= 1.0; return (Integer)(X - Xr - 1.0); }
inline Integer RoundTo(double X) { return X < 0.0? (Integer)(X - 0.5): X > 0.0? (Integer)(X + 0.5): 0; }
inline Integer RoundUp(double X) { double Xr = fmod(X, 1.0); if (Xr <= 0.0) Xr += 1.0; return (Integer)(X - Xr + 1.0); }

double GetNumber(void) {
   char Buf[0x20]; fgets(Buf, sizeof Buf - 1, stdin);
   return Buf[0] == '\0'? 0.0: atof(Buf);
}

// The next number expressible as a product of small primes.
inline Integer NextUnPrime(Integer X) {
   if (X < 1) X = 1;
   else if (X > 1) while (true) {
      Integer Y = X;
      for (; Y%2 == 0; Y /= 2);
      for (; Y%3 == 0; Y /= 3);
      if (Y == 1) break;
      X++;
   }
   return X;
}

double LogBase = 2.0; // Base for hybrid exponential/linear scaling. Anything other than 2 isn't fully supported yet.
inline double LogB(double X) {
// if (X == 0) fprintf(stderr, "Warning: log(0) returning -infinite\n");
   return LogBase == 1.0? X: log(X)/log(LogBase);
}

inline double RandNum(void) { // A uniformly-distributed random number ∈ [-1.0, +1.0).
#if RAND_MAX == 0x7fffffff
   int D = rand();
#elif RAND_MAX == 0x7fff
   int D = ((rand()&0xff) << 24) | ((rand()&0xff) << 16) | ((rand()&0xff) << 8) | (rand()&0xff);
#else
#   error Unhandled maximum random number value RAND_MAX. Please signal this error to the developer.
   int D = rand();
#endif
   return (double)D/0x40000000 - 1.0;
}

char *GetString(void) {
   char Buf[FILENAME_MAX]; fgets(Buf, sizeof Buf, stdin);
   size_t N = strlen(Buf);
   char *NewS = malloc(N*sizeof *NewS);
   for (size_t n = 0; n < N; n++) NewS[n] = Buf[n];
   NewS[N - 1] = '\0';
   return NewS;
}

// The DSP Routines.
// DFT Sub-Module: Begin.
const double Pi = 3.1415926535897932384626433832795029;
const double TwoPi = 6.28318530717958647692528676655900577;

// Support for complex numbers.
// Replace this with the compiler-native complex arithmetic, if available.
typedef struct Complex { double X, Y; } Complex;

static inline void ZeroC(Complex *CP) { CP->X = 0.0, CP->Y = 0.0; }
static inline void SetR(Complex *CP, double X) { CP->X = X, CP->Y = 0.0; }
static inline void SetC(Complex *CP, double X, double Y) { CP->X = X, CP->Y = Y; }
static inline void SetOne(Complex *CP, double Omega) { CP->X = cos(Omega), CP->Y = sin(Omega); }

static inline double AbsC(Complex A) {
   double Ax = A.X, Ay = A.Y;
   return sqrt(Ax*Ax + Ay*Ay);
}

static inline void AddMul(Complex *CP, Complex A, Complex B) {
   double Ax = A.X, Ay = A.Y, Bx = B.X, By = B.Y;
   CP->X += Ax*Bx - Ay*By, CP->Y += Ax*By + Ay*Bx;
}

static inline void Mul(Complex *CP, Complex A, Complex B) {
   double Ax = A.X, Ay = A.Y, Bx = B.X, By = B.Y;
   CP->X = Ax*Bx - Ay*By, CP->Y = Ax*By + Ay*Bx;
}

static unsigned Factor(unsigned *NP) { // Extract factors from *NP.
   unsigned N = *NP;
   if (N < 1) { *NP = 1; return N; }
   else if (N%4 == 0) { *NP = N/4; return 4; }
// else if (N%6 == 0) { *NP = N/6; return 6; }
   else if (N%2 == 0) { *NP = N/2; return 2; }
   else if (N%3 == 0) { *NP = N/3; return 3; }
   for (unsigned n = 5, n2 = n*n; n2 <= N; n++, n2 += 4*n++) if (N%n == 0) { *NP = N/n; return n; }
   *NP = 1; return N;
}

// In-situ transform (X[n]: n ∈ N) ↔ (X[k]: k ∈ N).
// Real layout: X[n] = Z_n, for n ∈ N; with Z_n* = Z_n.
// Co-Real layout: X[k] = Re(z_k) for 0 ≤ k ≤ N/2, X[N - k] = Im(z_k), for 0 < k < N/2; with z_n* = z_{N-n}.

// DFT encode the N-vector X: Real -> Co-Real.
void EnDFT(double *X, unsigned N) {
// fftw_plan Plan = fftw_plan_r2r_1d(N, X, X, 0, FFTW_ESTIMATE);
// fftw_execute(Plan), fftw_destroy_plan(Plan);
   if (N <= 1) return; // DFT's of size 1 or less are treated as the identity transform.
// Powers of unity (i.e. "Phasor"): W[n] = 1^{n/N}
   Complex *W = malloc(N*sizeof *W); if (W == NULL) fprintf(stderr, "Out of memory for W.\n"), exit(1);
// Decoded vector complexified in Real format
   Complex *Z = malloc(N*sizeof *Z); if (Z == NULL) fprintf(stderr, "Out of memory for Z.\n"), exit(1);
// Encoded vector complexified in Co-Real format
   Complex *z = malloc(N*sizeof *z); if (z == NULL) fprintf(stderr, "Out of memory for z.\n"), exit(1);
   double Theta = TwoPi/(double)N;
   for (unsigned n = 0; n < N; n++) SetOne(&W[n], n*Theta);
   for (unsigned n = 0; n < N; n++) SetR(&Z[n], X[n]);
   unsigned Q = 1;
   for (unsigned P = N; P > 1; ) {
      unsigned R = Factor(&P);
      if (Q > 1) {
         unsigned n = 0;
         for (unsigned q = 0; q < Q; q++) for (unsigned r = 0; r < R; r++) for (unsigned p = 0; p < P; n++, p++)
            Mul(&Z[n], W[(P*q*r)%N], z[n]);
      }
      unsigned PQ = P*Q, PR = P*R;
      unsigned n = 0;
      for (unsigned r = 0; r < R; r++) for (unsigned q = 0; q < Q; q++) for (unsigned p = 0; p < P; n++, p++) {
         ZeroC(&z[n]);
         for (unsigned s = 0; s < R; s++) AddMul(&z[n], W[(PQ*r*s)%N], Z[PR*q + P*s + p]);
      }
      Q *= R;
   }
// Load X from z here with the Co-Real format.
   unsigned Nq2 = N/2, Nr2 = N%2;
   X[0] = z[0].X;
   for (unsigned n = 1; n < Nq2; n++) X[n] = z[n].X, X[N - n] = z[n].Y;
   if (Nr2 == 0) X[Nq2] = z[Nq2].X;
// assert(Z[0].Y == 0.0);
// for (unsigned n = 1; n < Nq2; n++) assert(z[n].X == z[N - n].X && z[n].Y == -z[N - n].Y);
// if (Nr2 == 0) assert(z[Nq2].Y == 0.0);
   free(z), free(Z), free(W);
}

// DFT decode the N-vector X: Co-Real -> Real. This is normalized.
void DeDFT(double *X, unsigned N, bool Analytic) {
// fftw_plan Plan = fftw_plan_r2r_1d(N, X, X, 1, FFTW_ESTIMATE);
// fftw_execute(Plan), fftw_destroy_plan(Plan);
   if (N <= 1) return; // DFT's of size 1 or less are treated as the identity transform.
// Powers of unity (i.e. "Phasor"): W[n] = 1^{n/N}
   Complex *W = malloc(N*sizeof *W); if (W == NULL) fprintf(stderr, "Out of memory for W.\n"), exit(1);
// Encoded vector complexified in Co-Real format
   Complex *z = malloc(N*sizeof *z); if (z == NULL) fprintf(stderr, "Out of memory for z.\n"), exit(1);
// Decoded vector complexified in Real format
   Complex *Z = malloc(N*sizeof *Z); if (Z == NULL) fprintf(stderr, "Out of memory for Z.\n"), exit(1);
   double Theta = TwoPi/(double)N;
   for (unsigned n = 0; n < N; n++) SetOne(&W[n], -(int)n*Theta);
// Load z from X here with the Co-Real format.
   unsigned Nq2 = N/2, Nr2 = N%2;
   SetR(&z[0], X[0]);
   if (Analytic)
      for (unsigned n = 1; n < Nq2; n++) SetC(&z[n], 2.0*X[n], 2.0*X[N - n]), ZeroC(&z[N - n]);
   else
      for (unsigned n = 1; n < Nq2; n++) SetC(&z[n], X[n], X[N - n]), SetC(&z[N - n], X[n], -X[N - n]);
   if (Nr2 == 0) SetR(&z[Nq2], X[Nq2]);
   unsigned P = 1;
   for (unsigned Q = N; Q > 1; ) {
      unsigned R = Factor(&Q), PQ = P*Q;
      unsigned n = 0;
      for (unsigned q = 0; q < Q; q++) for (unsigned r = 0; r < R; r++) for (unsigned p = 0; p < P; n++, p++) {
         ZeroC(&Z[n]);
         for (unsigned s = 0; s < R; s++) AddMul(&Z[n], W[(PQ*r*s)%N], z[PQ*s + P*q + p]);
      }
      if (Q > 1) {
         unsigned n = 0;
         for (unsigned q = 0; q < Q; q++) for (unsigned r = 0; r < R; r++) for (unsigned p = 0; p < P; n++, p++)
            Mul(&z[n], W[(P*q*r)%N], Z[n]);
      }
      P *= R;
   }
   if (Analytic) for (unsigned n = 0; n < N; n++) X[n] = AbsC(Z[n])/N;
   else for (unsigned n = 0; n < N; n++) X[n] = Z[n].X/N;
   free(Z), free(z), free(W);
}
// DFT Sub-Module: End.

// Normalize the amplitude of XY to the range [-Norm, +Norm].
void Normalize(double **XY, Integer Xs, Integer Ys, double Norm) {
   double Amp = 0.0;
   for (Integer Y = 0; Y < Ys; Y++) for (Integer X = 0; X < Xs; X++) if (Amp < fabs(XY[Y][X])) Amp = fabs(XY[Y][X]);
// printf("norm : %.3f\n", Amp);
   if (Amp > 0.0) Amp = Norm/Amp; else return;
   for (Integer Y = 0; Y < Ys; Y++) for (Integer X = 0; X < Xs; X++) XY[Y][X] *= Amp;
// printf("ex : %.3f\n", XY[0][0]);
}

// Convert a relative scale coordinate P ∈ [0, 1] (i.e. band number/band count) into a frequency F ∈ [Lo, Hi].
double EnNote(double P, double Lo, double Hi) {
   if (LogBase == 1.0) return Lo + (Hi - Lo)*P; // Linear scale.
   double K = log(LogBase)*(log(Hi) - log(Lo))/log(2.0);
   return Lo + (Hi - Lo)*(exp(K*P) - 1.0)/(exp(K) - 1.0);
// For LogBase ≅ 1.0, the value is Lo + (Hi - Lo) P (1 + K (P - 1)/2 + O(K²)).
}

// Convert a frequency F ∈ [Lo, Hi] into a relative scale coordinate (i.e. band number/band count) P ∈ [0, 1].
double DeNote(double F, double Lo, double Hi) {
   double Z = (F - Lo)/(Hi - Lo); if (LogBase == 1.0) return Z;
   double K = log(LogBase)*(log(Hi) - log(Lo))/log(2.0);
   return log(1.0 + (exp(K) - 1.0)*Z)/K;
// return log(Lo*(1.0 + (exp(K) - 1.0)*Z)/log(LogBase))*log(LogBase)/K; // The original code up to 0.2.3.
// For LogBase ≅ 1.0, the value is F (1 + K (1 - F)/2 + O(K²)).
}

double *MakePTab(double LoP, double HiP, Integer Ps) {
// Tabulate the central frequencies of the bands.
   double *PTab = malloc(Ps*sizeof *PTab);
   for (Integer P = 0; P < Ps; P++) PTab[P] = EnNote((double)P/(Ps - 1), LoP, HiP);
// TODO: Change the sampling rate instead, if the upper frequency limit is too high.
   if (EnNote((double)Ps/(Ps - 1), LoP, HiP) > 0.5) printf("Warning: The upper frequency limit above the Nyquist frequency.\n");
   return PTab;
}

// The Blackman function and its (normalized) antideriative.
// const int Black0 = 7938, Black1 = 9240, Black2 = 1430;
const int Black0 = 21, Black1 = 25, Black2 = 4;

double Blackman(double X) {
   const int Black012 = Black0 + Black1 + Black2;
// return 0.42 + 0.50*cos(Pi*X) + 0.08*cos(TwoPi*X);
   return (Black0 + Black1*cos(Pi*X) + Black2*cos(TwoPi*X))/Black012;
}

// The anti-derivative of the Blackman function.
double BlackmanInt(double X) {
   double Omega = TwoPi*X, Cs = cos(Omega), Sn = sin(Omega);
// return X - Sn*(0.50 - 0.08*Cs)/(TwoPi*0.42);
   return X - Sn*(Black1 - Black2*Cs)/(TwoPi*Black0);
}

double BlackmanSquare(double X) {
   return
      - 0.6595044010905501*(cos(X) - 1.0)
      + 0.1601741366715479*(cos(2.0*X) - 1.0)
      - 0.0010709178680006*(cos(4.0*X) - 1.0)
      + 0.0001450093579222*(cos(5.0*X) - 1.0)
      + 0.0001008528049040*(cos(7.0*X) - 1.0)
      + 0.0000653092892874*(cos(8.0*X) - 1.0)
      + 0.0000293385615146*(cos(10.0*X) - 1.0)
      + 0.0000205351559060*(cos(11.0*X) - 1.0)
      + 0.0000108567682890*(cos(13.0*X) - 1.0)
      + 0.0000081549460136*(cos(14.0*X) - 1.0)
      + 0.0000048519309366*(cos(16.0*X) - 1.0)
      + 0.0000038284344102*(cos(17.0*X) - 1.0)
      + 0.0000024753630724*(cos(19.0*X) - 1.0);
}

double *MakeITab(Integer Is) { // A Blackman Square table.
   double Omega = TwoPi/Is;
   double *ITab = malloc(++Is*sizeof *ITab); // Make room for an extra 3.0.
   for (Integer I = 0; I < Is; I++) ITab[I] = BlackmanSquare(I*Omega);
   return ITab;
}

// Downsample the signal Xi:Ni to Xo:No using a Blackman function.
void DownSample(double *Xi, double *Xo, Integer Ni, Integer No) {
// This assumes Ni ≤ No.
   double I2O = (double)Ni/No, O2I = 1.0/I2O; // The scaling ratio (> 1.0) and its inverse.
   for (Integer no = 0; no < No; no++) {
      Integer niL = RoundUp(I2O*(no - 1)); if (niL < 0) niL = 0;
      Integer niH = RoundDn(I2O*(no + 1)); if (niH >= Ni) niH = Ni - 1;
      double SumG = 0.0, SumXG = 0.0; // The weighting sums.
   // Convolve over the interval [niL, niH].
      for (Integer ni = niL; ni <= niH; ni++) {
         double G = Blackman(O2I*ni - no); SumG += G, SumXG += Xi[ni]*G;
      }
      Xo[no] = SumXG/SumG;
   }
}

// Upsample the signal Xi:Ni to Xo:No using a Blackman square function (i.e. a Blackman function convolved with a square).
// It's like smoothing the result of a nearest neighbor interpolation with a Blackman finite impulse response window.
void UpSample(double *Xi, double *Xo, Integer Ni, Integer No, double *ITab, Integer Is) {
// This assumes Ni ≥ No.
   double Is3 = Is/3.0, I2O = (double)Ni/No;
   for (Integer no = 0; no < No; no++) {
   // The center position in Xi and the corresponding range [niL, niH].
      double n = no*I2O;
      Integer niL = RoundTo(n - 1); if (niL < 0) niL = 0;
      Integer niH = RoundTo(n + 1); if (niH >= Ni) niH = Ni - 1;
      double SumX = 0.0;
      for (Integer ni = niL; ni <= niH; ni++) {
      // The position X ⇒ the look-up table index Xq and fractional offset Xr.
         double X = Is3*(ni - n + 1.5), Xr = fmod(X, 1.0); Integer Xq = (Integer)X;
      // Convolve with the linearly interpolated windowing function.
         SumX += Xi[ni]*((1.0 - Xr)*ITab[Xq] + Xr*ITab[Xq + 1]);
      }
      Xo[no] = SumX;
   }
}

// Convert the sound W:Ws; scaled with sampling frequency NuW
// to a scalogram image XY:*XsP:Ys, with *XsP pixel width, Ys bands; scaled with pixel frequency NuX and Bpo bands/octave.
double **Analyze(double *W, Integer Ws, Integer NuW, Integer *XsP, Integer Ys, double Bpo, double NuX, double LoP, double HiP) {
   double *PTab = MakePTab(LoP, HiP, Ys); // Tabulate the bands' central frequencies.
// Size up and make the image.
   Integer Xs = RoundUp(Ws*NuX);
   printf("Image size: %ld x %ld\n", (long)Xs, (long)Ys);
   double **XY = malloc(Ys*sizeof *XY);
// Zero-pad W and align it to a larger even size NW, for simplicity's sake, and DFT encode it.
// Note: padding cannot be done with circular convolution.
   Integer NW = Ws - 1;
   if (LogBase == 1.0) NW += RoundTo(5.0/PTab[1] - PTab[0]); // Linear mode.
   else {
      double K = exp(-log(LogBase)/Bpo);
      NW += RoundTo(10.0/(PTab[0]*K*(1.0 - K)));
   }
   if (NW%2 == 1) NW++; // Even-align NW for the sake of simplicity.
   NW = RoundTo((double)NextUnPrime(RoundTo(NW*NuX))/NuX);
// In-place DFT of the original signal resized from Ws to NW and zero-padded over the interval [Ws, NW).
   printf("(EnDFT %ld) %0.2f Hz. - %0.2f Hz.\r", NW, LoP, HiP);
   W = realloc(W, NW*sizeof *W), memset(&W[Ws], 0, (NW - Ws)*sizeof *W), EnDFT(W, NW);
   Integer LoGs = RoundTo(NW*NuX); // The (band-independent) length of the down-sampled envelopes.
   for (Integer Y = 0; Y < Ys; Y++) {
   // The index range of the band in the frequency domain, in absolute and relative form.
      Integer LoF = RoundTo(NW*EnNote((double)(Y - 1)/(Ys - 1), LoP, HiP));
      Integer HiF = RoundTo(NW*EnNote((double)(Y + 1)/(Ys - 1), LoP, HiP));
      double LoL = DeNote((double)LoF/NW, LoP, HiP);
      double HiL = DeNote((double)HiF/NW, LoP, HiP);
   // Restrict F to the integer range (0, Nyquist] = [1, NW/2].
      if (HiF > NW/2) HiF = NW/2; if (LoF < 1) LoF = 1;
   // Size up the filtered signal: +1 for DC, *2 for complex values, but no +1 for Nyquist since the signal is odd-length.
      Integer Gs = 2*(HiF - LoF) + 1;
   // Widen narrow bands and make the length more composite.
      if (Gs < LoGs) Gs = LoGs; else if (Gs > LoGs) Gs = NextUnPrime(Gs);
      printf("%ld/%ld (DeDFT %ld) %0.2f Hz. - %0.2f Hz.\r", Y + 1, Ys, Gs, (double)LoF*NuW/NW, (double)HiF*NuW/NW);
      double *Z = calloc(Gs, sizeof *Z);
      for (Integer F = LoF; F < HiF; F++) {
      // Apply a cosine (Hann) window over the relative scales L ∈ [LoL, HiL] corresponding to F ∈ [LoF, HiF].
         double L = DeNote((double)F/NW, LoP, HiP);
         double G = 0.5*(1.0 - cos(TwoPi*(L - LoL)/(HiL - LoL)));
         Integer dF = F - LoF + 1;
         Z[dF] = G*W[F + 1], Z[Gs - dF] = G*W[NW - 1 - F];
      }
   // In-place DFT decode the bands and combine to get the analytic signal's amplitude (the "envelope").
      DeDFT(Z, Gs, true);
   // The net result of this is Z(q) = average W_S(q, p): p ∈ [LoF, HiF].
   // Trim or downsample the result (using Blackman down-sampling).
      if (Gs < LoGs) Z = realloc(Z, LoGs*sizeof *Z);
      else if (Gs > LoGs) { // If the band *has* to be downsampled.
         double *ExZ = Z; Z = malloc(LoGs*sizeof *Z), DownSample(ExZ, Z, Gs, LoGs), free(ExZ);
      }
      XY[Ys - Y - 1] = realloc(Z, Xs*sizeof *Z);
   }
   putchar('\n');
   Normalize(XY, Xs, Ys, 1.0);
   *XsP = Xs; return XY;
}

double *MakeGTab(Integer Gs, double Bw) {
   double *GTab = malloc(Gs*sizeof *GTab); // The kernel.
   GTab[0] = 0.0; for (Integer G = 1; G < Gs; G++) GTab[G] = 1.0;
   double dG = Bw*(double)(Gs - 1); Integer LoGs = RoundUp(dG); // Double transition bandwidth in real and integer form.
   for (Integer G = 0; G < LoGs; G++) GTab[Gs - 1 - G] = GTab[G + 1] = BlackmanInt(G/dG);
   return GTab;
}

// Convert the scalogram image XY:Xs:Ys with Xs pixel width, Ys bands; scaled with pixel frequency NuX and Bpo bands/octave
// to the sound WZ:*WsP; scaled with sampling frequency NuW.
double *Synthesize(double **XY, Integer Xs, Integer Ys, Integer *WsP, Integer NuW, double LoP, double HiP, double NuX, double Bpo) {
   const double TransitionBW = 16.0; // The transition bandwidth for the low-pass filter on the synthesis envelopes.
   double *PTab = MakePTab(LoP, HiP, Ys); // Tabulate the bands' central frequencies.
// The up-sampled, frequency-shifted band envelope and its length (made ≥ 2 Xs for circular convolution).
   Integer Zs = NextUnPrime(2*Xs); double *Z = malloc(Zs*sizeof *Z);
   printf("Sound duration: %0.3f seconds.\n", Xs/(NuX*NuW));
// The central frequency: indexed as Z[MidZ] + i Z[Zs - MidZ].
   Integer MidZ = RoundTo(0.25*Zs);
// The envelope frequency-domain filter and the (complex) length of its DFT.
// The DC element is included, but not the Nyquist element.
   Integer Gs = (Zs + 1) >> 1; double *GTab = MakeGTab(Gs, 1.0/TransitionBW);
// The shifted band for the output and its length.
// Changing the length stretches the envelopes.
   Integer Ws = NextUnPrime(RoundTo(0.5*Zs/NuX)); double *W = calloc(Ws, sizeof *W);
// The (complex) length of the output's DFT.
// The DC element is included, but not the Nyquist element.
   Integer NW = (Ws + 1) >> 1;
// The frequency domain filter.
   for (Integer Y = 0; Y < Ys; Y++) {
   // Reset and shift Z with a sine table randomly phased ∈ [-Pi, +Pi] and DFT-encode the result.
      memset(Z, 0, Zs*sizeof *Z);
   // Envelope sampling rate*2 and frequency shifting with random phase by 0.25.
      double Omega = Pi*RandNum(), Cs = cos(Omega), Sn = sin(Omega); // The band's random phase, as a unit complex number Cs + i Sn.
      for (Integer X = 0; X < Xs; X++) {
         double G = X&1? -XY[Ys - Y - 1][X]: +XY[Ys - Y - 1][X];
         Z[2*X] = Cs*G, Z[2*X + 1] = Sn*G; // Was -Sn*G, but it doesn't really make any difference either way.
      }
      EnDFT(Z, Zs);
      Integer MidG = RoundTo(PTab[Y]*Ws); // The central (DC) index of the band's envelope in the frequency domain of W.
      printf("%ld/%ld (EnDFT %ld) %0.2f Hz.\r", (long)Y + 1, (long)Ys, (long)Zs, (double)MidG*NuW/Ws);
   // Envelope (by complex multiplication) all frequencies ∈ (0.0, 0.5) other than MidG.
      for (Integer G = 1; G < Gs; G++) {
         Integer S = G + MidG - MidZ;
      // Frequencies ∈ (0.0, 0.5) of G except at MidG.
         if (S > 0 && S < NW) W[S] += Z[G]*GTab[G], W[Ws - S] += Z[Zs - G]*GTab[G];
      }
   }
// Convert back to the time domain, then trim the result (by ignoring the tail) and normalize it.
   DeDFT(W, Ws, false), Ws = RoundTo(Xs/NuX), Normalize(&W, Ws, 1, 1.0);
   printf("DeDFT %ld %0.2f Hz. - %0.2f Hz.\n", (long)Ws, LoP, HiP);
   *WsP = Ws; return W;
}

double LoopTime = 10.0; // The duration in seconds of noise loops.
Integer ITabN = 16000; // The size of the Blackman Square look-up table: kept small enough to be cachable.

double *Ramble(double **XY, Integer Xs, Integer Ys, Integer *WsP, Integer NuW, double LoP, double HiP, double NuX, double Bpo) {
   double *PTab = MakePTab(LoP, HiP, Ys); // Tabulate the bands' central frequencies.
   Integer Ws = RoundTo(Xs/NuX); double *W = calloc(Ws, sizeof *W); // The sound signal and its length.
   printf("Sound duration: %0.3f seconds.\n", (double)Ws/NuW);
   double *GW = malloc(Ws*sizeof *GW); // The interpolated envelope.
// Estimate the sample size of the longest windowed sinc; or for non-linear modes: the time-domain size of the longest FIR.
   Integer LoFs;
   if (LogBase == 1.0) LoFs = RoundTo(20.0/PTab[1] - PTab[0]);
   else {
      double K = exp(-log(LogBase)/Bpo); // Was originally K = exp(-log(2.0)/Bpo) up to 0.2.3.
      LoFs = RoundTo(10.0/(PTab[0]*K*(1.0 - K)));
   }
// Set up the noise loops, sizing them up to a composite number to ease up the inverse DFTs.
// The size is deduced from LoopTime, which itself is eventually meant to be taken from user input.
   Integer Fs = LoopTime*NuW; if (Fs < LoFs) Fs = LoFs;
   Fs = NextUnPrime(Fs);
   double *Pink = calloc(Fs, sizeof *Pink);
   for (Integer L = 1; L < (Fs + 1) >> 1; L++) {
   // TODO: The expression for Mag needs to be verified and (if necessary) fixed.
      double Mag = pow((double)L, 0.5*(1.0 - LogBase)), Omega = Pi*RandNum();
   // The band's randomly-phased value, as a complex number Mag*(cos(Omega) + i sin(Omega)).
      Pink[L] = Mag*cos(Omega), Pink[Fs - L] = Mag*sin(Omega);
   }
   double *White = malloc(Fs*sizeof *White);
// The interpolation table (using a Blackman square window).
   double *ITab = MakeITab(ITabN);
// The frequency-domain filter.
   for (Integer Y = 0; Y < Ys; Y++) {
   // The index range of the band in the frequency domain, in absolute and relative form.
      Integer LoF = RoundTo(Fs*EnNote((double)(Y - 1)/(Ys - 1), LoP, HiP));
      Integer HiF = RoundTo(Fs*EnNote((double)(Y + 1)/(Ys - 1), LoP, HiP));
      double LoL = DeNote((double)LoF/Fs, LoP, HiP);
      double HiL = DeNote((double)HiF/Fs, LoP, HiP);
   // Restrict F to the integer range (0, Nyquist] = [1, Fs/2].
      if (HiF > Fs/2) HiF = Fs/2; if (LoF < 1) LoF = 1;
      printf("%ld/%ld (DeDFT %ld) %0.2f Hz. - %0.2f Hz.\r", (long)Y + 1, (long)Ys, (long)Fs, (double)LoF*NuW/Fs, (double)HiF*NuW/Fs);
   // Create the noise loop in the frequency domain and DFT decode it.
      memset(White, 0, Fs*sizeof *White);
      for (Integer F = LoF; F < HiF; F++) {
      // Apply a cosine window for relative scales L ∈ [LoL, HiL] corresponding to F ∈ [LoF, HiF].
         double L = DeNote((double)F/Fs, LoP, HiP);
         double G = 0.5*(1.0 - cos(TwoPi*(L - LoL)/(HiL - LoL)));
         White[F + 1] = G*Pink[F + 1], White[Fs - 1 - F] = G*Pink[Fs - 1 - F];
      }
      DeDFT(White, Fs, false);
   // Interpolate the envelope.
      UpSample(XY[Ys - Y - 1], GW, Xs, Ws, ITab, ITabN);
   // Modulate.
      for (Integer w = 0, Gr = 0; w < Ws; w++) { // The noise loop index: Gr ≡ w%Fs.
         W[w] += GW[w]*White[Gr]; if (++Gr >= Fs) Gr = 0;
      }
   }
   putchar('\n');
   Normalize(&W, Ws, 1, 1.0);
   *WsP = Ws; return W;
}

// Similar to gamma correction, but on a logarithmic scale; e.g. Shine ≡ 2 means square root.
void Brighten(double **XY, Integer Xs, Integer Ys, double Shine) {
   for (Integer Y = 0; Y < Ys; Y++) for (Integer X = 0; X < Xs; X++) XY[Y][X] = pow(XY[Y][X], Shine);
}

// The Image I/O Routines.
double DeColor(Color C) {
   long B = C->B, G = C->G, R = C->R;
   const double Scale = 1.0/0xff;
   if (B == 0xff) {
      if (G == 0x00) return 2.0 - R*Scale; else if (R == 0x00) return 2.0 + G*Scale; else if (G == 0xff && R == 0xff) return 7.0;
   } else if (B == 0x00) {
      if (G == 0xff) return 4.0 + R*Scale; else if (R == 0xff) return 6.0 - G*Scale; else if (G == 0x00 && R == 0x00) return 0.0;
   } else if (G == 0xff && R == 0x00) return 4.0 - B*Scale;
   else if (G == 0x00 && B == R) return B*Scale;
   else if (R == 0xff && B == G) return 6.0 + B*Scale;
   return 0.0;
}

void EnColor(Color C, double Amp) {
   double b, g, r;
   if (Amp < 0.0) b = g = r = 0.0;
   else if (Amp < 1.0) b = r = Amp, g = 0.0;
   else if (Amp < 2.0) b = 1.0, g = 0.0, r = 2.0 - Amp;
   else if (Amp < 3.0) b = 1.0, g = Amp - 2.0, r = 0.0;
   else if (Amp < 4.0) b = 4.0 - Amp, g = 1.0, r = 0.0;
   else if (Amp < 5.0) b = 0.0, g = 1.0, r = Amp - 4.0;
   else if (Amp < 6.0) b = 0.0, g = 6.0 - Amp, r = 1.0;
   else if (Amp < 7.0) b = g = Amp - 6.0, r = 1.0;
   else b = g = r = 1.0;
   long B = RoundTo(b*0xff); if (B < 0) B = 0; else if (B >= 0x100) B = 0xff; C->B = B;
   long G = RoundTo(g*0xff); if (G < 0) G = 0; else if (G >= 0x100) G = 0xff; C->G = G;
   long R = RoundTo(r*0xff); if (R < 0) R = 0; else if (R >= 0x100) R = 0xff; C->R = R;
}

double **GetImage(char *InFile, bool Colored, double PerDb, long *XsP, long *YsP) {
   int Ok = 0;
   long Xs, Ys; FILE *InF = GetFile(InFile, &Xs, &Ys); if (InF == NULL) goto Exit0;
   Color Row = malloc(Xs*sizeof *Row); if (Row == NULL) goto Exit1;
// Allocate the image.
   double **XY = calloc(Ys, sizeof *XY); if (XY == NULL) goto Exit2;
   XY[0] = malloc(Xs*Ys*sizeof *XY[0]); if (XY[0] == NULL) goto Exit3;
   for (long Y = 1; Y < Ys; Y++) XY[Y] = XY[0] + Y*Xs;
// If !Colored, colors are to be converted to grey levels normalized to [0, 1] by averaging the three channels.
   double Scale = Colored? 1.0/7.0: 1.0/0xff/3.0;
   double k = PerDb <= 0.0? 0.0: log(10.0)/PerDb;
   for (long Y =  Ys - 1; Y >= 0; Y--) { // Read it upside down.
      if (!GetRow(InF, Row, Xs)) goto Exit4;
      for (long X = 0; X < Xs; X++) {
         Color C = &Row[X];
         double A = Scale*(Colored? DeColor(C): C->B + C->G + C->R); if (PerDb > 0.0) A = exp(k*(A - 1.0));
         XY[Y][X] = A;
      }
   }
   Ok = 1;
Exit4:
   if (!Ok) free(XY[0]), XY[0] = 0;
Exit3:
   if (!Ok) free(XY), XY = NULL;
Exit2:
   free(Row);
Exit1:
   fclose(InF);
Exit0:
   if (Ok) {
      if (XsP != NULL) *XsP = Xs;
      if (YsP != NULL) *YsP = Ys;
   }
   return XY;
}

bool PutImage(char *ExFile, bool Colored, double PerDb, double **XY, long Xs, long Ys) {
   bool Ok = false;
   FILE *ExF = PutFile(ExFile, Xs, Ys); if (ExF == NULL) goto Exit0;
   Color Row = malloc(Xs*sizeof *Row); if (Row == NULL) goto Exit1;
   if (PerDb > 0.0) {
      double AdB = PerDb/log(10.0);
      for (long Y = 0; Y < Ys; Y++) {
         double *RowG = XY[Y];
         for (long X = 0; X < Xs; X++) RowG[X] = 1.0 + AdB*log(RowG[X]);
      }
   }
   double Scale = (double)(Colored? 7: 0xff);
   for (long Y = Ys - 1; Y >= 0; Y--) { // Write it upside down.
      double *RowG = XY[Y];
      for (long X = 0; X < Xs; X++) {
         double Br = Scale*RowG[X]; if (Br > Scale) Br = Scale; else if (Br < 0.0) Br = 0.0;
         Color C = &Row[X]; if (Colored) EnColor(C, Br); else C->R = C->G = C->B = (int)Br;
      }
      if (!PutRow(ExF, Row, Xs)) goto Exit2;
   }
   Ok = true;
Exit2:
   free(Row);
Exit1:
   fclose(ExF);
Exit0:
   return Ok;
}

Integer GetSampleBits(void) {
   Integer Bits = 0;
   while (Bits != 8 && Bits != 16 && Bits != 32) {
      printf("Bits per sample (8/16/32) [16]: "), Bits = (Integer)GetNumber(); if (Bits == 0) Bits = 16; // The default value.
   // Originally, up to 0.2.3, a check against -0x7fffffff was also included for C90 compatibility.
   }
   return Bits;
}

// The Main Application.
// char *Version = "0.2.3", *Date = "2008 May 29", *Author = "Michel Rouzic";
// char *Version = "0.2.3.1", *Date = "2012 Dec 27", *Author = "Darth Ninja"; // First rewrite: redone layout and code normalization for analysis.
// char *Version = "0.2.3.2", *Date = "2014 Feb 09", *Author = "Darth Ninja"; // Second rewrite: for critical review and correction of errors.
// char *Version = "0.2.3.3", *Date = "2015 Mar 26", *Author = "Darth Ninja"; // Third rewrite: for recoding and removing system/library dependencies.
// char *Version = "0.3", *Date = "2015 Aug 06", *Author = "Darth Ninja"; // Self-contained (i.e. FFTW was removed) and thus portable.
// char *Version = "0.3.1", *Date = "2015 Nov 11", *Author = "Darth Ninja"; // Colored graphs.
// const char *Version = "0.3.2", *Date = "2015 Dec 13", *Author = "Darth Ninja"; // Decibel scale.
const char *Version = "0.3.2.1", *Date = "2020 Oct 12", *Author = "Lydia Marie Williamson"; // Re-integration with earlier versions and preparation inclusion in ALGLIB++.

bool Quiet = false;

void Configure(Integer *YsP, Integer Ws, Integer *NuWP, double *LoPP, double *HiPP, double *NuXP, double *PerDbP, double *BpoP, Integer Xs, Integer Synth) {
// printf("Configure...\n");
// Synth: 0 = Analysis, 1 = Synthesis.
#ifdef WIN32
   const char *CfgFile = "arss.conf";
#else
   char CfgFile[FILENAME_MAX]; sprintf(CfgFile, "%s/%s", getenv("HOME"), ".arss.conf"); // Path to the configuration file.
#endif
   FILE *CfgF = fopen(CfgFile, "rb"); // Open the configuration file.
   Integer NuW = *NuWP;
   if (NuW == 0) { // In synthesis mode, with no NuW yet defined.
      if (Quiet) CheckOut(true, "Define the output sample rate, using --sample-rate (-r).");
      printf("Sample rate [44100]: "), NuW = GetNumber(); if (NuW == 0) NuW = 44100; // The default value.
   // Originally, up to 0.2.3, a check against -0x7fffffff was also included for C90 compatibility.
   }
   double LoP = *LoPP, HiP = *HiPP, NuX = *NuXP, PerDb = *PerDbP, Bpo = *BpoP; Integer Ys = *YsP;
   bool LoPSet = LoP != 0.0, HiPSet = HiP != 0.0, NuXSet = NuX != 0.0, PerDbSet = PerDb >= 0.0, BpoSet = Bpo != 0.0, XsSet = Xs != 0, YsSet = Ys != 0;
   if (LoPSet && HiPSet && BpoSet && YsSet)
      CheckOut(true, "--min-freq (-min), --max-freq (-max), --bpo (-b)%s cannot all be used together.", Synth? " and --height (-y)": "");
   if (NuXSet && XsSet && !Synth) CheckOut(true, "--width (-x) and --pps (-p) cannot be used together.");
   size_t CfgN = 0; // Used to size-check the configuration file.
   if (CfgF != NULL) { // Load whatever settings from the file are present.
      for (Integer I = 0; I < 4*sizeof(double); I++) { char Ch; CfgN += fread(&Ch, sizeof Ch, 1, CfgF); }
      rewind(CfgF);
   }
   if (CfgN > 0) { // Read in whatever has not yet been set.
      double D;
      if (LoP == 0.0) fread(&LoP, sizeof LoP, 1, CfgF); else fread(&D, sizeof D, 1, CfgF);
      if (HiP == 0.0) fread(&HiP, sizeof HiP, 1, CfgF); else fread(&D, sizeof D, 1, CfgF);
      if (Bpo == 0.0) fread(&Bpo, sizeof Bpo, 1, CfgF); else fread(&D, sizeof D, 1, CfgF);
      if (NuX == 0.0) fread(&NuX, sizeof NuX, 1, CfgF); else fread(&D, sizeof D, 1, CfgF);
   } else { // Otherwise, load in the default values.
      if (LoP == 0.0) LoP = 27.5;
      if (HiP == 0.0) HiP = 22050.0; // Was originally 20000, up to 0.2.3.
      if (Bpo == 0.0) Bpo = 12.0;
      if (NuX == 0.0) NuX = 150.0;
   }
   if (CfgF != NULL) fclose(CfgF), CfgF = NULL;
   if (!LoPSet && !(HiPSet && BpoSet && YsSet)) {
      if (Quiet) CheckOut(true, "Define a minimum frequency, using --min-freq (-min).");
      printf("Minimum frequency (Hz.) [%0.3f]: ", LoP); double D = GetNumber(); if (D != 0.0) LoP = D;
      LoPSet = true;
   }
   if (!PerDbSet) {
      PerDb = 0.0;
      if (Quiet) CheckOut(true, "Set the number of dB sensitivity, using --per-dB (-d).");
      printf("dB sensitivity [%0.3f]: ", PerDb); double D = GetNumber(); if (D != 0.0) PerDb = D;
      PerDbSet = true;
   }
   if (!BpoSet && !(LoPSet && HiPSet && YsSet)) {
      if (Quiet) CheckOut(true, "Set the number of bands per octave, using --bpo (-b).");
      printf("Bands per octave [%0.3f]: ", Bpo); double D = GetNumber(); if (D != 0.0) Bpo = D;
      BpoSet = true;
   }
   if (!HiPSet && !(LoPSet && BpoSet && YsSet)) {
      Integer I = 0; double NuNy = 0.5*NuW, F;
      do I++, F = LoP*pow(LogBase, I/Bpo); while (F < NuNy);
      double MaxNu = LoP*pow(LogBase, (I - 2)/Bpo); // Maximum allowed frequency.
      if (HiP > MaxNu) HiP = MaxNu - fmod(MaxNu, 1.0); // Replaces the "Upper frequency limit above Nyquist frequency" warning.
      if (!Synth) {
         if (Quiet) CheckOut(true, "Define a maximum frequency, using --max-freq (-max).");
         printf("Maximum frequency (Hz.) (up to %0.3f) [%0.3f]: ", MaxNu, HiP); double D = GetNumber(); if (D != 0.0) HiP = D;
         if (HiP > MaxNu) HiP = MaxNu - fmod(MaxNu, 1.0); // Replaces the "Upper frequency limit above Nyquist frequency" warning.
      }
      HiPSet = true;
   }
// The frequency bounds in Hz.
   if (!LoPSet) LoP = HiP*pow(LogBase, -(Ys - 1)/Bpo), printf("Minimum frequency: %0.3f Hz.\n", LoP);
   if (!HiPSet) HiP = LoP*pow(LogBase, (Ys - 1)/Bpo), printf("Maximum frequency: %0.3f Hz.\n", HiP);
// The sound and image limits.
   if (!YsSet) Ys = 1 + RoundTo(Bpo*(LogB(HiP) - LogB(LoP))), printf("Bands: %d\n", Ys);
   if (!BpoSet) Bpo = LogBase == 1.0? HiP/NuW: (Ys - 1)/(LogB(HiP) - LogB(LoP)), printf("Bands per octave: %0.3f\n", Bpo);
// Calculate NuX, if we're in Analysis mode with Xs set (by the user).
   if (XsSet && !Synth) NuX = (double)Xs*NuW/Ws, printf("Pixels per second: %0.3f\n", NuX);
// If in Analysis mode with none set or NuX isn't set in Synthesis mode
   if (!NuXSet && !(!Synth && XsSet)) {
      if (Quiet) CheckOut(true, "Define a pixels per second setting, using --pps (-p).");
      printf("Pixels per second [%0.3f]: ", NuX); double D = GetNumber(); if (D != 0.0) NuX = D;
   }
// Saving the settings to a file.
   CfgF = fopen(CfgFile, "wb"); if (CfgF == NULL) CheckOut(true, "Cannot write the configuration to the file %s.", CfgFile);
   fwrite(&LoP, sizeof LoP, 1, CfgF), fwrite(&HiP, sizeof HiP, 1, CfgF);
   fwrite(&Bpo, sizeof Bpo, 1, CfgF), fwrite(&NuX, sizeof NuX, 1, CfgF);
   fclose(CfgF);
   *YsP = Ys, *NuWP = NuW, *PerDbP = PerDb, *BpoP = Bpo;
// Make LoP and NuX relative to the sampling rate, instead of in Hz..
   *LoPP = LoP /= NuW, *NuXP = NuX /= NuW;
// HiP is the central frequency of the last band.
// For linear mode, it is obtained from Bpo.
   *HiPP = LogBase == 1.0? Bpo: LoP*pow(LogBase, (double)(Ys - 1)/Bpo);
}

void ShowHelp(char *App) {
   printf("Usage: %s [options] input_file output_file [options].\n", App);
   printf("Example: %s -q in.bmp out.wav --noise --min-freq 55 -max 16000 --pps 100 -r 44100 -f 16\n", App);
   printf(
      "Basic Options:\n"
      "--help,         -h, /?       Display these options.\n"
   // --adv-help was originally an option separately displaying advaned help options.
      "--version,      -v           Display the program version.\n"
      "--quiet,        -q           No-prompt mode. Useful for scripting.\n"
      "--analysis,     -a           Analysis mode. Not needed if input file is *.wav.\n"
      "--sine,         -s           Sound synthesis.\n"
      "--noise,        -n           Noise synthesis.\n"
      "--min-freq,     -min [real]  Minimum frequency in Hertz.\n"
      "--max-freq,     -max [real]  Maximum frequency in Hertz.\n"
      "--per-dB        -d [real]    dB sensitivity (0 means: linear amplitude scale).\n"
      "--bpo,          -b [real]    Bands per octave frequency resolution.\n"
      "--pps,          -p [real]    Pixels per second time resolution.\n"
      "--height,       -y [integer] Graph height.\n"
      "--width,        -x [integer] Graph width.\n"
      "--sample-rate,  -r [integer] Output sample rate.\n"
      "--color,        -c           Color graphical display.\n"
      "--brightness,   -g [real]    Brightness inverse exponent: f(x) = x^{1/Brightness}.\n"
   // "--time,         -t [real]    Sound duration in seconds.\n"
      "--format-param, -f [integer] Bits per sound sample (8/16/32, default: 32).\n"
      "Advanced Options:\n"
      "--log-base         [real]    Frequency scale per 'octave' (default: 2).\n"
      "--linear,       -l           Linear frequency scale. Same as --log-base 1.\n"
      "--loop-size        [real]    Noise table size in seconds (default: 10).\n"
      "--bmsq-lut-size    [integer] Interpolation table size (default: 16000).\n"
   // --Pi [real] was originally an option permitting a value for π tp be set (default: 3.1415926535897932).
   );
}

int main(int AC, char *AV[]) {
   char *App = AC < 1? NULL: AV[0]; if (App == NULL || *App == '\0') App = "ARSS";
   int Status = EXIT_FAILURE;
   printf("The Analysis / Resynthesis Sound Scalograph %s.\n", Version);
   srand(time(NULL));
   Integer NuW = 0, Bits = 0, Xs = 0, Ys = 0;
   double LoP = 0.0, HiP = 0.0, PerDb = -1.0, NuX = 0.0, Bpo = 0.0, Shine = 1.0;
   bool Colored = false;
   char Mode = 0, *InFile = NULL, *ExFile = NULL;
   char *EndP = "";
   for (int A = 1; A < AC; ) {
      char *Arg = AV[A++], *Par;
      if (strcmp(Arg, "/?") == 0) Arg[0] = '-', Arg[1] = 'h'; // DOS-friendly help.
      if (Arg[0] != '-') // If the argument is not a parameter.
         if (InFile == NULL) InFile = Arg; else if (ExFile == NULL) ExFile = Arg;
         else { fprintf(stderr, "Only two files may be specified: \"%s\" was also listed.\n", Arg); goto UnDo0; }
      else if (strcmp(Arg, "--analysis") == 0 || strcmp(Arg, "-a") == 0) Mode = 1;
      else if (strcmp(Arg, "--sine") == 0 || strcmp(Arg, "-s") == 0) Mode = 2;
      else if (strcmp(Arg, "--noise") == 0 || strcmp(Arg, "-n") == 0) Mode = 3;
      else if (strcmp(Arg, "--color") == 0 || strcmp(Arg, "-c") == 0) Colored = true;
      else if (strcmp(Arg, "--quiet") == 0 || strcmp(Arg, "-q") == 0) Quiet = true;
      else if (strcmp(Arg, "--linear") == 0 || strcmp(Arg, "-l") == 0) LogBase = 1.0;
      else if (strcmp(Arg, "--sample-rate") == 0 || strcmp(Arg, "-r") == 0) {
         if (A >= AC) { fprintf(stderr, "Missing option for %s.\n", Arg); goto UnDo0; } else Par = AV[A++];
         NuW = strtol(Par, &EndP, 10); if (*EndP != '\0') { fprintf(stderr, "Invalid %s value (%s).\n", Arg, Par); goto UnDo0; }
      } else if (strcmp(Arg, "--min-freq") == 0 || strcmp(Arg, "-min") == 0) {
         if (A >= AC) { fprintf(stderr, "Missing option for %s.\n", Arg); goto UnDo0; } else Par = AV[A++];
         LoP = strtod(Par, &EndP); if (*EndP != '\0') { fprintf(stderr, "Invalid %s value (%s).\n", Arg, Par); goto UnDo0; }
         if (LoP == 0.0) LoP = DBL_MIN; // Make it infinitesimal so that it is treated as having been set (and to avoid DC on a logarithmic scale).
      } else if (strcmp(Arg, "--max-freq") == 0 || strcmp(Arg, "-max") == 0) {
         if (A >= AC) { fprintf(stderr, "Missing option for %s.\n", Arg); goto UnDo0; } else Par = AV[A++];
         HiP = strtod(Par, &EndP); if (*EndP != '\0') { fprintf(stderr, "Invalid %s value (%s).\n", Arg, Par); goto UnDo0; }
      } else if (strcmp(Arg, "--per-dB") == 0 || strcmp(Arg, "-d") == 0) {
         if (A >= AC) { fprintf(stderr, "Missing option for %s.\n", Arg); goto UnDo0; } else Par = AV[A++];
         PerDb = strtod(Par, &EndP); if (*EndP != '\0') { fprintf(stderr, "Invalid %s value (%s).\n", Arg, Par); goto UnDo0; }
      } else if (strcmp(Arg, "--bpo") == 0 || strcmp(Arg, "-b") == 0) {
         if (A >= AC) { fprintf(stderr, "Missing option for %s.\n", Arg); goto UnDo0; } else Par = AV[A++];
         Bpo = strtod(Par, &EndP); if (*EndP != '\0') { fprintf(stderr, "Invalid %s value (%s).\n", Arg, Par); goto UnDo0; }
      } else if (strcmp(Arg, "--pps") == 0 || strcmp(Arg, "-p") == 0) {
         if (A >= AC) { fprintf(stderr, "Missing option for %s.\n", Arg); goto UnDo0; } else Par = AV[A++];
         NuX = strtod(Par, &EndP); if (*EndP != '\0') { fprintf(stderr, "Invalid %s value (%s).\n", Arg, Par); goto UnDo0; }
      } else if (strcmp(Arg, "--height") == 0 || strcmp(Arg, "-y") == 0) {
         if (A >= AC) { fprintf(stderr, "Missing option for %s.\n", Arg); goto UnDo0; } else Par = AV[A++];
         Ys = strtol(Par, &EndP, 10); if (*EndP != '\0') { fprintf(stderr, "Invalid %s value (%s).\n", Arg, Par); goto UnDo0; }
      } else if (strcmp(Arg, "--width") == 0 || strcmp(Arg, "-x") == 0) {
         if (A >= AC) { fprintf(stderr, "Missing option for %s.\n", Arg); goto UnDo0; } else Par = AV[A++];
         Xs = strtol(Par, &EndP, 10); if (*EndP != '\0') { fprintf(stderr, "Invalid %s value (%s).\n", Arg, Par); goto UnDo0; }
      } else if (strcmp(Arg, "--loop-size") == 0) {
         if (A >= AC) { fprintf(stderr, "Missing option for %s.\n", Arg); goto UnDo0; } else Par = AV[A++];
         LoopTime = strtod(Par, &EndP); if (*EndP != '\0') { fprintf(stderr, "Invalid %s value (%s).\n", Arg, Par); goto UnDo0; }
      } else if (strcmp(Arg, "--log-base") == 0) {
         if (A >= AC) { fprintf(stderr, "Missing option for %s.\n", Arg); goto UnDo0; } else Par = AV[A++];
         LogBase = strtod(Par, &EndP); if (*EndP != '\0') { fprintf(stderr, "Invalid %s value (%s).\n", Arg, Par); goto UnDo0; }
      } else if (strcmp(Arg, "--bmsq-lut-size") == 0) {
         if (A >= AC) { fprintf(stderr, "Missing option for %s.\n", Arg); goto UnDo0; } else Par = AV[A++];
         ITabN = strtol(Par, &EndP, 10); if (*EndP != '\0') { fprintf(stderr, "Invalid %s value (%s).\n", Arg, Par); goto UnDo0; }
   // } else if (strcmp(Arg, "--pi") == 0) {
   //    if (A >= AC) { fprintf(stderr, "Missing option for %s.\n", Arg); goto UnDo0; } else Par = AV[A++];
   //    Pi = strtod(Par, &EndP); if (*EndP != '\0') { fprintf(stderr, "Invalid %s value (%s).\n", Arg, Par); goto UnDo0; }
   //    TwoPi = 2.0*Pi;
      } else if (strcmp(Arg, "--format-param") == 0 || strcmp(Arg, "-f") == 0) {
         if (A >= AC) { fprintf(stderr, "Missing option for %s.\n", Arg); goto UnDo0; } else Par = AV[A++];
         Bits = strtol(Par, &EndP, 10); if (*EndP != '\0') { fprintf(stderr, "Invalid %s value (%s).\n", Arg, Par); goto UnDo0; }
      } else if (strcmp(Arg, "--brightness") == 0 || strcmp(Arg, "-g") == 0) {
         if (A >= AC) { fprintf(stderr, "Missing option for %s.\n", Arg); goto UnDo0; } else Par = AV[A++];
         Shine = strtod(Par, &EndP); if (*EndP != '\0') { fprintf(stderr, "Invalid %s value (%s).\n", Arg, Par); goto UnDo0; }
   // } else if (strcmp(Arg, "--time") == 0 || strcmp(Arg, "-t") == 0) {
   //    if (A >= AC) { fprintf(stderr, "Missing option for %s.\n", Arg); goto UnDo0; } else Par = AV[A++];
   //    Duration = strtod(Par, &EndP); if (*EndP != '\0') { fprintf(stderr, "Invalid %s value (%s).\n", Arg, Par); goto UnDo0; }
      } else if (strcmp(Arg, "--version") == 0 || strcmp(Arg, "-v") == 0) {
         printf("Copyright (C) 2005-2008 Michel Rouzic\nProgram last modified by %s on %s.\n", Author, Date);
         Status = EXIT_SUCCESS; goto UnDo0;
      } else if (strcmp(Arg, "--help") == 0 || strcmp(Arg, "-h") == 0) {
         ShowHelp(App); Status = EXIT_SUCCESS; goto UnDo0;
      } else { fprintf(stderr, "%s: unrecognized option.\n", Arg); goto UnDo0; }
   }
   if (InFile != NULL) printf("Input file: %s.\n", InFile);
   else {
      if (Quiet) CheckOut(true, "Specify an input file.");
      printf("Input file: "), InFile = GetString();
   }
   if (ExFile != NULL) printf("Output file: %s.\n", ExFile);
   else {
      if (Quiet) CheckOut(true, "Specify an output file.");
      printf("Output file: "), ExFile = GetString();
   }
// The Synthesis mode is assumed if the input filename is of the form *.wav, otherwise the Analysis mode is assumed.
   if (Mode == 0 && strstr(InFile, ".wav") != NULL && strstr(InFile, ".wav")[4] == '\0') Mode = 1;
   if (Mode < 1 || Mode > 3) {
      if (Quiet) CheckOut(true, "Specify an operation mode, using --analysis (-a), --sine (-s) or --noise (-n).");
      else do
         printf("Choose the mode, 1 = Analysis, 2 = Sine synthesis, 3 = Noise synthesis: "), Mode = GetNumber();
      while (Mode < 1 || Mode > 3);
   }
   Integer StartTime;
   if (Mode == 1) {
      Integer Ws = 0, Zs; double **WZ = GetSound(InFile, &Ws, &Zs, &NuW);
      if (WZ == NULL) CheckOut(true, "The sound could not be read from the file %s.\n", InFile);
//    printf("%ld samples, %ld channels\n", (long)Ws, (long)Zs);
      Configure(&Ys, Ws, &NuW, &LoP, &HiP, &NuX, &PerDb, &Bpo, Xs, 0);
      StartTime = Now();
      double **XY = Analyze(WZ[0], Ws, NuW, &Xs, Ys, Bpo, NuX, LoP, HiP);
      if (Shine != 1.0) Brighten(XY, Xs, Ys, 1.0/Shine);
      if (!PutImage(ExFile, Colored, PerDb, XY, Xs, Ys)) CheckOut(true, "The image could not be written to the file %s.\n", ExFile);
   } else if (Mode == 2 || Mode == 3) {
      double **XY = GetImage(InFile, Colored, PerDb, &Xs, &Ys);
      if (XY == NULL) CheckOut(true, "The image could not be read from the file %s.\n", InFile);
      if (Bits == 0) Bits = Quiet? 32: GetSampleBits(); // If Bits is not yet defined.
      Integer Ws = 0; Configure(&Ys, Ws, &NuW, &LoP, &HiP, &NuX, &PerDb, &Bpo, Xs, 1);
      StartTime = Now();
      if (Shine != 1.0) Brighten(XY, Xs, Ys, Shine);
      double *W = Mode == 2?
         Synthesize(XY, Xs, Ys, &Ws, NuW, LoP, HiP, NuX, Bpo): // Sinusoidal output.
         Ramble(XY, Xs, Ys, &Ws, NuW, LoP, HiP, NuX, Bpo); // Noised output.
      if (!PutSound(ExFile, &W, Ws, 1, NuW, Bits)) CheckOut(true, "The sound could not be written to the file %s.\n", ExFile);
   }
   CheckOut(false, "Processing time: %0.3f seconds.", (double)(Now() - StartTime)/1000.0);
   Status = EXIT_SUCCESS;
UnDo0:
   return Status;
}
