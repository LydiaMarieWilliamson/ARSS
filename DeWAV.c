#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "DeWAV.h"

typedef uint8_t Num1; // typedef unsigned short Num2;
typedef uint16_t Num2; // typedef unsigned short Num2;
typedef uint32_t Num4; // typedef unsigned long Num4;

static inline Num1 GetNum1(FILE *InF) { // Read an 8-bit integer.
   Num1 A; fread(&A, 1, 1, InF);
   return A;
}

static inline Num2 GetNum2(FILE *InF) { // Read a 16-bit integer in little endian.
   Num1 A[2]; fread(A, 2, 1, InF);
   return (Num2)((A[1] << 8) | A[0]);
}

static inline Num4 GetNum4(FILE *InF) { // Read a 32-bit integer in little endian.
   Num1 A[4]; fread(A, 4, 1, InF);
   return (Num4)((A[3] << 24) | (A[2] << 16) | (A[1] << 8) | A[0]);
}

static void GetWZ1(FILE *InF, double **WZ, long Ws, long Zs) {
   const double Sc = 1.0/0x80;
   for (long w = 0; w < Ws; w++) for (long z = 0; z < Zs; z++) {
      double D = GetNum1(InF); if (D >= 0x80) D -= 0x100;
      WZ[z][w] = Sc*D;
   }
}

static void GetWZ2(FILE *InF, double **WZ, long Ws, long Zs) {
   const double Sc = 1.0/0x8000;
   for (long w = 0; w < Ws; w++) for (long z = 0; z < Zs; z++) {
      double D = GetNum2(InF); if (D >= 0x8000) D -= 0x10000;
      WZ[z][w] = Sc*D;
   }
}

static void GetWZ4(FILE *InF, double **WZ, long Ws, long Zs) {
   const double Sc = 1.0/0x80000000;
   for (long w = 0; w < Ws; w++) for (long z = 0; z < Zs; z++) {
      double D = GetNum4(InF); if (D >= 0x80000000) D -= 0x100000000;
      WZ[z][w] = Sc*D;
   }
}

static void CheckOut(int Fatal, char *Format, ...) {
   va_list AP; va_start(AP, Format); vfprintf(stderr, Format, AP); va_end(AP);
   fputc('\n', stderr);
   if (Fatal > 1) fprintf(stderr, "Exiting with error.\n");
   if (Fatal) exit(EXIT_FAILURE);
}

double **GetSound(char *Path, long *WsP, long *ZsP, long *NuWP) {
   FILE *InF = fopen(Path, "rb"); if (InF == NULL) goto UnDo0;
// File Header.
   Num4 Riff = GetNum4(InF); if (Riff != ('R' | 'I' << 8 | 'F' << 16 | 'F' << 24)) CheckOut(1, "This file is not in RIFF format\n");
   Num4 ChunkSize = GetNum4(InF);
   Num4 Wave = GetNum4(InF); if (Wave != ('W' | 'A' << 8 | 'V' << 16 | 'E' << 24)) CheckOut(1, "This is a RIFF file that is not in WAVE format.\n");
// Data Header.
   Num4 Fmt = GetNum4(InF); if (Fmt != ('f' | 'm' << 8 | 't' << 16 | ' ' << 24)) CheckOut(1, "The format used by this WAVE file is not recognized.\n");
   Num4 FmtN = GetNum4(InF); if (FmtN != 16) CheckOut(1, "This format segment in this WAVE file is not recognized.\n");
   Num2 Mode = GetNum2(InF), Zs = GetNum2(InF);
   Num4 NuW = GetNum4(InF), ByteRate = GetNum4(InF);
   Num2 SampleSize = GetNum2(InF), Bits = GetNum2(InF);
// Wave Header.
   Num4 Data = GetNum4(InF); if (Data != ('d' | 'a' << 8 | 't' << 16 | 'a' << 24)) CheckOut(1, "This data segment in this WAVE file is not recognized.\n");
   Num4 Ws = GetNum4(InF)/(Bits/8)/Zs;
// Wave Data.
   double **WZ = malloc(Zs*sizeof *WZ); if (WZ == NULL) goto UnDo1;
   if ((WZ[0] = malloc(Ws*Zs*sizeof *WZ[0])) == NULL) goto UnDo2;
   for (long z = 1; z < Zs; z++) WZ[z] = WZ[0] + Ws*z;
   switch (Bits) {
      case 8: GetWZ1(InF, WZ, Ws, Zs); break;
      case 16: GetWZ2(InF, WZ, Ws, Zs); break;
      case 32: GetWZ4(InF, WZ, Ws, Zs); break;
   }
   fclose(InF);
   *WsP = Ws, *ZsP = Zs, *NuWP = NuW; return WZ;
UnDo2:
   free(WZ);
UnDo1:
   fclose(InF);
UnDo0:
   return NULL;
}
