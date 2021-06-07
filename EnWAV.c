#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "EnWAV.h"

typedef uint8_t Num1; // typedef unsigned short Num2;
typedef uint16_t Num2; // typedef unsigned short Num2;
typedef uint32_t Num4; // typedef unsigned long Num4;

static inline void PutNum1(FILE *ExF, Num1 N) { // Write a 8-bit integer.
   Num1 A = (Num1)N; fwrite(&A, sizeof A, 1, ExF);
}

static inline void PutNum2(FILE *ExF, Num2 N) { // Write a 16-bit integer in little endian.
   Num1 A[2] = { N&0xff, (N >> 8)&0xff };
   fwrite(A, 2, 1, ExF);
}

static inline void PutNum4(FILE *ExF, Num4 N) { // Write a 32-bit integer in little endian.
   Num1 A[4] = { N&0xff, (N >> 8)&0xff, (N >> 16)&0xff, (N >> 24)&0xff };
   fwrite(A, 4, 1, ExF);
}

static void PutWZ1(FILE *ExF, double **WZ, long Ws, long Zs) {
   for (long w = 0; w < Ws; w++) for (long z = 0; z < Zs; z++) {
      long N = lround(0x80*WZ[z][w]);
      if (N >= 0x80) N = 0x7f; else if (N < -0x80) N = 0x80; else if (N < 0) N += 0x100;
      PutNum1(ExF, N);
   }
}

static void PutWZ2(FILE *ExF, double **WZ, long Ws, long Zs) {
   for (long w = 0; w < Ws; w++) for (long z = 0; z < Zs; z++) {
      long N = lround(0x8000*WZ[z][w]);
      if (N >= 0x8000) N = 0x7fff; else if (N < -0x8000) N = 0x8000; else if (N < 0) N += 0x10000;
      PutNum2(ExF, N);
   }
}

static void PutWZ4(FILE *ExF, double **WZ, long Ws, long Zs) {
   for (long w = 0; w < Ws; w++) for (long z = 0; z < Zs; z++) {
      long N = lround(0x80000000*WZ[z][w]);
      if (N >= 0x80000000) N = 0x7fffffff; else if (N < -0x80000000) N = 0x80000000; else if (N < 0) N += 0x100000000;
      PutNum4(ExF, N);
   }
}

bool PutSound(char *Path, double **WZ, long Ws, long Zs, long Nu, long Bits) {
   bool Status = false;
   FILE *ExF = fopen(Path, "wb"); if (ExF == NULL) goto UnDo0;
   const Num4 Riff = 'R' | 'I' << 8 | 'F' << 16 | 'F' << 24;
   const Num4 Wave = 'W' | 'A' << 8 | 'V' << 16 | 'E' << 24;
   const Num4 Fmt_ = 'f' | 'm' << 8 | 't' << 16 | ' ' << 24;
   const Num4 Data = 'd' | 'a' << 8 | 't' << 16 | 'a' << 24;
   Num4 Bytes = Bits/8, Size = Ws*Zs*Bytes, ChunkSize = Size + 36;
// File Header: "RIFF" (Size) "WAVE"
   PutNum4(ExF, Riff), PutNum4(ExF, ChunkSize), PutNum4(ExF, Wave);
// Data Header: "fmt " (Size: 16) (Bits: 3 or 1) Zs Nu Nu_Byte Size Nits
   PutNum4(ExF, Fmt_), PutNum4(ExF, 16);
   PutNum2(ExF, Bits == 32? 3: 1), PutNum2(ExF, Zs);
   PutNum4(ExF, Nu), PutNum4(ExF, Nu*Bytes);
   PutNum2(ExF, Bytes), PutNum2(ExF, Bits);
// Wave Header: "data" (Size) [B* if Bits == 8, W* if Bits == 16, L* if Bits == 32]
   PutNum4(ExF, Data), PutNum4(ExF, Size);
// Wave Data.
   switch (Bits) {
      case 8: PutWZ1(ExF, WZ, Ws, Zs); break;
      case 16: PutWZ2(ExF, WZ, Ws, Zs); break;
      case 32: PutWZ4(ExF, WZ, Ws, Zs); break;
   }
   fclose(ExF);
   Status = true;
UnDo0:
   return Status;
}
