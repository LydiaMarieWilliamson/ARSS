#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include "EnBMP.h"

typedef uint16_t Num2; // typedef unsigned short Num2;
typedef uint32_t Num4; // typedef unsigned long Num4;

// Generic BMP encoder.
static inline bool PutNum2(FILE *ExF, Num2 N) {
   int A = N&0xff, B = (N >> 8)&0xff;
   return fputc(A, ExF) != EOF && fputc(B, ExF) != EOF;
}

static inline bool PutNum4(FILE *ExF, Num4 N) {
   int A = N&0xff, B = (N >> 8)&0xff, C = (N >> 16)&0xff, D = (N >> 24)&0xff;
   return
      fputc(A, ExF) != EOF && fputc(B, ExF) != EOF &&
      fputc(C, ExF) != EOF && fputc(D, ExF) != EOF;
}

static bool PutHeader(FILE *ExF, long Xs, long Ys) {
   Num4 RowN = (3*Xs + 3)/4*4, Size = RowN*Ys, Offset = 54, FSize = Offset + Size, XX = 0;
   Num4 ILen = 40, Mode = 0, ResX = 0/*2834*/, ResY = 0/*2834*/, Colors = 0, MinColors = 0;
   Num2 Zs = 1, Bits = 24;
   return
// Header.
      PutNum2(ExF, ('M' << 8) | 'B') && PutNum4(ExF, FSize) &&
      PutNum4(ExF, XX) && PutNum4(ExF, Offset) &&
// Specs.
      PutNum4(ExF, ILen) && PutNum4(ExF, Xs) && PutNum4(ExF, Ys) &&
      PutNum2(ExF, Zs) && PutNum2(ExF, Bits) && PutNum4(ExF, Mode) && PutNum4(ExF, Size) &&
      PutNum4(ExF, ResX) && PutNum4(ExF, ResY) && PutNum4(ExF, Colors) && PutNum4(ExF, MinColors);
}

FILE *PutFile(char *Path, long Xs, long Ys) {
   FILE *ExF = Path == NULL? stdout: fopen(Path, "wb");
   if (ExF == NULL) fprintf(stderr, "Cannot put %s.\n", Path);
   else if (!PutHeader(ExF, Xs, Ys)) {
      if (Path != NULL) fclose(ExF);
      ExF = NULL;
   }
   return ExF;
}

bool PutRow(FILE *ExF, Color Row, long Xs) {
   for (long X = 0; X < Xs; X++) {
      Color C = &Row[X];
      if (fputc(C->B, ExF) == EOF) return false;
      if (fputc(C->G, ExF) == EOF) return false;
      if (fputc(C->R, ExF) == EOF) return false;
   }
   long RowN = (3*Xs + 3)/4*4;
   for (long I = 3*Xs; I < RowN; I++) if (fputc('\0', ExF) == EOF) return false;
   return true;
}
