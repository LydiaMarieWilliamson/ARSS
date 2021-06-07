#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include "DeBMP.h"

typedef uint16_t Num2; // typedef unsigned short Num2;
typedef uint32_t Num4; // typedef unsigned long Num4;

// Generic BMP decoder.
static bool GetNum2(FILE *InF, Num2 *NP) {
   int A = fgetc(InF); if (A == EOF) return false;
   int B = fgetc(InF); if (B == EOF) return false;
   if (NP != NULL) *NP = (B << 8) | A;
   return true;
}

static bool GetNum4(FILE *InF, Num4 *NP) {
   int A = fgetc(InF); if (A == EOF) return false;
   int B = fgetc(InF); if (B == EOF) return false;
   int C = fgetc(InF); if (C == EOF) return false;
   int D = fgetc(InF); if (D == EOF) return false;
   if (NP != NULL) *NP = (D << 24) | (C << 16) | (B << 8) | A;
   return true;
}

static bool GetHeader(FILE *InF, long *XsP, long *YsP) {
   Num2 Sig; Num4 FSize, XX, Offset;
   Num4 ILen, Xs, Ys; Num2 Zs, Bits; Num4 Mode, ISize;
   Num4 ResX, ResY, Colors, MinColors;
   if (!(
// Header.
      GetNum2(InF, &Sig) && GetNum4(InF, &FSize) &&
      GetNum4(InF, &XX) && GetNum4(InF, &Offset) &&
// Specs.
      GetNum4(InF, &ILen) && GetNum4(InF, &Xs) && GetNum4(InF, &Ys) &&
      GetNum2(InF, &Zs) && GetNum2(InF, &Bits) && GetNum4(InF, &Mode) && GetNum4(InF, &ISize) &&
      GetNum4(InF, &ResX) && GetNum4(InF, &ResY) && GetNum4(InF, &Colors) && GetNum4(InF, &MinColors)
   )) {
      fprintf(stderr, "Inconsistent header format.\n"); return false;
   }
// Consistency checks
   if (Sig != (('M' << 8) | 'B')) { fprintf(stderr, "Invalid signature (%c%c).\n", (char)(Sig&0xff), (char)((Sig >> 8)&0xff)); return false; }
   if (XX != 0) { fprintf(stderr, "Invalid reserved bytes (%ld != 0).\n", (long)XX); return false; }
   if (ILen != 40) { fprintf(stderr, "Invalid info length (%ld != 40).\n", (long)ILen); return false; }
   if (Xs <= 0) { fprintf(stderr, "Inconsistent image width %ld < 0.\n", (long)Xs); return false; }
   if (Ys == 0) {
      fprintf(stderr, "Inconsistent image height %ld.\n", (long)Ys); return false;
   } else if (Ys < 0) {
      fprintf(stderr, "Downward image format (%ld < 0) not supported.\n", (long)Ys); return false;
   }
   if (Zs != 1) { fprintf(stderr, "Colors planes %d != 1 not supported.\n", (int)Zs); return false; }
   if (Bits != 24) { fprintf(stderr, "Color bits %d != 24 not supported.\n", (int)Bits); return false; }
   Num4 RowN = (3*Xs + 3)/4*4, Size = RowN*Ys;
   switch (Mode) {
      case 0: break; // RGB is supported here.
      case 1: fprintf(stderr, "RLE 8 format not supported here.\n"); return false;
      case 2: fprintf(stderr, "RLE 4 format not supported here.\n"); return false;
      default: fprintf(stderr, "Inconsistent encoding mode (%ld).\n", Mode); return false;
   }
   if (FSize != 54L + Size) { fprintf(stderr, "Inconsistent file size spec (%ld != %ld).\n", FSize, 54L + Size); return false; }
   if (Colors < 0) {
      fprintf(stderr, "Inconsistent color spec (%ld colors < 0).\n", Colors); return false;
   } else if (Colors > 0) {
      fprintf(stderr, "No palette support provided here (%ld colors > 0).\n", Colors); return false;
   }
   if ((Offset -= 54L) < 0) { fprintf(stderr, "Inconsistent data offset (off by %ld).\n", Offset); return false; }
   if (MinColors < 0) { fprintf(stderr, "Inconsistent min colors (%ld colors < 0).\n", MinColors); return false; }
   if (Offset > 0) {
      fprintf(stderr, "Skipping unknown data.\n");
      while (Offset-- > 0) if (fgetc(InF) == EOF) {
         fprintf(stderr, "Inconsistent header length.\n"); return false;
      }
   }
   if (XsP != NULL) *XsP = Xs;
   if (YsP != NULL) *YsP = Ys;
   return true;
}

FILE *GetFile(char *Path, long *XsP, long *YsP) {
   FILE *InF = Path == NULL? stdin: fopen(Path, "rb");
   if (InF == NULL) fprintf(stderr, "Cannot get %s.\n", Path);
   else if (!GetHeader(InF, XsP, YsP)) {
      if (Path != NULL) fclose(InF);
      InF = NULL;
   }
   return InF;
}

bool GetRow(FILE *InF, Color Row, long Xs) {
   for (long X = 0; X < Xs; X++) {
      Color C = &Row[X];
      int B = fgetc(InF); if (B == EOF) return false; else C->B = B;
      int G = fgetc(InF); if (G == EOF) return false; else C->G = G;
      int R = fgetc(InF); if (R == EOF) return false; else C->R = R;
   }
   long RowN = (3*Xs + 3)/4*4;
   for (long I = 3*Xs; I < RowN; I++) if (fgetc(InF) == EOF) return false;
   return true;
}
