#if !defined BOOL_H
#   define BOOL_H
typedef enum { false, true } bool;
#endif
#if !defined COLOR_H
#   define COLOR_H
typedef struct Color { int B, G, R; } *Color;
#endif

// Generic BMP decoder.
FILE *GetFile(char *Path, long *XsP, long *YsP);
bool GetRow(FILE *InF, Color Row, long Xs);
