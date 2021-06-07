#if !defined BOOL_H
#   define BOOL_H
typedef enum { false, true } bool;
#endif
#if !defined COLOR_H
#   define COLOR_H
typedef struct Color { int B, G, R; } *Color;
#endif

// Generic BMP encoder.
FILE *PutFile(char *Path, long Xs, long Ys);
bool PutRow(FILE *ExF, Color Row, long Xs);
