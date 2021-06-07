#if !defined BOOL_H
#   define BOOL_H
typedef enum { false, true } bool;
#endif

// Generic output routine for WAV file format.
extern bool PutSound(char *Path, double **WZ, long Ws, long Zs, long Nu, long Bits);
