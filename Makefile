CC=gcc
CFLAGS=-O -g -std=gnu99
INSTALLDIR=${HOME}/bin
X=
O=.o
OBJS=arss$O EnBMP$O DeBMP$O EnWAV$O DeWAV$O
INSTALLDIR=/usr/local/bin
arss$X: ${OBJS}
	${CC} ${CFLAGS} ${OBJS} -lm -o arss$X
arss$0: EnBMP.h DeBMP.h EnWAV.h DeWAV.h
EnBMP$0: EnBMP.h
	$(CC) ${CFLAGS} -c EnBMP.c
DeBMP$0: DeBMP.h
	$(CC) ${CFLAGS} -c DeBMP.c
EnWAV$O: EnWAV.h
	$(CC) ${CFLAGS} -c EnWAV.c
DeWAV$O: DeWAV.h
	$(CC) ${CFLAGS} -c DeWAV.c
install: arss$X
	install -c -s arss$X ${INSTALLDIR}
clean:
	rm -f *$O core
clobber: clean
	rm -f arss$X
