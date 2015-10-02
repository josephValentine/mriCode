CC = gcc
CFLAGS = -g

dbthigh: extern.h init.o dbthigh.o global.o punwrap.o driver_punwrap2.o cost.o fft.o DicomWrite.o
	$(CC) ${CFLAGS} -o dbthigh dbthigh.o init.o global.o driver_punwrap2.o punwrap.o cost.o fft.o DicomWrite.o -lm -lglen -L/user/glen/lib

