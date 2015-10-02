CC = cc
CFLAGS = -fPIC

libglen.a: smooth.o resize.o fft.o fft2d.o vec_fft2d.o fft3d.o four1.o f_fft.o f_fft2d.o f_fft3d.o f_four1.o spline.o \
	ran1.o splint.o clubksb.o cludcmp.o nrutil.o dft.o dft2.o f_dft.o fftshift.o clubksb_gm.o cludcmp_gm.o dsvdcmp.o     \
	dpythag.o dsvdcmp_gm.o gasdev.o m_gasdev.o gaussinv.o
	ar rv ../libglen.a smooth.o resize.o fft.o fft2d.o vec_fft2d.o fft3d.o four1.o f_fft.o f_fft2d.o f_fft3d.o \
		f_four1.o spline.o splint.o clubksb.o cludcmp.o dft.o dft2.o f_dft.o fftshift.o nrutil.o \
	        clubksb_gm.o cludcmp_gm.o dsvdcmp.o dpythag.o dsvdcmp_gm.o ran1.o gasdev.o \
		gaussinv.o m_gasdev.o; \
		ranlib ../libglen.a; rm *.o

spline.o splint.o gaussinv.o: nrutil.o
nrutil.o: nrutil.h
resize.o: fft2d.o fftshift.o 
smooth.o: fft2d.o
fft2d.o: fft.o
vec_fft2d.o: fft.o
fft.o: four1.o
f_fft2d.o: f_fft.o
f_fft.o: f_four1.o
dsvdcmp.o: dpythag.o
dsvdcmp_gm.o: dpythag.o

