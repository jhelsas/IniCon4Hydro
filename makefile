domain-spliting: domain-spliting.o splitandfit.o trial-functions.o src/splitandfit.h src/trial-functions.h
	g++ -o domain-spliting obj/domain-spliting.o obj/splitandfit.o obj/trial-functions.o -lm -lgsl -lgslcblas

gauss_whs: gauss_whs.o splitandfit.o src/splitandfit.h 
	g++ -o gauss_whs obj/gauss_whs.o obj/splitandfit.o -lm -lgsl -lgslcblas

gauss_e2s: gauss_e2s.o splitandfit.o src/splitandfit.h 
	g++ -o gauss_e2s obj/gauss_e2s.o obj/splitandfit.o -lm -lgsl -lgslcblas

gauss_qgp: gauss_qgp.o splitandfit.o src/splitandfit.h 
	g++ -o gauss_qgp obj/gauss_qgp.o obj/splitandfit.o -lm -lgsl -lgslcblas

gauss_zoltan: gauss_zoltan.o splitandfit.o conv_funct.o src/splitandfit.h 
	g++ -o gauss_zoltan obj/gauss_zoltan.o obj/splitandfit.o obj/conv_funct.o -lm -lgsl -lgslcblas

gubser: gubser.o splitandfit.o src/splitandfit.h 
	g++ -o gubser obj/gubser.o obj/splitandfit.o -lm -lgsl -lgslcblas

phsd_ico: phsd_ico.o splitandfit.o src/splitandfit.h 
	g++ -o phsd_ico obj/phsd_ico.o obj/splitandfit.o -lm -lgsl -lgslcblas `root-config --libs`

example: example.o splitandfit.o src/splitandfit.h
	g++ -o example obj/example.o obj/splitandfit.o -lm -lgsl -lgslcblas

cubic_dome: cubic_dome.o splitandfit.o src/splitandfit.h 
	g++ -o cubic_dome obj/cubic_dome.o obj/splitandfit.o -lm -lgsl -lgslcblas

sph-dens: splitandfit.o trial-functions.o sph-dens.o src/splitandfit.h src/trial-functions.h
	g++ -o sph-dens obj/sph-dens.o obj/splitandfit.o obj/trial-functions.o -lm -lgsl -lgslcblas

plot_densities: splitandfit.o trial-functions.o plot_densities.o src/splitandfit.h src/trial-functions.h
	g++ -o plot_densities obj/plot_densities.o obj/splitandfit.o obj/trial-functions.o -lm -lgsl -lgslcblas

gauss_whs.o: src/splitandfit.h src/trial-functions.h
	g++ -c src/examples/gauss_whs.cpp -O3 -o obj/gauss_whs.o

gauss_e2s.o: src/splitandfit.h src/trial-functions.h
	g++ -c src/examples/gauss_e2s.cpp -O3 -o obj/gauss_e2s.o

gauss_qgp.o: src/splitandfit.h src/trial-functions.h
	g++ -c src/examples/gauss_qgp.cpp -O3 -o obj/gauss_qgp.o

gauss_zoltan.o: src/splitandfit.h src/trial-functions.h
	g++ -c src/examples/gauss_zoltan.cpp -O3 -o obj/gauss_zoltan.o

gubser.o: src/splitandfit.h 
	g++ -c src/examples/gubser.cpp -O3 -o obj/gubser.o

phsd_ico.o: src/splitandfit.h 
	g++ -c src/examples/phsd_ico.cpp -O3 -o obj/phsd_ico.o `root-config --cflags`

example.o: src/splitandfit.h 
	g++ -c src/examples/example.cpp -O3 -o obj/example.o

cubic_dome.o: src/splitandfit.h 
	g++ -c src/examples/cubic_dome.cpp -O3 -o obj/cubic_dome.o

domain-spliting.o: src/splitandfit.h src/trial-functions.h
	g++ -c src/domain-spliting.cpp -O3 -o obj/domain-spliting.o

splitandfit.o: src/splitandfit.cpp src/splitandfit.h
	g++ -c src/splitandfit.cpp -O3 -o obj/splitandfit.o

trial-functions.o: src/splitandfit.h
	g++ -c src/trial-functions.cpp -O3 -o obj/trial-functions.o

sph-dens.o: src/splitandfit.h 
	g++ -c src/sph-dens.cpp -o obj/sph-dens.o

conv_funct.o: src/conv_funct.h
	g++ -c src/conv_funct.cpp -O3 -o obj/conv_funct.o

plot_densities.o:
	g++ -c src/plot_densities.cpp -o obj/plot_densities.o
