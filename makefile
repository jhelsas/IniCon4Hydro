domain-spliting: domain-spliting.o splitandfit.o trial-functions.o src/splitandfit.h src/trial-functions.h
	g++ -o domain-spliting obj/domain-spliting.o obj/splitandfit.o obj/trial-functions.o -lm -lgsl -lgslcblas

sph-dens: splitandfit.o trial-functions.o sph-dens.o src/splitandfit.h src/trial-functions.h
	g++ -o sph-dens obj/sph-dens.o obj/splitandfit.o obj/trial-functions.o -lm -lgsl -lgslcblas

domain-spliting.o: src/splitandfit.h src/trial-functions.h
	g++ -c src/domain-spliting.cpp -O3 -o obj/domain-spliting.o

splitandfit.o: src/splitandfit.cpp src/splitandfit.h
	g++ -c src/splitandfit.cpp -O3 -o obj/splitandfit.o

trial-functions.o: src/splitandfit.h
	g++ -c src/trial-functions.cpp -O3 -o obj/trial-functions.o

sph-dens.o: src/splitandfit.h 
	g++ -c src/sph-dens.cpp -o obj/sph-dens.o
