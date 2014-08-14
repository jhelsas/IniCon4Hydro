domain-spliting: splitandfit.o domain-spliting.o trial-functions.o src/splitandfit.h src/trial-functions.h
	g++ -o domain-spliting obj/splitandfit.o obj/domain-spliting.o obj/trial-functions.o -lm -lgsl -lgslcblas

domain-spliting.o: src/splitandfit.h src/trial-functions.h
	g++ -c src/domain-spliting.cpp -O3 -o obj/domain-spliting.o

splitandfit.o: src/splitandfit.cpp src/splitandfit.h
	g++ -c src/splitandfit.cpp -O3 -o obj/splitandfit.o

trial-functions.o: src/splitandfit.h
	g++ -c src/trial-functions.cpp -O3 -o obj/trial-functions.o

sph-dens:
	g++ -o sph-dens src/sph-dens.cpp -lm
