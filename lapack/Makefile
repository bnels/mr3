
.c.o:
	g++ -03 -std=c++11 -c $<

dlaneg: dlaneg.c
	g++ -Wall -std=c++11 -lm -pedantic -o dlaneg dlaneg.c -llapack -lblas

clean:
	rm -f *.o *~ dlaneg
