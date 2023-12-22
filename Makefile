Projet_C : GJB_Extraction.o Generation.o main.o
	gcc -o Projet_C GJB_Extraction.o Generation.o main.o
Generation.o : Generation.c
	gcc -c Generation.c
GJB_Extraction.o : 	GJB_Extraction.c
	gcc -c GJB_Extraction.c
clean :
	rm -f Projet_C *.o
