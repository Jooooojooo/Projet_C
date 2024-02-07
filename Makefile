# Compilateur
CC = gcc

# Options de compilation
CFLAGS = -O3 -fopenmp 
### Si sur MacOS avec homebrew, faire brew install libomp puis récupérer l'emplacement de téléchargement pour ensuite compiler.
### Exemple ci-dessous:
### CFLAGS = -O3 -Xpreprocessor -fopenmp -lomp -I/opt/homebrew/Cellar/libomp/17.0.6/include -L/opt/homebrew/Cellar/libomp/17.0.6/lib

# Fichiers source
SRCS = Version_NonParallelisee.c Version_Parallelisee.c Version_SBox.c

# Fichiers objets pour chaque exécutable
OBJS_EXEC1 = Version_NonParallelisee.o
OBJS_EXEC2 = Version_Parallelisee.o
OBJS_EXEC3 = Version_SBox.o

# Noms des exécutables
EXEC1 = Version_NonParallelisee
EXEC2 = Version_Parallelisee
EXEC3 = Version_SBox

# Règle par défaut pour construire tous les exécutables
all: $(EXEC1) $(EXEC2) $(EXEC3)

# Règle pour créer le fichier objet et créer le premier exécutable
$(EXEC1): $(OBJS_EXEC1)
	$(CC) $(CFLAGS) $(OBJS_EXEC1) -o $(EXEC1)

# Règle pour créer le fichier objet et créer le deuxième exécutable
$(EXEC2): $(OBJS_EXEC2)
	$(CC) $(CFLAGS) $(OBJS_EXEC2) -o $(EXEC2)

# Règle pour créer le fichier objet et créer le troisième exécutable
$(EXEC3): $(OBJS_EXEC3)
	$(CC) $(CFLAGS) $(OBJS_EXEC3) -o $(EXEC3)

# Règle pour compiler Version_NonParallelisee.c en Version_NonParallelisee.o
Version_NonParallelisee.o: Version_NonParallelisee.c
	$(CC) $(CFLAGS) -c Version_NonParallelisee.c -o Version_NonParallelisee.o

# Règle pour compiler Version_Parallelisee.c en Version_Parallelisee.o
Version_Parallelisee.o: Version_Parallelisee.c
	$(CC) $(CFLAGS) -c Version_Parallelisee.c -o Version_Parallelisee.o

# Règle pour compiler Version_SBox.c en Version_SBox.o
Version_SBox.o: Version_SBox.c
	$(CC) $(CFLAGS) -c Version_SBox.c -o Version_SBox.o

# Règle de nettoyage pour supprimer les fichiers objets et les exécutables
clean:
	rm -f $(OBJS_EXEC1) $(OBJS_EXEC2) $(OBJS_EXEC3) $(EXEC1) $(EXEC2) $(EXEC3)
