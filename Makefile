SRC_FILES=main.c rng.c SFMT-src-1.5.1/SFMT.c

main: $(SRC_FILES)
	mpicc -DSFMT_MEXP=607 -Wall -Wextra -pedantic -std=c99 -o $@ $(SRC_FILES) -lm -O3

prof: $(SRC_FILES)
	mpicc -DSFMT_MEXP=607 -p -g -std=c99 -o $@ $(SRC_FILES) -lm -O0
	./prof 800 4 600 1e6
	gprof prof | head -24

clean:
	rm -f main prof gmon.out

run: main
	mpirun -np 2 main 800 4 600 1e8
