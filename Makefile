all:
	mpicc serial.c -o serial -lm -O3 -std=c11
low-version:
	mpicc serial.c -o serial -lm -O3
check:
	sh ./check.sh
