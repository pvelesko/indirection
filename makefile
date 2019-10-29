all: intel

intel: *.cpp *.hpp
	icpc -O2 -qopenmp -qopt-report=5 *.cpp -o intel

clean:
	rm -f *.o intel *.optrpt
