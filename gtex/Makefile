CXX?=g++

all: recount2csr recountcsr
recount2csr: recount2csr.cpp
	$(CXX) -I. recount2csr.cpp -o recount2csr -std=c++17 -march=native -I../bonsai/hll -O3 -march=native -lz
recountcsr: recount2csr
	cp $< $@
