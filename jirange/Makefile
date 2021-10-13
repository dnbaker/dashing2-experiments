CXX?=g++

all: select_buckets safecat

%: %.cpp
	$(CXX) -O3 $< -march=native -o $@ -std=c++17
