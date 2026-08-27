CXX = hipcc
CXXFLAGS = -w -O3 --offload-arch=gfx1012
LDFLAGS = -lraylib

sim: src/barnes.cpp 
		$(CXX) $(CXXFLAGS) src/barnes.cpp -o sim $(LDFLAGS)
