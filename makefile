CXX = hipcc
CXXFLAGS = -w -O3 --offload-arch=gfx1012
LDFLAGS = -lraylib

sim: barnes.cpp 
		$(CXX) $(CXXFLAGS) barnes.cpp -o sim $(LDFLAGS)
