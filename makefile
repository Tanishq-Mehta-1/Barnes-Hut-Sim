sim: barnes.cpp 
		g++ -w -O3 barnes.cpp -o sim -lraylib

q_tree: quadtree.cpp	
		g++ -w -O3 quadtree.cpp -o q_tree