# GNBG-Instance-C
C++ implementation of Generalized Numerical Benchmark Generator (GNBG) instance for GECCO 2024 competition
Includes implementation of simple Differential Evolution (DE) with rand/1 strategy and binomial crossover
Problem parameters are read from f#.txt files which should be prepared with python script convert.py from f#.mat
Competition page: https://competition-hub.github.io/GNBG-Competition/
Reference:
D. Yazdani, M. N. Omidvar, D. Yazdani, K. Deb, and A. H. Gandomi, "GNBG: A Generalized
  and Configurable Benchmark Generator for Continuous Numerical Optimization," arXiv prepring	arXiv:2312.07083, 2023.
A. H. Gandomi, D. Yazdani, M. N. Omidvar, and K. Deb, "GNBG-Generated Test Suite for Box-Constrained Numerical Global
  Optimization," arXiv preprint arXiv:2312.07034, 2023.
MATLAB and Python versions: https://github.com/Danial-Yazdani/GNBG-Instances

# Compilation and usage
First convert f#.mat files into f#.txt with a specific format using python script:

python convert.py

Next compile the gnbg-c++.cpp file with any compiler (e.g. GCC):

g++ -std=c++11 -O3  gnbg-c++.cpp -o gnbg-c++

Running the executable will generate files "Res_DE_f#_r#.txt" for every function and run for further analysis.
