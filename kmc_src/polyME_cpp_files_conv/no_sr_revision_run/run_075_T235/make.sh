module reset
module load gcc/10.2.0

g++ -O2 -g -std=c++17 -o program main.cpp reaction_functions.cpp save_variables.cpp
