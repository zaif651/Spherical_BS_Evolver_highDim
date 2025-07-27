mkdir obj
mkdir obj/Debug
g++-15 -Wall -fexceptions -g  -c "mathutils.cpp" -o obj/Debug/mathutils.o
g++-15 -Wall -fexceptions -g  -c "BosonStar.cpp" -o obj/Debug/BosonStar.o
g++-15 -Wall -fexceptions -g  -c "LinearPerturbation.cpp" -o obj/Debug/LinearPerturbation.o
g++-15 -Wall -fexceptions -g  -c "EvolutionVariables.cpp" -o obj/Debug/EvolutionVariables.o
g++-15 -Wall -fexceptions -g  -c "main.cpp" -o obj/Debug/main.o
g++-15  -o BosonStarStability.o obj/Debug/BosonStar.o obj/Debug/LinearPerturbation.o obj/Debug/EvolutionVariables.o obj/Debug/main.o obj/Debug/mathutils.o   -lm -lquadmath
