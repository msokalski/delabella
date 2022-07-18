#CPP=clang++
CPP=g++

# adding -I for clang
OPT="-std=c++17 -g -O1 -I/usr/include/SDL2 -ffast-math"

if [ -d "delaunator" ]; then
    OPT="$OPT -DDELAUNATOR"
fi

if [ -d "CDT" ]; then
    OPT="$OPT -DCdt"
fi

if [[ $OSTYPE == 'darwin'* ]]; then
    OPT="$OPT -framework OpenGL"
fi    

g++ $OPT delabella.cpp delabella-sdl2.cpp -lSDL2 -lGL -o delabella-sdl2
