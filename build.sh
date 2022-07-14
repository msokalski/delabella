#CPP=clang++
CPP=g++

# adding -I for clang
# adding -frounding-math to enable effect of calling fesetround 
OPT="-std=c++17 -O3 -I/usr/include/SDL2 -frounding-math"

if [ -d "delaunator" ]; then
    OPT="$OPT -DDELAUNATOR"
fi

if [ -d "CDT" ]; then
    OPT="$OPT -DCdt"
fi

if [[ $OSTYPE == 'darwin'* ]]; then
    OPT="$OPT -framework OpenGL"
fi    

# g++ $OPT delabella.cpp delabella-sdl2.cpp -lSDL2 -lGL -o delabella-sdl2

if [ -d "crude-xa" ]; then
    $CPP $OPT -DCRUDE_XA delabella.cpp delabella-sdl2.cpp crude-xa/src/crude-xa.c -lSDL2 -lGL -o delabella-xa-sdl2
fi
