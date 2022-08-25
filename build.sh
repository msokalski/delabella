#CPP=clang++
CPP=g++

# fp-math must conform to strict IEEE behaviour, ie: don't use -ffast-math
# adding -I for clang
OPT="-std=c++17 -O3 -I/usr/include/SDL2"
LNK=""

if [ -d "fade" ]; then
    OPT="$OPT -DWITH_FADE"
    LNK="-Lfade/lib_ubuntu20.04_x86_64 -lfade2d -Wl,-rpath=fade/lib_ubuntu20.04_x86_64"
fi

if [ -d "delaunator-cpp" ]; then
    OPT="$OPT -DWITH_DELAUNATOR"
fi

if [ -d "CDT" ]; then
    OPT="$OPT -DWITH_CDT"
fi

if [[ $OSTYPE == 'darwin'* ]]; then
    OPT="$OPT -framework OpenGL"
fi    

$CPP $OPT delabella.cpp delabella-sdl2.cpp $LNK -lSDL2 -lGL -o delabella-sdl2
