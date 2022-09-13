#CPP=clang++
CPP=g++

# fp-math must conform to strict IEEE behaviour, ie: don't use -ffast-math
# adding -I for clang
OPT="-std=c++20 -O3 -I/usr/include/SDL2"
LNK="-lSDL2 -lGL"
SRC="delabella.cpp delabella-sdl2.cpp"
OUT="delabella-sdl2"

if [ -d "triangle" ]; then
    OPT="$OPT -DWITH_TRIANGLE"
    SRC="$SRC triangle/triangle.c"
fi

if [ -d "fade" ]; then
    OPT="$OPT -DWITH_FADE"
    LNK="-Lfade/lib_ubuntu20.04_x86_64 -lfade2d -Wl,-rpath=fade/lib_ubuntu20.04_x86_64 $LNK"
    # requires libgmp-dev
fi

if [ -d "delaunator-cpp" ]; then
    OPT="$OPT -DWITH_DELAUNATOR"
fi

if [ -d "CDT" ]; then
    OPT="$OPT -DWITH_CDT -DCDT_USE_BOOST"
fi

if [[ $OSTYPE == 'darwin'* ]]; then
    OPT="$OPT -framework OpenGL"
fi    

$CPP $OPT $SRC $LNK -o $OUT
