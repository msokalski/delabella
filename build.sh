g++ -O3 delabella.cpp delabella-sdl2.cpp -lSDL2 -lGL -o delabella-sdl2

if [ -d "crude-xa" ]; then
g++ -O3 -DCRUDE_XA delabella.cpp delabella-sdl2.cpp crude-xa/src/crude-xa.c -lSDL2 -lGL -o delabella-xa-sdl2
fi
