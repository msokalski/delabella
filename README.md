# delabella
## 2D Delaunay triangulation (dela) - super stable (bella!)

- Bunch of credits must go to David for inventing such beautiful algorithm (Newton Apple Wrapper):

  http://www.s-hull.org/

- It is pretty well described in this paper:

  https://arxiv.org/ftp/arxiv/papers/1602/1602.04707.pdf

## delabella usage:

```c++

    int POINTS = 1000000;
    int* abc = new int[3 * (2 * POINTS - 4/*5*/)];
    double* xy = new double[2 * POINTS];

    srand(36341);

    for (int i = 0; i < 2 * POINTS; i++)
    {
      double quot = (rand() & 0x3FFF) - 8192; // +-8k (14bits)
      unsigned int frac = ((rand() & 0x7FFF) << 15) | (rand() & 0x7FFF); // 30 bits
      xy[i] = quot + frac / (double)(1<<30);
    }
  
    int verts = DelaBella(POINTS, xy, abc);

    // if positive, all ok 
    if (verts>0)
    {
      int tris = verts/3;
      // do something useful with abc array
      // ...
    }
    else
    {
      // no points given or all points are colinear
      // make emergency call ...
    }

    delete[] xy;
    delete[] abc;

```
