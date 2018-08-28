# delabella
## 2D Delaunay triangulation (dela) - super stable (bella!)

- Bunch of credits must go to David for inventing such beautiful algorithm (Newton Apple Wrapper):

  http://www.s-hull.org/

- It is pretty well described in this paper:

  https://arxiv.org/ftp/arxiv/papers/1602/1602.04707.pdf
  
## Example triangulation of signed distance field around some fancy font glyph

![Example](sample.png)

## delabella usage:

```c++

#include "delabella.h"
// ...

	// somewhere in your code ...

	int POINTS = 1000000;

	struct MyPoint
	{
		char something;
		float x;
		int something_else;
		float y;
		float foo[5];
	};

	MyPoint* cloud = new MyPoint[POINTS];

	srand(36341);

	// gen some random input
	for (int i = 0; i < POINTS; i++)
	{
		cloud[i].x = rand();
		cloud[i].y = rand();
	}

	IDelaBella* idb = IDelaBella::Create();

	int verts = idb->Triangulate(POINTS, &cloud->x, &cloud->y, sizeof(MyPoint));

	// if positive, all ok 
	if (verts>0)
	{
		int tris = verts / 3;
		const DelaBella_Triangle* dela = idb->GetFirstDelaunayTriangle();
		for (int i = 0; i<tris; i++)
		{
			// do something with dela triangle 
			// ...
			dela = dela->next;
		}
	}
	else
	{
		// no points given or all points are colinear
		// make emergency call ...
	}

	delete[] cloud;
	idb->Destroy();

	// ...


```
## performance comparison:
Output from speedtest.cpp

|         points |        QHULL |       S-HULL |    S-HULL-3D |    DELABELLA |
| --------------:| ------------:| ------------:| ------------:| ------------:|
|             10 |     0.040 ms |     0.013 ms |     0.017 ms |     0.003 ms |
|             25 |     0.098 ms |     0.034 ms |     0.033 ms |     0.007 ms |
|             50 |     0.193 ms |     0.069 ms |     0.063 ms |     0.017 ms |
|            100 |     0.400 ms |     0.144 ms |     0.136 ms |     0.038 ms |
|            250 |     1.076 ms |     0.366 ms |     0.339 ms |     0.113 ms |
|            500 |     2.169 ms |     0.761 ms |     0.665 ms |     0.209 ms |
|          1,000 |     4.274 ms |     1.627 ms |     1.419 ms |     0.541 ms |
|          2,500 |    12.269 ms |     4.455 ms |     4.093 ms |     1.550 ms |
|          5,000 |    23.248 ms |     9.392 ms |     8.585 ms |     3.253 ms |
|         10,000 |    49.383 ms |    20.740 ms |    19.684 ms |     7.072 ms |
|         25,000 |   153.298 ms |    66.618 ms |    58.248 ms |    19.122 ms |
|         50,000 |   286.513 ms |   130.383 ms |   117.455 ms |    42.017 ms |
|        100,000 |   614.543 ms |   292.067 ms |   243.510 ms |    90.936 ms |
|        250,000 |  1695.421 ms |   832.388 ms |   677.082 ms |   251.337 ms |
|        500,000 |  3687.803 ms |  1795.950 ms |  1417.546 ms |   539.161 ms |
|      1,000,000 |  7801.639 ms |  3971.826 ms |  3007.223 ms |  1160.063 ms |
|      2,500,000 | 20803.077 ms | 11282.402 ms |  8043.722 ms |  3159.174 ms |
|      5,000,000 | 72522.570 ms | 24885.967 ms | 20205.864 ms |  8888.995 ms |

![speedtest results in Log,Log scale](speedtest.png)
