/*
DELABELLA - Delaunay triangulation library
Copyright (C) 2018 GUMIX - Marcin Sokalski
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <tuple>
#include "delaunator-cpp.h"

#include "delabella.h"

#include "NewtonApple_hull3D.h"

extern "C"
{
#include "qhull/libqhull_r.h"
#include "qhull/qhull_ra.h"
}

#include <algorithm>
#include "s_hull_pro.h"

#include <windows.h>

bool pointSortPredicate(const Shx& a, const Shx& b)
{
	if (a.r < b.r)
		return true;
	else if (a.r > b.r)
		return false;
	else if (a.c < b.c)
		return true;
	else
		return false;
};

bool pointComparisonPredicate(const Shx& a, const Shx& b)
{
	return a.r == b.r && a.c == b.c;
}

void print_points(int n)
{
	if (n < 1000)
		printf("% 14d", n);
	else
	if (n < 1000000)
		printf("% 10d,%03d", n/1000, n%1000);
	else
	if (n < 1000000000)
		printf("% 6d,%03d,%03d", n / 1000000, n/1000 % 1000, n%1000);
}

struct MyPoint
{
	int x, y;
	void get_xy(double xy[2], void* cookie)
	{
		xy[0] = x;
		xy[1] = y;
	}
};

int main(int argc, char* argv[])
{
	LARGE_INTEGER t0, t1, fr;
	QueryPerformanceFrequency(&fr);
	double dt;

	int POINTS;

	const static int test[18] =
	{
		10,25,50,
		100,250,500,
		1000,2500,5000,
		10000,25000,50000,
		100000,250000,500000,
		1000000,2500000,5000000
	};

	int passes = sizeof(test) / sizeof(int);

	printf("|         points |        QHULL |       S-HULL |    S-HULL-3D |   DELAUNATOR |    DELABELLA |\n");
	printf("| --------------:| ------------:| ------------:| ------------:| ------------:| ------------:|\n");

	for (int pass = -5; pass < passes; pass++)
	{
		if (pass<0)
			POINTS = test[0];
		else
			POINTS = test[pass];

		int max_tris = 2 * POINTS - 5;
		int* abc = new int[3 * max_tris];
		double* xy = new double[2 * POINTS];

		double dt_qhull, dt_shull, dt_shull3d, dt_delaunator, dt_delabella;
		dt_qhull = dt_shull = dt_shull3d = dt_delaunator = dt_delabella = -1;

		srand(pass + 32109);

		for (int i = 0; i < 2 * POINTS; i++)
		{
			double quot = (rand() & 0x3FFF) - 8192; // +-8k (14bits)
			unsigned int frac = ((rand() & 0x7FFF) << 15) | (rand() & 0x7FFF); // 30 bits
			xy[i] = quot + frac / (double)(1 << 30);
		}

		if (1)
		{
			char hidden_options[] = " d n v H U Qb QB Qc Qf Qg Qi Qm Qr QR Qv Qx TR E V FC Fi Fo Ft Fp FV Q0 Q1 Q2 Q3 Q4 Q5 Q6 Q7 Q8 Q9 ";

			// int curlong, totlong; /* used !qh_NOmem */
			int numpoints = POINTS, dim = 3;

			coordT *points = (coordT*)malloc(sizeof(coordT)*POINTS * 3);
			for (int i = 0; i < POINTS; i++)
			{
				points[3 * i + 0] = xy[2 * i + 0];
				points[3 * i + 1] = xy[2 * i + 1];
				points[3 * i + 2] = xy[2 * i + 0] * xy[2 * i + 0] + xy[2 * i + 1] * xy[2 * i + 1];
			}


			boolT ismalloc = false;

			qhT qh_qh;
			qhT *qh = &qh_qh;

			qh_init_A(qh, stdin, stdout, stderr, argc, argv);  /* sets qh->qhull_command */

			qh->NOerrexit = False;
			qh_option(qh, "delaunay  Qbbound-last", NULL, NULL);
			qh->DELAUNAY = True;     /* 'd'   */
			qh->TRIangulate = True; /* Qt*/
			qh->SCALElast = True;    /* 'Qbb' */
			qh->KEEPcoplanar = True; /* 'Qc', to keep coplanars in 'p' */
			qh_checkflags(qh, qh->qhull_command, hidden_options);
			qh_initflags(qh, qh->qhull_command);
			// points = qh_readpoints(qh, &numpoints, &dim, &ismalloc);
			if (dim >= 5) {
				qh_option(qh, "Qxact_merge", NULL, NULL);
				qh->MERGEexact = True; /* 'Qx' always */
			}
			qh_init_B(qh, points, numpoints, dim, ismalloc);

			QueryPerformanceCounter(&t0);
			qh_qhull(qh);
			qh_findgood_all(qh, qh->facet_list);
			QueryPerformanceCounter(&t1);

			/*
			qh_check_output(qh);
			qh_produce_output(qh);
			if (qh->VERIFYoutput && !qh->FORCEoutput && !qh->STOPpoint && !qh->STOPcone)
			qh_check_points(qh);
			*/

			int tris = qh->num_good;
			dt = (double)(t1.QuadPart - t0.QuadPart) / (double)fr.QuadPart;
			//printf("QHL tris = %d, dt = %f ms\n", tris, dt * 1000);

			free(points);

			dt_qhull = dt;
		}

		if (1)
		{
			std::vector<Shx> pts;
			std::vector<Triad> triads;

			for (int i = 0; i < POINTS; i++)
			{
				Shx pt;
				pt.id = i;
				pt.r = (float)xy[2 * i + 0];
				pt.c = (float)xy[2 * i + 1];
				pts.push_back(pt);
			}

			QueryPerformanceCounter(&t0);
			std::sort(pts.begin(), pts.end(), pointSortPredicate);
			std::vector<Shx>::iterator newEnd = std::unique(pts.begin(), pts.end(), pointComparisonPredicate);
			pts.resize(newEnd - pts.begin());
			s_hull_pro(pts, triads);
			QueryPerformanceCounter(&t1);

			double dt = (double)(t1.QuadPart - t0.QuadPart) / (double)fr.QuadPart;
			//printf("SHL tris = %d, dt = %f ms\n", (int)triads.size(), dt * 1000);

			dt_shull = dt;
		}

		if (1)
		{
			std::vector<R3> pts;
			std::vector<Tri> triads;

			for (int i = 0; i < POINTS; i++)
			{
				R3 pt;
				pt.id = i;
				pt.r = xy[2 * i + 0];
				pt.c = xy[2 * i + 1];
				pt.z = pt.r*pt.r + pt.c*pt.c;
				pts.push_back(pt);
			}

			QueryPerformanceCounter(&t0);
			NewtonApple_Delaunay(pts, triads);
			QueryPerformanceCounter(&t1);

			dt = (double)(t1.QuadPart - t0.QuadPart) / (double)fr.QuadPart;
			//printf("SH3 tris = %d, dt = %f ms\n", (int)triads.size(), dt * 1000);

			dt_shull3d = dt;
		}

		if (1)
		{

			std::vector<double> coords;
			for (int i = 0; i < POINTS; i++)
			{
				coords.push_back(xy[2 * i + 0]);
				coords.push_back(xy[2 * i + 1]);
			}

			QueryPerformanceCounter(&t0);
			delaunator::Delaunator d(coords);
			size_t verts = d.triangles.size();
			QueryPerformanceCounter(&t1);


			dt = (double)(t1.QuadPart - t0.QuadPart) / (double)fr.QuadPart;
			//printf("DTR tris = %d, dt = %f ms\n", verts/3, dt * 1000);

			dt_delaunator = dt;
		}

		if (1)
		{
			IDelaBella* idb = IDelaBella::Create();
			
			struct MyErrLogStream
			{
				MyErrLogStream(FILE* f) : file(f) {}
				static int errlog(void* stream, const char* fmt, ...)
				{
					int ret;
					MyErrLogStream* s = (MyErrLogStream*)stream;
					va_list arglist;
					va_start(arglist, fmt);
					ret = vfprintf(s->file, fmt, arglist);
					va_end(arglist);
					return ret;
				}
				FILE* file;
			};

			MyErrLogStream myerrlog(stderr);
			idb->SetErrLog(MyErrLogStream::errlog, &myerrlog);

			QueryPerformanceCounter(&t0);
			int verts = idb->Triangulate(POINTS, xy, xy+1, sizeof(double[2]));
			QueryPerformanceCounter(&t1);

			idb->Destroy();

			dt = (double)(t1.QuadPart - t0.QuadPart) / (double)fr.QuadPart;
			//printf("DEL tris = %d, dt = %f ms\n", verts / 3, dt * 1000);

			dt_delabella = dt;
		}

		if (pass >= 0)
		{
			printf("| ");
			print_points(POINTS);
			printf(" | %9.3f ms | %9.3f ms | %9.3f ms | %9.3f ms | %9.3f ms |\n",
				dt_qhull * 1000,
				dt_shull * 1000,
				dt_shull3d * 1000,
				dt_delaunator * 1000,
				dt_delabella * 1000);
		}

		delete[] xy;
		delete[] abc;
	}

	return 0;
}



void a()
{
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
}