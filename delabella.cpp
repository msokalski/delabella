/*
DELABELLA - Delaunay triangulation library
Copyright (C) 2018 GUMIX - Marcin Sokalski
*/

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include "delabella.h"

#ifdef IA_FAST
#include <fenv.h>
//anyone support this?
//#pragma STDC FENV_ACCESS ON
struct IA_FastRound
{
	const int restore;

	IA_FastRound() : restore(fegetround()) 
	{ 
		fesetround(FE_DOWNWARD); 
	}

	~IA_FastRound() 
	{ 
		fesetround(restore); 
	}
};
#else
struct IA_FastRound { IA_FastRound() {} };
#endif

static Unsigned28 s14sqr(const Signed14& s)
{
//	Signed29 m = 0x1p-1019;
//	Signed29 t = (Signed29)s * m;
//	Signed29 tt = t * t;
//	return tt;
	return (Unsigned28)((Signed29)s*(Signed29)s);
}

typedef DelaBella_Circumcenter Norm;

struct Vect
{
	Signed15 x, y;
	Signed29 z;

	Norm cross (const Vect& v) const // cross prod
	{
		Norm n;
		n.x = (Signed45)y*(Signed45)v.z - (Signed45)z*(Signed45)v.y;
		n.y = (Signed45)z*(Signed45)v.x - (Signed45)x*(Signed45)v.z;
		n.z = (Signed31)x*(Signed31)v.y - (Signed31)y*(Signed31)v.x;
		return n;
	}
};

IDelaBella::~IDelaBella()
{
}

struct CDelaBella : IDelaBella
{
	struct Face;

	struct Vert : DelaBella_Vertex
	{
		Unsigned28 z;

		Vect operator - (const Vert& v) const // diff
		{
			Vect d;
			d.x = (Signed15)x - (Signed15)v.x;
			d.y = (Signed15)y - (Signed15)v.y;
			d.z = (Signed29)z - (Signed29)v.z;
			return d;
		}

		static bool overlap(const Vert* v1, const Vert* v2)
		{
			return v1->x == v2->x && v1->y == v2->y;
		}

		bool operator < (const Vert& v) const
		{
			if (z < v.z)
				return true;
			if (z > v.z || x == v.x && y == v.y)
				return false;

			XA_REF ax = x;
			XA_REF ay = y;
			XA_REF az = ax * ax + ay * ay;
			XA_REF bx = v.x;
			XA_REF by = v.y;
			XA_REF bz = bx * bx + by * by;

			if (az < bz || az == bz && (x < v.x || x == v.x && y < v.y))
				return true;

			return false;
		}

		static int u28cmp(const void* a, const void* b)
		{
			const Vert* va = (const Vert*)a;
			const Vert* vb = (const Vert*)b;
			if (va->x* va->x + va->y * va->y < vb->x * vb->x + vb->y * vb->y)
				return -1;
			if (va->x * va->x + va->y * va->y > vb->x * vb->x + vb->y * vb->y)
				return 1;
			if (va->y < vb->y)
				return -1;
			if (va->y > vb->y)
				return 1;
			if (va->x < vb->x)
				return -1;
			if (va->x > vb->x)
				return 1;
			if (va->i < vb->i)
				return -1;
			if (va->i > vb->i)
				return 1;
			return 0;
		}
	};

	struct Face : DelaBella_Triangle
	{
		static Face* Alloc(Face** from)
		{
			Face* f = *from;
			*from = (Face*)f->next;
			f->next = 0;
			return f;
		}

		void Free(Face** to)
		{
			next = *to;
			*to = this;
		}

		Face* Next(const Vert* p) const
		{
			// return next face in same direction as face vertices are (cw/ccw)

			if (v[0] == p)
				return (Face*)f[1];
			if (v[1] == p)
				return (Face*)f[2];
			if (v[2] == p)
				return (Face*)f[0];
			return 0;
		}

		int sign() const
		{
			// try aprox
			if (n.z < 0)
				return -1;
			if (n.z > 0)
				return +1;

			// calc exact
			XA_REF x0 = v[0]->x, y0 = v[0]->y;
			XA_REF x1 = v[1]->x, y1 = v[1]->y;
			XA_REF x2 = v[2]->x, y2 = v[2]->y;

			XA_REF v1x = x1 - x0, v1y = y1 - y0;
			XA_REF v2x = x2 - x0, v2y = y2 - y0;

			return (v1x * v2y).cmp(v1y * v2x);
		}


		bool dot0(const Vert& p)
		{
			Vect pv = p - *(Vert*)v[0];
			Signed62 approx_dot =
				(Signed62)n.x * (Signed62)pv.x +
				(Signed62)n.y * (Signed62)pv.y +
				(Signed62)n.z * (Signed62)pv.z;

			if (approx_dot < 0 || approx_dot > 0)
				return false;

			XA_REF px = p.x;
			XA_REF py = p.y;

			XA_REF adx = (XA_REF)v[0]->x - px;
			XA_REF ady = (XA_REF)v[0]->y - py;
			XA_REF bdx = (XA_REF)v[1]->x - px;
			XA_REF bdy = (XA_REF)v[1]->y - py;
			XA_REF cdx = (XA_REF)v[2]->x - px;
			XA_REF cdy = (XA_REF)v[2]->y - py;

			return
				(adx * adx + ady * ady) * (cdx * bdy - bdx * cdy) +
				(bdx * bdx + bdy * bdy) * (adx * cdy - cdx * ady) ==
				(cdx * cdx + cdy * cdy) * (adx * bdy - bdx * ady);
		}

		bool dotP(const Vert& p)
		{
			Vect pv = p - *(Vert*)v[0];
			Signed62 approx_dot =
				(Signed62)n.x * (Signed62)pv.x +
				(Signed62)n.y * (Signed62)pv.y +
				(Signed62)n.z * (Signed62)pv.z;

			if (approx_dot > 0)
				return true;
			if (approx_dot < 0)
				return false;

			XA_REF px = p.x;
			XA_REF py = p.y;

			XA_REF adx = (XA_REF)v[0]->x - px;
			XA_REF ady = (XA_REF)v[0]->y - py;
			XA_REF bdx = (XA_REF)v[1]->x - px;
			XA_REF bdy = (XA_REF)v[1]->y - py;
			XA_REF cdx = (XA_REF)v[2]->x - px;
			XA_REF cdy = (XA_REF)v[2]->y - py;

			return
				(adx * adx + ady * ady) * (cdx * bdy - bdx * cdy) +
				(bdx * bdx + bdy * bdy) * (adx * cdy - cdx * ady) >
				(cdx * cdx + cdy * cdy) * (adx * bdy - bdx * ady); // re-swapped
		}

		bool dotN(const Vert& p)
		{
			Vect pv = p - *(Vert*)v[0];
			Signed62 approx_dot =
				(Signed62)n.x * (Signed62)pv.x +
				(Signed62)n.y * (Signed62)pv.y +
				(Signed62)n.z * (Signed62)pv.z;

			if (approx_dot < 0)
				return true;
			if (approx_dot > 0)
				return false;

			XA_REF px = p.x;
			XA_REF py = p.y;

			XA_REF adx = (XA_REF)v[0]->x - px;
			XA_REF ady = (XA_REF)v[0]->y - py;
			XA_REF bdx = (XA_REF)v[1]->x - px;
			XA_REF bdy = (XA_REF)v[1]->y - py;
			XA_REF cdx = (XA_REF)v[2]->x - px;
			XA_REF cdy = (XA_REF)v[2]->y - py;

			return
				(adx * adx + ady * ady) * (cdx * bdy - bdx * cdy) +
				(bdx * bdx + bdy * bdy) * (adx * cdy - cdx * ady) <
				(cdx * cdx + cdy * cdy) * (adx * bdy - bdx * ady); // re-swapped
		}

		// test
		bool dotNP(const Vert& p)
		{
			Vect pv = 
			{ 
				(Signed15)p.x - (Signed15)v[0]->x,
				(Signed15)p.y - (Signed15)v[0]->y,
				(Signed29)(*(Vert*)(v[0])).z - (Signed29)p.z // z filp!
			};
			Signed62 approx_dot_xy = (Signed62)n.x * (Signed62)pv.x + (Signed62)n.y * (Signed62)pv.y;
			Signed62 approx_dot_z = (Signed62)n.z * (Signed62)pv.z;

			#ifndef DB_AUTO_TEST
			if (approx_dot_xy < approx_dot_z)
				return true;
			if (approx_dot_xy > approx_dot_z)
				return false;
			#endif

			XA_REF px = p.x;
			XA_REF py = p.y;

			XA_REF adx = (XA_REF)v[0]->x - px;
			XA_REF ady = (XA_REF)v[0]->y - py;
			XA_REF bdx = (XA_REF)v[1]->x - px;
			XA_REF bdy = (XA_REF)v[1]->y - py;
			XA_REF cdx = (XA_REF)v[2]->x - px;
			XA_REF cdy = (XA_REF)v[2]->y - py;

			bool ret = 
				(adx * adx + ady * ady) * (cdx * bdy - bdx * cdy) +
				(bdx * bdx + bdy * bdy) * (adx * cdy - cdx * ady) <=
				(cdx * cdx + cdy * cdy) * (adx * bdy - bdx * ady); // re-swapped

			#ifdef DB_AUTO_TEST
			if (approx_dot_xy < approx_dot_z)
				assert(ret);
			else
			if (approx_dot_xy > approx_dot_z)
				assert(!ret);
			#endif

			return ret;
		}

		void cross() // cross of diffs
		{
			n = (*(Vert*)v[1] - *(Vert*)v[0]).cross(*(Vert*)v[2] - *(Vert*)v[0]);
		}
	};

	Vert* vert_alloc;
	Face* face_alloc;
	int* vert_map;
	int max_verts;
	int max_faces;

	Face* first_dela_face;
	Face* first_hull_face;
	Vert* first_boundary_vert;
	Vert* first_internal_vert;

	int inp_verts;
	int out_verts;
	int out_hull_faces;
	int out_boundary_verts;
	int unique_points;

	int(*errlog_proc)(void* file, const char* fmt, ...);
	void* errlog_file;

	int Prepare(int* start, Face** hull, int* out_hull_faces, Face** cache)
	{
		if (errlog_proc)
			errlog_proc(errlog_file, "[PRO] sorting vertices\n");
		int points = inp_verts;

		std::sort(vert_alloc, vert_alloc + points);

		// rmove dups
		{
			// it's fast!
			//if (errlog_proc)
			//	errlog_proc(errlog_file, "[PRO] looking for duplicates\n");

			int w = 0, r = 1; // skip initial no-dups block
			while (r < points && !Vert::overlap(vert_alloc + r, vert_alloc + w))
			{
				w++;
				r++;
			}

			int d = w; // dup map
			w++;

			while (r < points)
			{
				// fill map with dups only
				// unique verts will be filled after initial hull creation
				// which may require additional sorting on unique verts
				vert_map[vert_alloc[r].i] = d; // add first dup in run
				r++;

				// skip dups
				while (r < points && Vert::overlap(vert_alloc + r, vert_alloc + r - 1))
				{
					vert_map[vert_alloc[r].i] = d; // add next dup in run
					r++;
				}

				// copy next no-dups block (in percent chunks?)
				while (r < points && !Vert::overlap(vert_alloc + r, vert_alloc + r - 1))
					vert_alloc[w++] = vert_alloc[r++];

				d = w - 1;
			}

			if (points - w)
			{
				if (errlog_proc)
					errlog_proc(errlog_file, "[WRN] detected %d duplicates in xy array!\n", points - w);
				points = w;
			}
		}

		if (points < 3)
		{
			vert_map[vert_alloc[0].i] = 0;

			if (points == 2)
			{
				vert_map[vert_alloc[1].i] = 1;

				if (errlog_proc)
					errlog_proc(errlog_file, "[WRN] all input points are colinear, returning single segment!\n");
				first_boundary_vert = vert_alloc + 0;
				first_internal_vert = 0;
				vert_alloc[0].next = (DelaBella_Vertex*)vert_alloc + 1;
				vert_alloc[1].next = 0;
			}
			else
			{
				if (errlog_proc)
					errlog_proc(errlog_file, "[WRN] all input points are identical, returning signle point!\n");
				first_boundary_vert = vert_alloc + 0;
				first_internal_vert = 0;
				vert_alloc[0].next = 0;
			}

			out_boundary_verts = points;
			return -points;
		}

		// now with progress
		//if (errlog_proc)
		//	errlog_proc(errlog_file, "[PRO] preparing initial hull\n");

		int i;
		Face f; // tmp
		f.v[0] = vert_alloc + 0;

		Signed14 lo_x = vert_alloc[0].x, hi_x = lo_x;
		Signed14 lo_y = vert_alloc[0].y, hi_y = lo_y;
		int lower_left = 0;
		int upper_right = 0;
		for (i = 1; i < 3; i++)
		{
			Vert* v = vert_alloc + i;
			f.v[i] = v;

			if (v->x < lo_x)
			{
				lo_x = v->x;
				lower_left = i;
			}
			else
			if (v->x == lo_x && v->y < vert_alloc[lower_left].y)
				lower_left = i;

			if (v->x > hi_x)
			{
				hi_x = v->x;
				upper_right = i;
			}
			else
			if (v->x == hi_x && v->y > vert_alloc[upper_right].y)
				upper_right = i;

			lo_y = v->y < lo_y ? v->y : lo_y;
			hi_y = v->y > hi_y ? v->y : hi_y;
		}
		
		f.cross();


		int pro = 0;
		// skip until points are coplanar
		while (i < points && f.dot0(vert_alloc[i]))
		{
			if (i >= pro)
			{
				int p = (int)((uint64_t)100 * i / points);
				pro = (int)((uint64_t)(p + 1) * points / 100);
				if (pro >= points)
					pro = points - 1;
				if (i == points - 1)
					p = 100;
				if (errlog_proc)
					errlog_proc(errlog_file, "\r[%2d%s] convex hull triangulation%s", p, p >= 100 ? "" : "%", p >= 100 ? "\n" : "");
			}

			Vert* v = vert_alloc + i;

			if (v->x < lo_x)
			{
				lo_x = v->x;
				lower_left = i;
			}
			else
			if (v->x == lo_x && v->y < vert_alloc[lower_left].y)
				lower_left = i;

			if (v->x > hi_x)
			{
				hi_x = v->x;
				upper_right = i;
			}
			else
			if (v->x == hi_x && v->y > vert_alloc[upper_right].y)
				upper_right = i;

			lo_y = v->y < lo_y ? v->y : lo_y;
			hi_y = v->y > hi_y ? v->y : hi_y;

			i++;
		}

		bool colinear = f.sign() == 0; // hybrid
		if (colinear)
		{
			// choose x or y axis to sort verts (no need to be exact)
			if (hi_x - lo_x > hi_y - lo_y)
			{
				struct { bool operator()(const Vert& a, const Vert& b) const { return a.x < b.x; } } c;
				std::sort(vert_alloc, vert_alloc + i, c);
			}
			else
			{
				struct { bool operator()(const Vert& a, const Vert& b) const { return a.y < b.y; } } c;
				std::sort(vert_alloc, vert_alloc + i, c);
			}
		}
		else
		{
			// split to lower (below diag) and uppert (above diag) parts, 
			// sort parts separately in opposite directions
			// mark part with Vert::sew temporarily

			Signed15 part_x = (Signed15)vert_alloc[lower_left].y - (Signed15)vert_alloc[upper_right].y;
			Signed15 part_y = (Signed15)vert_alloc[upper_right].x - (Signed15)vert_alloc[lower_left].x;

			for (int j = 0; j < i; j++)
			{
				if (j == lower_left)
				{
					// lower
					vert_alloc[j].sew = &f;
				}
				else
				if (j == upper_right)
				{
					// upper
					vert_alloc[j].sew = 0;
				}
				else
				{
					Signed15 vx = (Signed15)vert_alloc[j].x - (Signed15)vert_alloc[lower_left].x;
					Signed15 vy = (Signed15)vert_alloc[lower_left].y - (Signed15)vert_alloc[j].y; // flipped!

					Signed31 dot_x = part_x * vx;
					Signed31 dot_y = part_y * vy;

					if (dot_x < dot_y)
					{
						// lower
						vert_alloc[j].sew = &f;
					}
					else
					if (dot_x > dot_y)
					{
						// upper
						vert_alloc[j].sew = 0;
					}
					else
					{
						// requires XA 
						// calc everything everytime from scatch?
						XA_REF llx = vert_alloc[lower_left].x;
						XA_REF lly = vert_alloc[lower_left].y;
						XA_REF urx = vert_alloc[upper_right].x;
						XA_REF ury = vert_alloc[upper_right].y;
						XA_REF xa_part_x = lly - ury;
						XA_REF xa_part_y = urx - llx;

						XA_REF vx = (XA_REF)vert_alloc[j].x - llx;
						XA_REF vy = lly - (XA_REF)vert_alloc[j].y; // flipped!

						XA_REF dot_x = xa_part_x * vx;
						XA_REF dot_y = xa_part_y * vy;

						if (dot_x < dot_y)
						{
							// lower
							vert_alloc[j].sew = &f;
						}
						else
						{
							#ifdef DB_AUTO_TEST
							assert(dot_x > dot_y);
							#endif
							// upper
							vert_alloc[j].sew = &f;
						}
					}
				}
			}

			part_x = (Signed15)vert_alloc[upper_right].x - (Signed15)vert_alloc[lower_left].x;
			part_y = (Signed15)vert_alloc[upper_right].y - (Signed15)vert_alloc[lower_left].y;

			struct
			{
				// default to CW order (unlikely, if wrong, we will reverse later)
				bool operator()(const Vert& a, const Vert& b) const
				{
					// if a is lower and b is upper, return true
					if (a.sew && !b.sew)
						return false;

					// if a is upper and b is lower, return false
					if (!a.sew && b.sew)
						return true;

					// calc a and b dot prods with (ur.x-ll.x,ur.y-ll.y)
					Signed15 ax = (Signed15)a.x - (Signed15)ll->x;
					Signed15 ay = (Signed15)a.y - (Signed15)ll->y;
					Signed15 dot_a = part_x * ax + part_y * ay;

					Signed15 bx = (Signed15)b.x - (Signed15)ll->x;
					Signed15 by = (Signed15)b.y - (Signed15)ll->y;
					Signed31 dot_b = part_x * bx + part_y * by;

					if (dot_a < dot_b)
						return a.sew == 0;
					if (dot_a > dot_b)
						return a.sew != 0;

					// xa needed
					{
						// calc everything everytime from scatch?
						XA_REF llx = ll->x;
						XA_REF lly = ll->y;
						XA_REF urx = ur->x;
						XA_REF ury = ur->y;
						XA_REF xa_part_x = urx - llx;
						XA_REF xa_part_y = ury - lly;

						XA_REF ax = (XA_REF)a.x - llx;
						XA_REF ay = (XA_REF)a.y - lly;
						XA_REF dot_a = xa_part_x * ax + xa_part_y * ay;

						XA_REF bx = (XA_REF)b.x - llx;
						XA_REF by = (XA_REF)b.y - lly;
						XA_REF dot_b = xa_part_x * bx + xa_part_y * by;

						if (dot_a < dot_b)
							return a.sew == 0;
						else
						{
							#ifdef DB_AUTO_TEST
							assert(dot_a > dot_b);
							#endif
							return a.sew != 0;
						}
					}

					// otherwise
					#ifdef DB_AUTO_TEST
					assert(0);
					#endif
					return false;
				}
				const Vert* ll;
				const Vert* ur;
				Signed15 part_x;
				Signed15 part_y;
			} c{ vert_alloc + lower_left, vert_alloc + upper_right, part_x, part_y };
			std::sort(vert_alloc, vert_alloc + i, c);
		}

		// fill map with unique verts
		for (int j = 0; j < points; j++)
			vert_map[vert_alloc[j].i] = j;

		*out_hull_faces = 0;

		// alloc faces only if we're going to create them
		if (i < points || !colinear)
		{
			int hull_faces = 2 * points - 4;
			*out_hull_faces = hull_faces;

			if (max_faces < hull_faces)
			{
				if (max_faces)
				{
					//free(face_alloc);
					delete[] face_alloc;
				}
				max_faces = 0;
				//face_alloc = (Face*)malloc(sizeof(Face) * hull_faces);
				face_alloc = new Face[hull_faces];
				if (face_alloc)
					max_faces = hull_faces;
				else
				{
					if (errlog_proc)
						errlog_proc(errlog_file, "[ERR] Not enough memory, shop for some more RAM. See you!\n");
					return 0;
				}
			}

			for (int i = 1; i < hull_faces; i++)
				face_alloc[i - 1].next = face_alloc + i;
			face_alloc[hull_faces - 1].next = 0;

			*cache = face_alloc;
		}

		if (i == points)
		{
			if (errlog_proc)
			{
				if (colinear)
					errlog_proc(errlog_file, "[WRN] all input points are colinear\n");
				else
					errlog_proc(errlog_file, "[WRN] all input points are cocircular\n");
			}

			if (colinear)
			{
				// link verts into open list
				first_boundary_vert = vert_alloc + 0;
				first_internal_vert = 0;
				out_boundary_verts = points;

				for (int j = 1; j < points; j++)
					vert_alloc[j - 1].next = (DelaBella_Vertex*)vert_alloc + j;
				vert_alloc[points - 1].next = 0;

				return -points;
			}

			// we're almost done!
			// time to build final flat hull
			Face* next_p = Face::Alloc(cache);
			Face* next_q = Face::Alloc(cache);
			Face* prev_p = 0;
			Face* prev_q = 0;

			for (int j = 2; j < i; j++)
			{
				Face* p = next_p;
				p->v[0] = vert_alloc + 0;
				p->v[1] = vert_alloc + j-1;
				p->v[2] = vert_alloc + j;
				p->cross();

				// mirrored
				Face* q = next_q;
				q->v[0] = vert_alloc + 0;
				q->v[1] = vert_alloc + j;
				q->v[2] = vert_alloc + j-1;
				q->cross();

				if (j < i - 1)
				{
					next_p = Face::Alloc(cache);
					next_q = Face::Alloc(cache);
				}
				else
				{
					next_p = 0;
					next_q = 0;
				}

				p->f[0] = q;
				p->f[1] = next_p ? next_p : q;
				p->f[2] = prev_p ? prev_p : q;

				q->f[0] = p;
				q->f[1] = prev_q ? prev_q : p;
				q->f[2] = next_q ? next_q : p;

				prev_p = p;
				prev_q = q;
			}

			*hull = prev_q;
		}
		else
		{
			// time to build cone hull with i'th vertex at the tip and 0..i-1 verts in the base
			// build cone's base in direction it is invisible to cone's tip!

			f.v[0] = vert_alloc + 0;
			f.v[1] = vert_alloc + 1;
			f.v[2] = vert_alloc + 2;
			f.cross();

			int one = 1, two = 2;
			if (!f.dotNP(vert_alloc[i]))
			{
				// if i-th vert can see the contour we will flip every face
				one = 2;
				two = 1;
			}

			Face* next_p = Face::Alloc(cache);
			Face* prev_p = 0;

			Face* first_q = 0;
			Face* next_q = Face::Alloc(cache);
			Face* prev_q = 0;

			for (int j = 2; j < i; j++)
			{
				Face* p = next_p;
				p->v[0] = vert_alloc + 0;
				p->v[one] = vert_alloc + j - 1;
				p->v[two] = vert_alloc + j;
				p->cross();

				Face* q;

				if (j == 2)
				{
					// first base triangle also build extra tip face
					q = Face::Alloc(cache);
					q->v[0] = vert_alloc + i;
					q->v[one] = vert_alloc + 1;
					q->v[two] = vert_alloc + 0;
					q->cross();

					q->f[0] = p;
					q->f[one] = 0; // LAST_Q;
					q->f[two] = next_q;

					first_q = q;
					prev_q = q;
				}

				q = next_q;
				q->v[0] = vert_alloc + i;
				q->v[one] = vert_alloc + j;
				q->v[two] = vert_alloc + j-1;
				q->cross();

				next_q = Face::Alloc(cache);

				q->f[0] = p;
				q->f[one] = prev_q;
				q->f[two] = next_q;
				prev_q = q;

				p->f[0] = q;

				if (j < i - 1)
				{
					next_p = Face::Alloc(cache);
				}
				else
				{
					// last base triangle also build extra tip face
					q = next_q;
					q->v[0] = vert_alloc + i;
					q->v[one] = vert_alloc + 0;
					q->v[two] = vert_alloc + i-1;
					q->cross();

					q->f[0] = p;
					q->f[one] = prev_q;
					q->f[two] = first_q;

					first_q->f[one] = q;

					next_p = 0;
					prev_q = q;
				}

				p->f[one] = next_p ? next_p : q;
				p->f[two] = prev_p ? prev_p : first_q;

				prev_p = p;
			}

			*hull = prev_q;

			i++;

			if (i == points)
			{
				if (errlog_proc)
					errlog_proc(errlog_file, "\r[%2d%s] convex hull triangulation%s", 100, "", "\n");
			}
		}

		*start = i;
		return points;
	}

	int Triangulate(int* other_faces)
	{
		int i = 0;
		Face* hull = 0;
		int hull_faces = 0;
		Face* cache = 0;

		int points = Prepare(&i, &hull, &hull_faces, &cache);
		unique_points = points < 0 ? -points : points;
		if (points <= 0)
			return points;

		/////////////////////////////////////////////////////////////////////////
		// ACTUAL ALGORITHM

		int pro = 0;
		for (; i < points; i++)
		{
			if (i >= pro)
			{
				int p = (int)((uint64_t)100 * i / points);
				pro = (int)((uint64_t)(p+1) * points / 100);
				if (pro >= points)
					pro = points - 1;
				if (i == points - 1)
					p = 100;
				if (errlog_proc)
					errlog_proc(errlog_file, "\r[%2d%s] convex hull triangulation%s", p, p>=100 ? "":"%", p>=100?"\n":"");
			}

			//ValidateHull(alloc, 2 * i - 4);
			Vert* q = vert_alloc + i;
			Vert* p = vert_alloc + i - 1;
			Face* f = hull;

			// 1. FIND FIRST VISIBLE FACE
			//    simply iterate around last vertex using last added triange adjecency info
			//while (f->dot(*q) <= 0)
			while (f->dotNP(*q))
			//while (f->dotN(*q)) // we want to consume coplanar faces
			{
				f = f->Next(p);
				if (f == hull)
				{
					//printf(".");
					// if no visible face can be located at last vertex,
					// let's run through all faces (approximately last to first),
					// yes this is emergency fallback and should not ever happen.
					f = face_alloc + 2 * i - 4 - 1;
					//while (f->dot(*q) <= 0)
					while (f->dotNP(*q))
					//while (f->dotN(*q)) // we want to consume coplanar faces
					{
						assert(f != face_alloc); // no face is visible? you must be kidding!
						f--;
					}
				}
			}

			// 2. DELETE VISIBLE FACES & ADD NEW ONES
			//    (we also build silhouette (vertex loop) between visible & invisible faces)

			int del = 0;
			int add = 0;

			// push first visible face onto stack (of visible faces)
			Face* stack = f;
			f->next = f; // old trick to use list pointers as 'on-stack' markers
			while (stack)
			{
				// pop, take care of last item ptr (it's not null!)
				f = stack;
				stack = (Face*)f->next;
				if (stack == f)
					stack = 0;
				f->next = 0;

				// copy parts of old face that we still need after removal
				Vert* fv[3] = { (Vert*)f->v[0],(Vert*)f->v[1],(Vert*)f->v[2] };
				Face* ff[3] = { (Face*)f->f[0],(Face*)f->f[1],(Face*)f->f[2] };

				// delete visible face
				f->Free(&cache);
				del++;

				// check all 3 neighbors
				for (int e = 0; e < 3; e++)
				{
					Face* n = ff[e];
					if (n && !n->next) // ensure neighbor is not processed yet & isn't on stack
					{
						// if neighbor is not visible we have slihouette edge
						//if (n->dot(*q) <= 0) 
						if (n->dotNP(*q))
						//if (n->dotN(*q)) // consuming coplanar faces
						{
							// build face
							add++;

							// ab: given face adjacency [index][],
							// it provides [][2] vertex indices on shared edge (CCW order)
							const static int ab[3][2] = { { 1,2 },{ 2,0 },{ 0,1 } };

							Vert* a = fv[ab[e][0]];
							Vert* b = fv[ab[e][1]];

							Face* s = Face::Alloc(&cache);
							s->v[0] = a;
							s->v[1] = b;
							s->v[2] = q;

							s->cross();
							s->f[2] = n;

							// change neighbour's adjacency from old visible face to cone side
							if (n->f[0] == f)
								n->f[0] = s;
							else
							if (n->f[1] == f)
								n->f[1] = s;
							else
							if (n->f[2] == f)
								n->f[2] = s;
							else
								assert(0);

							// build silhouette needed for sewing sides in the second pass
							a->sew = s;
							a->next = b;
						}
						else
						{
							// disjoin visible faces
							// so they won't be processed more than once

							if (n->f[0] == f)
								n->f[0] = 0;
							else
								if (n->f[1] == f)
									n->f[1] = 0;
								else
									if (n->f[2] == f)
										n->f[2] = 0;
									else
										assert(0);

							// push neighbor face, it's visible and requires processing
							n->next = stack ? stack : n;
							stack = n;
						}
					}
				}
			}

			// if add<del+2 hungry hull has consumed some point
			// that means we can't do delaunay for some under precission reasons
			// althought convex hull would be fine with it
			assert(add == del + 2);

			// 3. SEW SIDES OF CONE BUILT ON SLIHOUTTE SEGMENTS

			hull = face_alloc + 2 * i - 4 + 1; // last added face

										  // last face must contain part of the silhouette
										  // (edge between its v[0] and v[1])
			Vert* entry = (Vert*)hull->v[0];

			Vert* pr = entry;
			do
			{
				// sew pr<->nx
				Vert* nx = (Vert*)pr->next;
				pr->sew->f[0] = nx->sew;
				nx->sew->f[1] = pr->sew;
				pr = nx;
			} while (pr != entry);
		}

		assert(2 * i - 4 == hull_faces);
		//ValidateHull(alloc, hull_faces);

		for (int j = 0; j < points; j++)
		{
			vert_alloc[j].next = 0;
			vert_alloc[j].sew = 0;
		}

		if (errlog_proc)
			errlog_proc(errlog_file, "[PRO] postprocessing triangles\n");

		int others = 0;

		i = 0;
		Face** prev_dela = &first_dela_face;
		Face** prev_hull = &first_hull_face;
		pro = 0;
		for (int j = 0; j < hull_faces; j++)
		{
			if (j >= pro)
			{
				int p = (int)((uint64_t)100 * j / hull_faces);
				pro = (int)((uint64_t)(p + 1) * hull_faces / 100);
				if (pro >= hull_faces)
					pro = hull_faces - 1;
				if (j == hull_faces - 1)
					p = 100;
				if (errlog_proc)
					errlog_proc(errlog_file, "\r[%2d%s] postprocessing triangles%s", p, p >= 100 ? "" : "%", p >= 100 ? "\n" : "");
			}

			Face* f = face_alloc + j;

			// back-link all verts to some_face
			// yea, ~6x times, sorry
			((Vert*)f->v[0])->sew = f;
			((Vert*)f->v[1])->sew = f;
			((Vert*)f->v[2])->sew = f;

			bool nz_neg;
			// recalc triangle normal as accurately as possible on ALL triangles
			#if 0
			{
				XA_REF x0 = f->v[0]->x, y0 = f->v[0]->y, z0 = x0 * x0 + y0 * y0;
				XA_REF x1 = f->v[1]->x, y1 = f->v[1]->y/*, z1 = x1 * x1 + y1 * y1*/;
				XA_REF x2 = f->v[2]->x, y2 = f->v[2]->y/*, z2 = x2 * x2 + y2 * y2*/;

				XA_REF v1x = x1 - x0, v1y = y1 - y0, v1z = /*z1*/x1 * x1 + y1 * y1 - z0;
				XA_REF v2x = x2 - x0, v2y = y2 - y0, v2z = /*z2*/x2 * x2 + y2 * y2 - z0;

				// maybe attach these values temporarily to f
				// so we could sort triangles and build polygons?
				XA_REF nx = v1y * v2z - v1z * v2y;
				XA_REF ny = v1z * v2x - v1x * v2z;
				XA_REF nz = v1x * v2y - v1y * v2x;

				int exponent = (nx->quot + ny->quot + 2 * nz->quot + 2) / 4;
				nx->quot -= exponent;
				ny->quot -= exponent;
				nz->quot -= exponent;

				f->n.x = (double)nx;
				f->n.y = (double)ny;

				/*
				// needed?
				if (f->n.z < 0 && nz >= 0 || f->n.z >= 0 && nz < 0)
					printf("fixed!\n");
				*/

				f->n.z = (double)nz;

				/*
				// needed?
				if (nz != 0 && f->n.z == 0)
					f->n.z = nz->sign ? -FLT_MIN : FLT_MIN;

				if (nx != 0 && f->n.x == 0)
					f->n.x = nx->sign ? -FLT_MIN: FLT_MIN;

				if (ny != 0 && f->n.y == 0)
					f->n.y = ny->sign ? -FLT_MIN : FLT_MIN;
				*/

				nz_neg = nz < 0;
			}
			#else
			{
				// calc xa only where necessary
				nz_neg = f->sign() < 0; // hybrid
			}
			#endif


			if (nz_neg)
			{
				f->index = i;  // store index in dela list
				*prev_dela = f;
				prev_dela = (Face**)&f->next;
				i++;
			}
			else
			{
				f->index = ~others; // store index in hull list (~ to mark it is not dela)
				*prev_hull = f;
				prev_hull = (Face**)&f->next;
				others++;
			}
		}

		if (other_faces)
			*other_faces = others;

		*prev_dela = 0;
		*prev_hull = 0;

		if (errlog_proc)
			errlog_proc(errlog_file, "[PRO] postprocessing edges\n");

		// let's trace boudary contour, at least one vertex of first_hull_face
		// must be shared with dela face, find that dela face
		DelaBella_Iterator it;
		DelaBella_Vertex* v = first_hull_face->v[0];
		const DelaBella_Triangle* t = v->StartIterator(&it);
		const DelaBella_Triangle* e = t; // end
		
		first_boundary_vert = (Vert*)v;
		out_boundary_verts = 1;

		Signed14 wrap_x = v->x;
		Signed14 wrap_y = v->y;

		while (1)
		{
			if (t->index >= 0)
			{
				int pr = it.around-1; if (pr<0) pr = 2;
				int nx = it.around+1; if (nx>2) nx = 0;

				if (t->f[pr]->index < 0)
				{
					bool wrap = t->v[nx] == first_boundary_vert;
					// change x,y to the edge's normal as accurately as possible
					/*
					// THIS IS EVIL!
					#if 0
					{
						// exact edge normal
						XA_REF px = v->x, py = v->y;
						XA_REF vx = wrap ? wrap_x : t->v[nx]->x, vy = wrap ? wrap_y : t->v[nx]->y;
						XA_REF nx = py - vy;
						XA_REF ny = vx - px;
						int exponent = (nx->quot + ny->quot + 1) / 2;
						nx->quot -= exponent;
						ny->quot -= exponent;
						v->x = (double)nx;
						v->y = (double)ny;
					}
					#else
					{
						// approx edge normal
						Signed14 vx = wrap ? wrap_x : t->v[nx]->x, vy = wrap ? wrap_y : t->v[nx]->y;
						Signed14 nx = v->y - vy;
						Signed14 ny = vx - v->x;
						v->x = nx;
						v->y = ny;
					}
					#endif
					*/

					// let's move from: v to t->v[nx]
					v->next = t->v[nx];
					v = v->next;
					if (v == first_boundary_vert)
						break; // lap finished
					out_boundary_verts++;
					t = t->StartIterator(&it,nx);
					e = t;
					continue;
				}
			}
			t = it.Next();
			assert(t!=e);
		}

		if (errlog_proc)
			errlog_proc(errlog_file, "[PRO] postprocessing vertices\n");

		// link all other verts into internal list
		first_internal_vert = 0;
		Vert** prev_inter = &first_internal_vert;
		for (int j=0; j<points; j++)
		{
			if (!vert_alloc[j].next)
			{
				Vert* next = vert_alloc+j;
				*prev_inter = next;
				prev_inter = (Vert**)&next->next;
			}
		}

		return 3*i;
	}

	bool ReallocVerts(int points)
	{
		inp_verts = points;
		out_verts = 0;

		first_dela_face = 0;
		first_hull_face = 0;
		first_boundary_vert = 0;

		if (max_verts < points)
		{
			if (max_verts)
			{
				free(vert_map);
				vert_map = 0;

				//free(vert_alloc);
				delete [] vert_alloc;
				vert_alloc = 0;
				max_verts = 0;
			}

			vert_map = (int*)malloc(sizeof(int) * points);

			//vert_alloc = (Vert*)malloc(sizeof(Vert)*points);
			vert_alloc = new Vert[points];
			if (vert_alloc)
				max_verts = points;
			else
			{
				if (errlog_proc)
					errlog_proc(errlog_file, "[ERR] Not enough memory, shop for some more RAM. See you!\n");
				return false;
			}
		}

		return true;
	}

	// private
	Face* FindConstraintOffenders(Vert* va, Vert* vb, Face*** ptail, Vert** restart)
	{
		static const int rotate[3][3] = { {0,1,2},{1,2,0},{2,0,1} };
		static const int other_vert[3][2] = { {1,2},{2,0},{0,1} };

		// returns list of faces!
		// with Face::index replaced with vertex indice opposite to the offending edge
		// (we are about to change triangulation after all so indexes will change as well)
		// and it's safe for detrimination of dela/hull (sign is preserved)

		Face* list = 0;
		Face** tail = &list;

		DelaBella_Iterator it;

		int xa_ready = 0; // 0:none, 1:xa_b only, 2:all
		XA_REF xa_b[2];

	//restart:

		const DelaBella_Triangle* first = va->StartIterator(&it);
		const DelaBella_Triangle* face = first;

		if (xa_ready == 2)
			xa_ready = 1;

		IA_VAL ia_a[2] = { (IA_VAL)va->x , (IA_VAL)va->y };
		IA_VAL ab[2] = { (IA_VAL)vb->x - ia_a[0], (IA_VAL)vb->y - ia_a[1] };

		XA_REF xa_a[2];
		XA_REF xa_ab[2];

		// find first face around va, containing offending edge
		Vert* v0;
		Vert* v1;
		Face* N = 0;
		int a, b, c;

		while (1)
		{
			if (face->index < 0)
			{
				face = it.Next();
				assert(face != first);
				continue;
			}

			a = it.around;
			b = other_vert[a][0];
			v0 = (Vert*)(face->v[b]);
			if (v0 == vb)
			{
				*tail = 0;
				*ptail = list ? tail : 0;
				return list; // ab is already there
			}

			c = other_vert[a][1];
			v1 = (Vert*)(face->v[c]);
			if (v1 == vb)
			{
				*tail = 0;
				*ptail = list ? tail : 0;
				return list; // ab is already there
			}

			IA_VAL a0[2] = { (IA_VAL)v0->x - ia_a[0], (IA_VAL)v0->y - ia_a[1] };
			IA_VAL a1[2] = { (IA_VAL)v1->x - ia_a[0], (IA_VAL)v1->y - ia_a[1] };

			IA_VAL l_a0b = a0[0] * ab[1];
			IA_VAL r_a0b = a0[1] * ab[0];
			IA_VAL l_a1b = a1[0] * ab[1];
			IA_VAL r_a1b = a1[1] * ab[0];


			if (l_a0b < r_a0b && l_a1b > r_a1b)
			{
				// offending edge!
				N = (Face*)face;
				break;
			}

			if (!(l_a0b > r_a0b || l_a1b < r_a1b))
			{
				if (xa_ready < 2)
				{
					if (xa_ready < 1)
					{
						xa_b[0] = (XA_REF)vb->x;
						xa_b[1] = (XA_REF)vb->y;
					}
					xa_a[0] = (XA_REF)va->x;
					xa_a[1] = (XA_REF)va->y;
					xa_ab[0] = xa_b[0] - xa_a[0];
					xa_ab[1] = xa_b[1] - xa_a[1];
					xa_ready = 2;
				}

				XA_REF xa_a0[2] = { (XA_REF)v0->x - xa_a[0], (XA_REF)v0->y - xa_a[1] };
				XA_REF xa_a1[2] = { (XA_REF)v1->x - xa_a[0], (XA_REF)v1->y - xa_a[1] };

				XA_REF xa_a0b = xa_a0[0] * xa_ab[1] - xa_a0[1] * xa_ab[0];
				XA_REF xa_a1b = xa_a1[0] * xa_ab[1] - xa_a1[1] * xa_ab[0];

				if (xa_a0b == 0)
				{
					// v0 overlap
					
					//va = v0;
					//goto restart;
					
					*restart = v0;
					*tail = 0;
					*ptail = list ? tail : 0;
					return list;
				}

				if (xa_a1b == 0)
				{
					// v1 overlap

					//va = v1;
					//goto restart;

					*restart = v1;
					*tail = 0;
					*ptail = list ? tail : 0;
					return list;
				}

				if (xa_a0b < 0 && xa_a1b > 0)
				{
					// offending edge!
					N = (Face*)face;
					break;
				}
			}

			face = it.Next();
			assert(face != first);
		}

		while (1)
		{
			if (a)
			{
				// rotate N->v[] and N->f 'a' times 'backward' such offending edge appears opposite to v[0]
				const int* r = rotate[a];

				DelaBella_Vertex* v[3] = { N->v[0] ,N->v[1] ,N->v[2] };
				N->v[0] = v[r[0]];
				N->v[1] = v[r[1]];
				N->v[2] = v[r[2]];

				DelaBella_Triangle* f[3] = { N->f[0] ,N->f[1] ,N->f[2] };
				N->f[0] = f[r[0]];
				N->f[1] = f[r[1]];
				N->f[2] = f[r[2]];
			}

			// add edge
			*tail = N;
			tail = (Face**)&N->next;

			// what is our next face?
			Face* F = (Face*)(N->f[0]);
			int d, e, f;

			if (F->f[0] == N)
			{
				d = 0;
				e = 1;
				f = 2;
			}
			else
			if (F->f[1] == N)
			{
				d = 1;
				e = 2;
				f = 0;
			}
			else
			{
				d = 2;
				e = 0;
				f = 1;
			}

			Vert* vr = (Vert*)(F->v[d]);

			if (vr == vb)
			{
				*restart = 0;
				*tail = 0;
				*ptail = list ? tail : 0;
				return list;
			}

			// is vr above or below ab ?

			IA_VAL ar[2] = { (IA_VAL)vr->x - ia_a[0], (IA_VAL)vr->y - ia_a[1] };

			IA_VAL l_abr = ab[0] * ar[1];
			IA_VAL r_abr = ab[1] * ar[0];

			if (l_abr > r_abr)
			{
				// above: de edge (a' = f vert)
				a = f;
				v0 = (Vert*)F->v[d];
				v1 = (Vert*)F->v[e];
				N = F;
			}
			else
			if (l_abr < r_abr)
			{
				// below: fd edge (a' = e vert)
				a = e;
				v0 = (Vert*)F->v[f];
				v1 = (Vert*)F->v[d];
				N = F;
			}
			else
			{
				// XA NEEDED, could be above, below or overlap!
				// if overlap: do restart on vr
				if (xa_ready < 2)
				{
					if (xa_ready < 1)
					{
						xa_b[0] = (XA_REF)vb->x;
						xa_b[1] = (XA_REF)vb->y;
					}
					xa_a[0] = (XA_REF)va->x;
					xa_a[1] = (XA_REF)va->y;
					xa_ab[0] = xa_b[0] - xa_a[0];
					xa_ab[1] = xa_b[1] - xa_a[1];
					xa_ready = 2;
				}

				XA_REF xa_ar[2] = { (XA_REF)vr->x - xa_a[0], (IA_VAL)vr->y - xa_a[1] };
				XA_REF xa_abr = xa_ab[0] * xa_ar[1] - xa_ab[1] * xa_ar[0];

				if (xa_abr == 0)
				{
					// if overlap
					
					//va = vr;
					//goto restart;

					*restart = vr;
					*tail = 0;
					*ptail = list ? tail : 0;
					return list;
				}

				if (xa_abr > 0)
				{
					// above: de edge (a' = f vert)
					a = f;
					v0 = (Vert*)F->v[d];
					v1 = (Vert*)F->v[e];
					N = F;
				}
				else
				{
					// below: fd edge (a' = e vert)
					a = e;
					v0 = (Vert*)F->v[f];
					v1 = (Vert*)F->v[d];
					N = F;
				}
			}
		}

		assert(0);
		*restart = 0;
		*ptail = 0;
		return 0;
	}

	virtual int Constrain(int num, const int* pa, const int* pb, int advance_bytes)
	{
		IA_FastRound round;

		if (advance_bytes <= 0)
			advance_bytes = 2*sizeof(int);

		int flips = 0;

		int pro = 0;
		for (int con = 0; con < num; con++)
		{
			if (con >= pro)
			{
				int p = (int)((uint64_t)100 * con / num);
				pro = (int)((uint64_t)(p + 1) * num / 100);
				if (pro >= num)
					pro = num - 1;
				if (con == num - 1)
					p = 100;
				if (errlog_proc)
					errlog_proc(errlog_file, "\r[%2d%s] constraining%s", p, p >= 100 ? "" : "%", p >= 100 ? "\n" : "");
			}

			int a = *(const int*)((const char*)pa + con * advance_bytes);
			int b = *(const int*)((const char*)pb + con * advance_bytes);

			if (!first_dela_face || a == b)
				continue;

			// current
			Vert* va = (Vert*)GetVertexByIndex(a);
			if (!va)
				continue;

			// destination
			Vert* vc = (Vert*)GetVertexByIndex(b);
			if (!vc)
				continue;

			do {
			// 1. Make list of offenders
			Face** tail = 0;
			Vert* restart = 0;
			Face* list = FindConstraintOffenders(va, vc, &tail, &restart);
			//if (!list && restart)
			//{
			//	va = restart;
			//	continue;
			//}

			Face* flipped = 0;

			Vert* vb = restart ? restart : vc;

			// will we relink'em back to dela and replace indexes?

			// 2. Repeatedly until there are no offenders on the list
			// - remove first edge from the list of offenders
			// - if it forms concave quad diagonal, push it back to the list of offenders
			// - otherwise do flip then:
			//   - if flipped diagonal still intersects ab, push it back to the list of offenders
			//   - otherwise push it to the list of flipped edges
			while (list)
			{
				const int a = 0, b = 1, c = 2;

				Face* N = list;
				list = (Face*)N->next;
				N->next = 0;

				Vert* v0 = (Vert*)(N->v[b]);
				Vert* v1 = (Vert*)(N->v[c]);

				Face* F = (Face*)(N->f[a]);
				int d, e, f;

				if (F->f[0] == N)
				{
					d = 0;
					e = 1;
					f = 2;
				}
				else
				if (F->f[1] == N)
				{
					d = 1;
					e = 2;
					f = 0;
				}
				else
				{
					d = 2;
					e = 0;
					f = 1;
				}

				Vert* v = (Vert*)N->v[0]; // may be not same as global va (if we have had a skip)
				Vert* vr = (Vert*)(F->v[d]);

				// is va,v0,vr,v1 a convex quad?
				IA_VAL a0[2] = { (IA_VAL)v0->x - (IA_VAL)v->x, (IA_VAL)v0->y - (IA_VAL)v->y };
				IA_VAL a1[2] = { (IA_VAL)v1->x - (IA_VAL)v->x, (IA_VAL)v1->y - (IA_VAL)v->y };
				IA_VAL ar[2] = { (IA_VAL)vr->x - (IA_VAL)v->x, (IA_VAL)vr->y - (IA_VAL)v->y };

				IA_VAL l_a0r = a0[0] * ar[1];
				IA_VAL r_a0r = a0[1] * ar[0];
				IA_VAL l_a1r = a1[0] * ar[1];
				IA_VAL r_a1r = a1[1] * ar[0];

				if (l_a0r > r_a0r || l_a1r < r_a1r)
				{
					// CONCAVE CUNT!
					*tail = N;
					tail = (Face**)&N->next;
					continue;
				}

				if (!(l_a0r < r_a0r && l_a1r > r_a1r))
				{
					// XA NEEDED

					XA_REF xa_a[2] = { (XA_REF)v->x , (XA_REF)v->y };
					XA_REF xa_a0[2] = { (XA_REF)v0->x - xa_a[0], (XA_REF)v0->y - xa_a[1] };
					XA_REF xa_a1[2] = { (XA_REF)v1->x - xa_a[0], (XA_REF)v1->y - xa_a[1] };
					XA_REF xa_ar[2] = { (XA_REF)vr->x - xa_a[0], (XA_REF)vr->y - xa_a[1] };

					XA_REF xa_a0r = xa_a0[0] * xa_ar[1] - xa_a0[1] * xa_ar[0];
					XA_REF xa_a1r = xa_a1[0] * xa_ar[1] - xa_a1[1] * xa_ar[0];

					if (xa_a0r >= 0 || xa_a1r <= 0)
					{
						// CONCAVE CUNT!
						*tail = N;
						tail = (Face**)&N->next;
						continue;
					}
				}

				// it's convex, xa already checked
				// if (l_a0r < r_a0r && l_a1r > r_a1r)
				{
					if (f == 0)
					{
						/*           *                             *
							   va*  / \                      va*  / \
								  \/   \                        \/   \
								  /\ O  \                       /\ O  \
							   v /  \    \v0                 v /  \    \v0
								*----\----*---------*	      *----\----*---------*
							   / \a   \ b/f\      q/	     / \'-,c\   a\      q/
							  /   \  N \/   \  Q  /   -->   /   \f '-, N  \  Q  /
							 /  P  \   /\ F  \   /		   /  P  \  F '-, b\   /
							/p      \c/e \   d\ /		  /p      \e   \d'-,\ /
						   *---------*----+----*		 *---------*----\----*
									v1     \     vr		           v1    \    vr
											*vb                           *vb
						*/

						Face* O = (Face*)N->f[c];

						Face* P = (Face*)(N->f[b]);
						int p = P->f[0] == N ? 0 : P->f[1] == N ? 1 : 2;

						Face* Q = (Face*)(F->f[e]);
						int q = Q->f[0] == F ? 0 : Q->f[1] == F ? 1 : 2;

						// do flip
						N->v[a] = v0;
						N->v[b] = vr;
						N->v[c] = v;
						N->f[a] = F;
						N->f[b] = O;
						N->f[c] = Q;
						F->v[f] = v;  // from v0
						F->f[d] = P;  // from N
						F->f[e] = N;  // from Q
						P->f[p] = F;  // from N
						Q->f[q] = N;  // from F

						v0->sew = N;
						v1->sew = F;
					}
					else // e==0, (or d==0 but we don't care about it)
					{
						/*           *						       *
									/p\						      /p\
								   /   \					     /   \
								  /  P  \					    /  P  \
							   v /       \ v0                v /       \ v0
								*---------*          	      *---------*
							   / \a  N  b/f\         	     / \'-,e    f\
						 va*__/___\_____/___\__*vb --> va*__/___\b '-, F__\__*vb
							 /     \   /  F  \             /     \  N '-, d\
							/   O   \c/e     d\  		  /   O   \a    c'-,\
						   *---------*---------*		 *---------*---------*
								   v1 \       / vr		         v1 \       / vr
									   \  Q  /                       \  Q  /
										\   /						  \   /
										 \q/						   \q/
										  *							    *
						*/

						Face* O = (Face*)N->f[b];

						Face* P = (Face*)(N->f[c]);
						int p = P->f[0] == N ? 0 : P->f[1] == N ? 1 : 2;

						Face* Q = (Face*)(F->f[f]);
						int q = Q->f[0] == F ? 0 : Q->f[1] == F ? 1 : 2;

						// do flip
						N->v[a] = v1;
						N->v[b] = v;
						N->v[c] = vr;
						N->f[a] = F;
						N->f[b] = Q;
						N->f[c] = O;
						F->v[e] = v;  // from v1
						F->f[d] = P;  // from N
						F->f[f] = N;  // from Q
						P->f[p] = F;  // from N
						Q->f[q] = N;  // from F

						v0->sew = F;
						v1->sew = N;
					}

					flips++;

					// if v and vr are on strongly opposite sides of the edge
					// push N's edge back to offenders otherwise push to the new edges

					if (va == v || vr == vb)
					{

						// resolved!
						N->next = flipped;
						flipped = N;
					}
					else
					{
						IA_VAL ab[2] = { (IA_VAL)vb->x - (IA_VAL)va->x, (IA_VAL)vb->y - (IA_VAL)va->y };
						IA_VAL av[2] = { (IA_VAL)v->x - (IA_VAL)va->x, (IA_VAL)v->y - (IA_VAL)va->y };
						IA_VAL ar[2] = { (IA_VAL)vr->x - (IA_VAL)va->x, (IA_VAL)vr->y - (IA_VAL)va->y };

						IA_VAL l_ab_av = ab[0] * av[1];
						IA_VAL r_ab_av = ab[1] * av[0];
						IA_VAL l_ab_ar = ab[0] * ar[1];
						IA_VAL r_ab_ar = ab[1] * ar[0];

						if (l_ab_av < r_ab_av && l_ab_ar < r_ab_ar || l_ab_av > r_ab_av && l_ab_ar > r_ab_ar)
						{
							// resolved
							N->next = flipped;
							flipped = N;
						}
						else
						if (l_ab_av < r_ab_av && l_ab_ar > r_ab_ar || l_ab_av > r_ab_av && l_ab_ar < r_ab_ar)
						{
							// unresolved
							*tail = N;
							tail = (Face**)&N->next;
						}
						else
						{
							// XA NEEDED
							XA_REF xa_ab[2] = { (XA_REF)vb->x - (XA_REF)va->x, (XA_REF)vb->y - (XA_REF)va->y };
							XA_REF xa_av[2] = { (XA_REF)v->x - (XA_REF)va->x, (XA_REF)v->y - (XA_REF)va->y };
							XA_REF xa_ar[2] = { (XA_REF)vr->x - (XA_REF)va->x, (XA_REF)vr->y - (XA_REF)va->y };

							XA_REF xa_ab_av = xa_ab[0] * xa_av[1] - xa_ab[1] * xa_av[0];
							XA_REF xa_ab_ar = xa_ab[0] * xa_ar[1] - xa_ab[1] * xa_ar[0];

							if (xa_ab_av < 0 && xa_ab_ar < 0 || xa_ab_av > 0 && xa_ab_ar > 0)
							{
								// resolved
								N->next = flipped;
								flipped = N;
							}
							else
							if (xa_ab_av < 0 && xa_ab_ar > 0 || xa_ab_av > 0 && xa_ab_ar < 0)
							{
								// unresolved
								*tail = N;
								tail = (Face**)&N->next;
							}
						}
					}
				}
			}

			// update cross prods
			Face* N = flipped;
			while (N)
			{
				N->cross();
				Face* F = (Face*)N->f[0];
				F->cross(); // fix!
				N = (Face*)N->next;
			}

			// 3. Repeatedly until no flip occurs
			// for every edge from new edges list,
			// if 2 triangles sharing the edge violates delaunay criterion
			// do diagonal flip

			while (1)
			{
				bool no_flips = true;
				Face* N = flipped;
				while (N)
				{
					if (N->v[1] == va && N->v[2] == vb || N->v[1] == vb && N->v[2] == va)
					{
						N = (Face*)N->next;
						continue;
					}

					const int a = 0, b = 1, c = 2;

					Vert* v0 = (Vert*)(N->v[b]);
					Vert* v1 = (Vert*)(N->v[c]);

					Face* F = (Face*)(N->f[a]);
					int d, e, f;

					if (F->f[0] == N)
					{
						d = 0;
						e = 1;
						f = 2;
					}
					else
					if (F->f[1] == N)
					{
						d = 1;
						e = 2;
						f = 0;
					}
					else
					{
						d = 2;
						e = 0;
						f = 1;
					}

					Vert* v = (Vert*)N->v[0];
					Vert* vr = (Vert*)(F->v[d]);

					// fixed by calling F->cross()
					// bool np = N->dotP(*vr);
					// bool fp = F->dotP(*v);
					// assert(np && fp || !np && !fp);
					if (N->dotP(*vr) /* || F->dotP(*v)*/) 
					{
						no_flips = false;
						flips++;

						if (f == 0)
						{
							Face* O = (Face*)N->f[c];
							Face* P = (Face*)(N->f[b]);
							int p = P->f[0] == N ? 0 : P->f[1] == N ? 1 : 2;
							Face* Q = (Face*)(F->f[e]);
							int q = Q->f[0] == F ? 0 : Q->f[1] == F ? 1 : 2;

							// do flip
							N->v[a] = v0;
							N->v[b] = vr;
							N->v[c] = v;
							N->f[a] = F;
							N->f[b] = O;
							N->f[c] = Q;
							F->v[f] = v;  // from v0
							F->f[d] = P;  // from N
							F->f[e] = N;  // from Q
							P->f[p] = F;  // from N
							Q->f[q] = N;  // from F

							v0->sew = N;
							v1->sew = F;
						}
						else
						{
							Face* O = (Face*)N->f[b];
							Face* P = (Face*)(N->f[c]);
							int p = P->f[0] == N ? 0 : P->f[1] == N ? 1 : 2;
							Face* Q = (Face*)(F->f[f]);
							int q = Q->f[0] == F ? 0 : Q->f[1] == F ? 1 : 2;

							// do flip
							N->v[a] = v1;
							N->v[b] = v;
							N->v[c] = vr;
							N->f[a] = F;
							N->f[b] = Q;
							N->f[c] = O;
							F->v[e] = v;  // from v1
							F->f[d] = P;  // from N
							F->f[f] = N;  // from Q
							P->f[p] = F;  // from N
							Q->f[q] = N;  // from F

							v0->sew = F;
							v1->sew = N;
						}

						N->cross();
						F->cross();
					}

					N = (Face*)N->next;
				}

				if (no_flips)
					break;
			}

			va = restart;
			} while (va);
		}

		// clean up the mess we've made with dela faces list !!!
		int hull_faces = 2 * unique_points - 4;
		Face** tail = &first_dela_face;
		int index = 0;
		for (int i = 0; i < hull_faces; i++)
		{
			if (face_alloc[i].index >= 0)
			{
				*tail = face_alloc + i;
				tail = (Face**)&face_alloc[i].next;
				face_alloc[i].index = index++;
			}
		}
		*tail = 0;

		return flips;
	}

	virtual int Polygonize(const DelaBella_Triangle* poly[])
	{
		const DelaBella_Triangle** buf = 0;
		if (!poly)
		{
			buf = (const DelaBella_Triangle**)malloc(sizeof(const DelaBella_Triangle*) * out_verts / 3);
			poly = buf;
			if (!poly)
				return -1;
		}

		// clear poly indices;
		Face* f = first_dela_face;
		while (f)
		{
			f->index = 0x40000000;
			f = (Face*)f->next;
		}

		int num = 0;
		f = first_dela_face;
		while (f)
		{
			bool new_poly = true;
			Face* next = (Face*)f->next;
			for (int i = 0; i < 3; i++)
			{
				Face* a = (Face*)f->f[i];
				int index = a->index;

				if (index >= 0 && index < 0x40000000)
				{
					int j = 0;
					for (; j < 3; j++)
					{
						Vert* v = (Vert*)a->v[j];
						if (v!=f->v[0] && v!=f->v[1] && v!=f->v[2] && !f->dot0(*v))
							break;
					}

					if (j == 3)
					{
						int dest = f->index;
						if (dest < 0x40000000)
						{
							// merging polys !!!
							Face* m = (Face*)poly[index];
							while (m)
							{
								Face* n = (Face*)m->next;
								m->index = dest;
								m->next = (Face*)poly[dest];
								poly[dest] = m;
								m = n;
							}

							// fill the gap
							if (index < num - 1)
							{
								dest = index;
								index = num - 1;
								poly[dest] = 0;

								m = (Face*)poly[index];
								while (m)
								{
									Face* n = (Face*)m->next;
									m->index = dest;
									m->next = (Face*)poly[dest];
									poly[dest] = m;
									m = n;
								}
							}

							num--;
						}
						else
						{
							new_poly = false;
							// merge f with existing poly
							f->index = index;
							f->next = (Face*)poly[index];
							poly[index] = f;
						}
					}
				}
			}

			if (new_poly)
			{
				// start new poly
				f->next = 0;
				f->index = num;
				poly[num++] = f;
			}

			f = next;
		}

		#if 1
		// ALTER POST PROC:
		// re-order triangles in every polygon and re-order indices in faces:
		// - first face defines first 3 verts in contour with v[0],v[1],v[2]
		// - every next face define just 1 additional vertex at v[0]
		// thay can form fan or strip or arbitraty mix of both
		// without changing existing edges!

		for (int p = num-1; p >= 0; p--)
		{
			Face* f = (Face*)(poly[p]);
			if (!f->next)
			{
				// single triangle
				// just link into list of polys
				if (p < num - 1)
					f->next = (Face*)poly[p + 1];
				continue;
			}

			Face* first = 0;
			// break list appart
			// so we can unmark all poly faces

			int dbg_num = 0;
			while (f)
			{
				Face* n = (Face*)f->next;

				if (!first)
				{
					// lookup good starting face (it will appear as last in the poly after all)
					// having exactly 2 poly boundary edges, exactly 1 inner edge

					int inner_edges =
						(f->f[0]->index == p) +
						(f->f[1]->index == p) +
						(f->f[2]->index == p);

					if (inner_edges == 1)
						first = f;
				}

				f->next = 0; // unmark
				f = n;

				dbg_num++;
			}

			assert(first);

			Face* list = first; // here list is empty, first is used as sentinel

			f = first; // current face

			Face* last = 0; // will be one inserted right after first
			
			bool step_on = false; // is current vertex inserted

			DelaBella_Iterator it;
			f->StartIterator(&it, f->f[0]->index == p ? 2 : f->f[1]->index == p ? 0 : 1);

			int dbg = 0;

			while (1)
			{
				if (!step_on && !f->next)
				{
					if (list == first)
						last = f;
					step_on = true;
					f->next = list;
					list = f;
					dbg++;

					// rotate such it.around becomes 0
					if (it.around != 0)
					{
						Face* fr = (Face*)f->f[0];
						Vert* vr = (Vert*)f->v[0];

						if (it.around == 1)
						{
							f->f[0] = f->f[1];
							f->f[1] = f->f[2];
							f->f[2] = fr;
							f->v[0] = f->v[1];
							f->v[1] = f->v[2];
							f->v[2] = vr;
						}
						else // it.around == 2
						{
							f->f[0] = f->f[2];
							f->f[2] = f->f[1];
							f->f[1] = fr;
							f->v[0] = f->v[2];
							f->v[2] = f->v[1];
							f->v[1] = vr;
						}

						// adjust iterator after rot
						it.around = 0;
					}
				}

				static const int next_probe[] = { 1,2,0 };
				if (f->f[next_probe[it.around]]->index != p)
				{
					// step on other leg:
					// and restart iterator
					static const int other_leg[3] = { 2,0,1 };
					f->StartIterator(&it, other_leg[it.around]);
					step_on = false;
					continue;
				}

				f = (Face*)it.Next();
				if (f == first)
					break;
			}

			// rotate list so first will be wrapped back to head of the face list
			first->next = list;
			list = first;

			// link last face with first face in next poly
			last->next = (p < num - 1) ? (Face*)poly[p + 1] : 0;

			// store ordered list in poly
			poly[p] = list;

			assert(dbg == dbg_num);
		}

		#else
		// merge polys into single list
		// they are separated by indexes

		for (int i = 0; i < num-1; i++)
		{
			f = (Face*)poly[i];
			while (f->next)
				f = (Face*)f->next;
			f->next = (Face*)poly[i + 1];
		}
		#endif

		first_dela_face = (Face*)poly[0];

		if (buf)
			free(buf);

		return num;
	}


	virtual int Triangulate(int points, const float* x, const float* y = 0, int advance_bytes = 0)
	{
		if (!x)
			return 0;

		if (!y)
			y = x + 1;

		if (advance_bytes < static_cast<int>(sizeof(float) * 2))
			advance_bytes = sizeof(float) * 2;

		if (!ReallocVerts(points))
			return 0;

		IA_FastRound round;

		for (int i = 0; i < points; i++)
		{
			Vert* v = vert_alloc + i;
			v->i = i;
			v->x = (Signed14)*(const float*)((const char*)x + i*advance_bytes);
			v->y = (Signed14)*(const float*)((const char*)y + i*advance_bytes);
			v->z = s14sqr(v->x) + s14sqr(v->y);
			//v->e = fabs(v->x) + fabs(v->y) + v->z;
		}
		
		out_hull_faces = 0;
		unique_points = 0;
		out_verts = Triangulate(&out_hull_faces);
		return out_verts;
	}

	virtual int Triangulate(int points, const double* x, const double* y, int advance_bytes)
	{
		if (!x)
			return 0;
		
		if (!y)
			y = x + 1;
		
		if (advance_bytes < static_cast<int>(sizeof(double) * 2))
			advance_bytes = sizeof(double) * 2;

		if (!ReallocVerts(points))
			return 0;

		IA_FastRound round;

		for (int i = 0; i < points; i++)
		{
			Vert* v = vert_alloc + i;
			v->i = i;
			v->x = (Signed14)*(const double*)((const char*)x + i*advance_bytes);
			v->y = (Signed14)*(const double*)((const char*)y + i*advance_bytes);
			v->z = s14sqr(v->x) + s14sqr(v->y);
			//v->e = fabs(v->x) + fabs(v->y) + v->z;
		}

		out_hull_faces = 0;
		unique_points = 0;
		out_verts = Triangulate(&out_hull_faces);
		return out_verts;
	}

	virtual int Triangulate(int points, const long double* x, const long double* y, int advance_bytes)
	{
		if (!x)
			return 0;
		
		if (!y)
			y = x + 1;
		
		if (advance_bytes < static_cast<int>(sizeof(double) * 2))
			advance_bytes = sizeof(long double) * 2;

		if (!ReallocVerts(points))
			return 0;

		IA_FastRound round;

		for (int i = 0; i < points; i++)
		{
			Vert* v = vert_alloc + i;
			v->i = i;
			v->x = (Signed14)*(const long double*)((const char*)x + i*advance_bytes);
			v->y = (Signed14)*(const long double*)((const char*)y + i*advance_bytes);
			v->z = s14sqr(v->x) + s14sqr(v->y);
			//v->e = fabs(v->x) + fabs(v->y) + v->z;
		}

		out_hull_faces = 0;
		unique_points = 0;
		out_verts = Triangulate(&out_hull_faces);
		return out_verts;
	}	

	virtual void Destroy()
	{
		if (vert_map)
			free(vert_map);

		if (face_alloc)
		{
			//free(face_alloc);
			delete [] face_alloc;
		}
		if (vert_alloc)
		{
			//free(vert_alloc);
			delete [] vert_alloc;
		}
		delete this;
	}

	// num of points passed to last call to Triangulate()
	virtual int GetNumInputPoints() const
	{
		return inp_verts;
	}

	// num of verts returned from last call to Triangulate()
	virtual int GetNumOutputIndices() const
	{
		return out_verts;
	}

	virtual int GetNumOutputHullFaces() const
	{
		return out_hull_faces;
	}

	virtual int GetNumBoundaryVerts() const
	{
		return out_verts < 0 ? -out_verts : out_boundary_verts;
	}

	virtual int GetNumInternalVerts() const
	{
		return out_verts < 0 ? 0 : unique_points - out_boundary_verts;
	}

	virtual const DelaBella_Triangle* GetFirstDelaunayTriangle() const
	{
		return first_dela_face;
	}

	virtual const DelaBella_Triangle* GetFirstHullTriangle() const
	{
		return first_hull_face;
	}

	virtual const DelaBella_Vertex* GetFirstBoundaryVertex() const
	{
		return first_boundary_vert;
	}

	virtual const DelaBella_Vertex* GetFirstInternalVertex() const
	{
		return first_internal_vert;
	}

	virtual const DelaBella_Vertex* GetVertexByIndex(int i) const
	{
		if (i < 0 || i >= inp_verts)
			return 0;
		return vert_alloc + vert_map[i];
	}

	virtual void SetErrLog(int(*proc)(void* stream, const char* fmt, ...), void* stream)
	{
		errlog_proc = proc;
		errlog_file = stream;
	}
};

IDelaBella* IDelaBella::Create()
{
	CDelaBella* db = new CDelaBella;
	if (!db)
		return 0;

	db->vert_map = 0;
	db->vert_alloc = 0;
	db->face_alloc = 0;
	db->max_verts = 0;
	db->max_faces = 0;

	db->first_dela_face = 0;
	db->first_hull_face = 0;
	db->first_boundary_vert = 0;

	db->inp_verts = 0;
	db->out_verts = 0;
	db->out_hull_faces = 0;
	db->unique_points = 0;

	db->errlog_proc = 0;
	db->errlog_file = 0;

	return db;
}

