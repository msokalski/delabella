/*
DELABELLA - Delaunay triangulation library
Copyright (C) 2018 GUMIX - Marcin Sokalski
*/

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include "delabella.h"

static Unsigned28 s14sqr(const Signed14& s)
{
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

static int cmp_calls = 0;
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
			cmp_calls++;
			return u28cmp(this, &v) < 0;
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
			//always exact?
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
		printf("cmp_calls=%d\n", cmp_calls);
		cmp_calls = 0;

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

			w++;

			while (r < points)
			{
				r++;

				// skip dups
				while (r < points && Vert::overlap(vert_alloc + r, vert_alloc + r - 1))
					r++;

				// copy next no-dups block (in percent chunks?)
				while (r < points && !Vert::overlap(vert_alloc + r, vert_alloc + r - 1))
					vert_alloc[w++] = vert_alloc[r++];
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
			if (points == 2)
			{
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
			// todo: choose direction based on which side of contour it-h vertex appears at (so we won't have to reverse)
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

				/*
				// if i-th vert can see the contour we reverse entire base
				int h = i / 2;
				for (int j = 0; j < h; j++)
				{
					int k = i - 1 - j;
					Vert v = vert_alloc[j];
					vert_alloc[j] = vert_alloc[k];
					vert_alloc[k] = v;
				}
				*/
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
		if (vert_map)
			free(vert_map);
		vert_map = 0;

		inp_verts = points;
		out_verts = 0;

		first_dela_face = 0;
		first_hull_face = 0;
		first_boundary_vert = 0;

		if (max_verts < points)
		{
			if (max_verts)
			{
				//free(vert_alloc);
				delete [] vert_alloc;
				vert_alloc = 0;
				max_verts = 0;
			}

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

	virtual bool Constrain(int a, int b)
	{
		// current
		Vert* va = vert_alloc + a;// (Vert*)GetVertexByIndex(a);
		if (!va)
			return false;

		// destination
		Vert* vb = vert_alloc + b;//(Vert*)GetVertexByIndex(b);
		if (!vb)
			return false;

		DelaBella_Iterator it;
		const DelaBella_Triangle* f = va->StartIterator(&it);
		const DelaBella_Triangle* e = f;

		Face* face = 0;
		int edge;

		while (1)
		{
			if (f->index >= 0)
			{
				for (int e = 0; e < 3; e++)
				{
					if (f->v[e] == va)
					{
						// now heck if va-vb passes throu this face
						static const int other[3][2] = {{1,2},{2,0},{0,1}};
						Vert* v0 = (Vert*)(f->v[other[e][0]]);
						Vert* v1 = (Vert*)(f->v[other[e][1]]);
						IA_VAL ab[2] = { (IA_VAL)vb->x - (IA_VAL)va->x, (IA_VAL)vb->y - (IA_VAL)va->y };
						IA_VAL a0[2] = { (IA_VAL)v0->x - (IA_VAL)va->x, (IA_VAL)v0->y - (IA_VAL)va->y };
						IA_VAL a1[2] = { (IA_VAL)v1->x - (IA_VAL)va->x, (IA_VAL)v1->y - (IA_VAL)va->y };

						IA_VAL a0b = a0[0] * ab[1] - a0[1] * ab[0];
						IA_VAL a1b = a1[0] * ab[1] - a1[1] * ab[0];

						if (a0b < 0 && a1b > 0)
						{
							// clean
							face = (Face*)f;
							edge = e;
							break;
						}

						if (a0b > 0 || a1b < 0)
							break;

						// XA needed
						{
							int a = 0;
						}
					}
				}
				if (face)
					break;
			}
			f = it.Next();
			if (f == e)
				return false;
		}


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

	virtual const DelaBella_Vertex* GetVertexByIndex(int i)
	{
		if (i < 0 || i >= inp_verts)
			return 0;

		if (!vert_map)
		{
			vert_map = (int*)malloc(sizeof(int) * inp_verts);
			for (int j = 0; j < inp_verts; j++)
				vert_map[vert_alloc[j].i] = j;
		}

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

#ifndef __cplusplus
void* DelaBella_Create()
{
	return IDelaBella::Create();
}

void  DelaBella_Destroy(void* db)
{
	((IDelaBella*)db)->Destroy();
}

void  DelaBella_SetErrLog(void* db, int(*proc)(void* stream, const char* fmt, ...), void* stream)
{
	((IDelaBella*)db)->SetErrLog(proc, stream);
}

int   DelaBella_TriangulateFloat(void* db, int points, float* x, float* y, int advance_bytes)
{
	return ((IDelaBella*)db)->Triangulate(points, x, y, advance_bytes);
}

int   DelaBella_TriangulateDouble(void* db, int points, double* x, double* y, int advance_bytes)
{
	return ((IDelaBella*)db)->Triangulate(points, x, y, advance_bytes);
}

int   DelaBella_TriangulateLongDouble(void* db, int points, long double* x, long double* y, int advance_bytes)
{
	return ((IDelaBella*)db)->Triangulate(points, x, y, advance_bytes);
}


int   DelaBella_GetNumInputPoints(void* db)
{
	return ((IDelaBella*)db)->GetNumInputPoints();
}

int   DelaBella_GetNumOutputIndices(void* db)
{
	return ((IDelaBella*)db)->GetNumOutputIndices();
}

int   Delabella_GetNumOutputHullFaces(void* db);
{
	return ((IDelaBella*)db)->GetNumOutputHullFaces();
}


const DelaBella_Triangle* GetFirstDelaunayTriangle(void* db)
{
	return ((IDelaBella*)db)->GetFirstDelaunayTriangle();
}

const DelaBella_Triangle* GetFirstHullTriangle(void* db)
{
	return ((IDelaBella*)db)->GetFirstHullTriangle();
}

const DelaBella_Vertex*   GetFirstHullVertex(void* db)
{
	return ((IDelaBella*)db)->GetFirstHullVertex();
}
#endif
