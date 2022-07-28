/*
DELABELLA - Delaunay triangulation library
Copyright (C) 2018-2022 GUMIX - Marcin Sokalski
*/

#define DELABELLA_AUTOTEST
// in case of troubles, allows to see if any assert pops up.
// define it globally (like with -DDELABELLA_AUTOTEST)

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include "delabella.h"
#include "predicates.h"

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif

static const bool classify = true;

static uint64_t uSec()
{
#ifdef _WIN32
	LARGE_INTEGER c;
	static LARGE_INTEGER f;
	static BOOL bf = QueryPerformanceFrequency(&f);
	QueryPerformanceCounter(&c);
	uint64_t n = c.QuadPart;
	uint64_t d = f.QuadPart;
	uint64_t m = 1000000;
	// calc microseconds = n*m/d carefully!
	// naive mul/div would work only for upto 5h on 1GHz freq
	// we exploit fact that m*d fits in uint64 (upto 18THz freq)
	// so n%d*m fits as well,
	return n / d * m + n % d * m / d;
#else
	timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return (uint64_t)ts.tv_sec * 1000000 + ts.tv_nsec / 1000;
#endif
}


template <typename T, typename I>
IDelaBella2<T,I>::~IDelaBella2()
{
}

template <typename T, typename I>
struct CDelaBella2 : IDelaBella2<T,I>
{
	CDelaBella2() : 
		vert_map(0),
		vert_alloc(0),
		face_alloc(0),
		max_verts(0),
		max_faces(0),
		first_dela_face(0),
		first_hull_face(0),
		first_boundary_vert(0),
		first_internal_vert(0),
		inp_verts(0),
		out_verts(0),
		out_boundary_verts(0),
		polygons(0),
		out_hull_faces(0),
		unique_points(0),
		errlog_proc(0),
		errlog_file(0)
	{
	}

	struct Face;

	struct Iter : IDelaBella2<T,I>::Iterator {};

	struct Vert : IDelaBella2<T,I>::Vertex
	{
		static bool overlap(const Vert* v1, const Vert* v2)
		{
			return v1->x == v2->x && v1->y == v2->y;
		}

		bool operator < (const Vert& v) const
		{
			T dif = predicates::adaptive::sqrlendif2d(this->x, this->y, v.x, v.y);

			if (dif < 0)
				return true;
			if (dif > 0)
				return false;
			if (this->x < v.x || this->x == v.x && this->y < v.y)
				return true;
			return false;
		}
	};

	struct Face : IDelaBella2<T,I>::Simplex
	{
		void RotateEdgeFlagsCW() // <<
		{
			uint8_t f = this->flags;
			this->flags = ((f << 1) & 0b00110110) | ((f >> 2) & 0b00001001) | (f & 0b11000000);
		}

		void RotateEdgeFlagsCCW() // >>
		{
			uint8_t f = this->flags;
			this->flags = ((f >> 1) & 0b00011011) | ((f << 2) & 0b00100100) | (f & 0b11000000);
		}

		bool IsDelaFace() const
		{
			return !(this->flags & 0b10000000);
		}

		bool IsEdgeFixed(int at) const
		{
			return this->flags & (1 << at);
		}

		void ToggleEdgeFixed(int at)
		{
			this->flags = (this->flags ^ (0b00001000 << at)) | (0b00000001 << at);
		}

		uint8_t GetEdgeBits(int at) const
		{
			return (this->flags >> at) & 0b00001001;
		}

		void SetEdgeBits(int at, uint8_t bits)
		{
			this->flags = (bits << at) | (this->flags & ~(0b00001001 << at));
		}

		static Face* Alloc(Face** from)
		{
			Face* f = *from;
			*from = (Face*)f->next;
			f->next = 0;
			return f;
		}

		void Free(Face** to)
		{
			this->next = *to;
			*to = this;
		}

		Face* Next(const Vert* p) const
		{
			if (this->v[0] == p)
				return (Face*)this->f[1];
			if (this->v[1] == p)
				return (Face*)this->f[2];
			if (this->v[2] == p)
				return (Face*)this->f[0];
			return 0;
		}

		bool signN() const
		{
			return 0 > predicates::adaptive::orient2d(
				this->v[0]->x, this->v[0]->y, 
				this->v[1]->x, this->v[1]->y, 
				this->v[2]->x, this->v[2]->y);
		}

		bool sign0() const
		{
			return 0 == predicates::adaptive::orient2d(
				this->v[0]->x, this->v[0]->y, 
				this->v[1]->x, this->v[1]->y, 
				this->v[2]->x, this->v[2]->y);
		}

		bool dot0(const Vert& p) const
		{
			return
				predicates::adaptive::incircle(
					p.x, p.y,
					this->v[0]->x, this->v[0]->y,
					this->v[1]->x, this->v[1]->y,
					this->v[2]->x, this->v[2]->y) == 0;
		}

		bool dotP(const Vert& p) const
		{
			return
				predicates::adaptive::incircle(
					p.x, p.y,
					this->v[0]->x, this->v[0]->y,
					this->v[1]->x, this->v[1]->y,
					this->v[2]->x, this->v[2]->y) > 0;
		}

		bool dotN(const Vert& p) const
		{
			return
				predicates::adaptive::incircle(
					p.x, p.y,
					this->v[0]->x, this->v[0]->y,
					this->v[1]->x, this->v[1]->y,
					this->v[2]->x, this->v[2]->y) < 0;
		}

		bool dotNP(const Vert& p) const
		{
			return
				predicates::adaptive::incircle(
					p.x,p.y, 
					this->v[0]->x, this->v[0]->y, 
					this->v[1]->x, this->v[1]->y, 
					this->v[2]->x, this->v[2]->y) <= 0;
		}
	};

	Vert* vert_alloc;
	Face* face_alloc;
	I* vert_map;
	I max_verts;
	I max_faces;

	Face* first_dela_face;
	Face* first_hull_face;
	Vert* first_boundary_vert;
	Vert* first_internal_vert;

	I inp_verts;
	I out_verts;
	I polygons;
	I out_hull_faces;
	I out_boundary_verts;
	I unique_points;

	int(*errlog_proc)(void* file, const char* fmt, ...);
	void* errlog_file;

	I Prepare(I* start, Face** hull, I* out_hull_faces, Face** cache, uint64_t* sort_stamp)
	{
		uint64_t time0 = uSec();

		if (errlog_proc)
			errlog_proc(errlog_file, "[...] sorting vertices");
		I points = inp_verts;

		std::sort(vert_alloc, vert_alloc + points);

		// rmove dups
		{
			vert_map[vert_alloc[0].i] = 0;

			I w = 0, r = 1; // skip initial no-dups block
			while (r < points && !Vert::overlap(vert_alloc + r, vert_alloc + w))
			{
				vert_map[vert_alloc[r].i] = r;
				w++;
				r++;
			}

			I d = w; // dup map
			w++;

			while (r < points)
			{
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
				{
					vert_map[vert_alloc[r].i] = w;
					vert_alloc[w++] = vert_alloc[r++];
				}

				d = w - 1;
			}

			uint64_t time1 = uSec();
			if (sort_stamp)
				*sort_stamp = time1;

			if (errlog_proc)
				errlog_proc(errlog_file, "\r[100] sorting vertices (%lld ms)\n", (time1 - time0) / 1000);
			time0 = time1;

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
				vert_alloc[0].next = vert_alloc + 1;
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

		I i;
		Face f; // tmp
		f.v[0] = vert_alloc + 0;

		T lo_x = vert_alloc[0].x, hi_x = lo_x;
		T lo_y = vert_alloc[0].y, hi_y = lo_y;
		I lower_left = 0;
		I upper_right = 0;
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
		
		int pro = 0;
		// skip until points are coplanar
		while (i < points && f.dot0(vert_alloc[i]))
		{
			if (i >= pro)
			{
				uint64_t p = (int)((uint64_t)100 * i / points);
				pro = (int)((p + 1) * points / 100);
				if (pro >= points)
					pro = (int)points - 1;
				if (i == points - 1)
				{
					p = 100;
				}
				if (errlog_proc)
					errlog_proc(errlog_file, "\r[%2d%s] convex hull triangulation ", p, p >= 100 ? "" : "%");
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

		I* vert_sub = 0;

		bool colinear = f.sign0(); // hybrid
		if (colinear)
		{
			vert_sub = (I*)malloc(sizeof(I) * ((size_t)i+1));
			if (!vert_sub)
			{
				if (errlog_proc)
					errlog_proc(errlog_file, "[ERR] Not enough memory, shop for some more RAM. See you!\n");
				return 0;
			}

			for (I s = 0; s <= i; s++)
				vert_sub[s] = s;

			// choose x or y axis to sort verts (no need to be exact)
			if (hi_x - lo_x > hi_y - lo_y)
			{
				struct 
				{ 
					bool operator() (const I& a, const I& b)
					{
						return less(vert_alloc[a], vert_alloc[b]);
					}

					bool less(const Vert& a, const Vert& b) const 
					{ 
						return a.x < b.x; 
					}

					Vert* vert_alloc;
				} c;

				c.vert_alloc = vert_alloc;
				std::sort(vert_sub, vert_sub + i, c);
			}
			else
			{
				struct 
				{ 
					bool operator() (const I& a, const I& b)
					{
						return less(vert_alloc[a], vert_alloc[b]);
					}

					bool less(const Vert& a, const Vert& b) const 
					{ 
						return a.y < b.y; 
					} 

					Vert* vert_alloc;
				} c;

				c.vert_alloc = vert_alloc;
				std::sort(vert_sub, vert_sub + i, c);
			}
		}
		else
		{
			// split to lower (below diag) and uppert (above diag) parts, 
			// sort parts separately in opposite directions
			// mark part with Vert::sew temporarily

			for (I j = 0; j < i; j++)
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
					Vert* ll = vert_alloc + lower_left;
					Vert* ur = vert_alloc + upper_right;
					T dot = predicates::adaptive::orient2d(ll->x, ll->y, ur->x, ur->y, vert_alloc[j].x, vert_alloc[j].y);
					if (dot < 0)
					{
						// lower
						vert_alloc[j].sew = &f;
					}
					else
					{
						#ifdef DELABELLA_AUTOTEST
						assert(dot > 0);
						#endif
						// upper
						vert_alloc[j].sew = 0;
					}
				}
			}

			struct
			{
				// default to CW order (unlikely, if wrong, we will reverse using one=2 and two=1)
				
				bool operator()(const I& a, const I& b) const
				{
					return less(vert_alloc[a], vert_alloc[b]);
				}

				bool less(const Vert& a, const Vert& b) const
				{
					// if a is lower and b is upper, return true
					if (a.sew && !b.sew)
						return false;

					// if a is upper and b is lower, return false
					if (!a.sew && b.sew)
						return true;

					// actually we can compare coords directly
					if (a.sew)
					{
						// lower
						if (a.x > b.x)
							return true;
						if (a.x == b.x)
							return a.y > b.y;
						return false;
					}
					else
					{
						// upper
						if (a.x < b.x)
							return true;
						if (a.x == b.x)
							return a.y < b.y;
						return false;
					}

					// otherwise
					#ifdef DELABELLA_AUTOTEST
					assert(0);
					#endif
					return false;
				}

				Vert* vert_alloc;
			} c;

			vert_sub = (I*)malloc(sizeof(I) * ((size_t)i+1));
			if (!vert_sub)
			{
				if (errlog_proc)
					errlog_proc(errlog_file, "[ERR] Not enough memory, shop for some more RAM. See you!\n");
				return 0;
			}

			for (I s = 0; s <= i; s++)
				vert_sub[s] = s;

			c.vert_alloc = vert_alloc;
			std::sort(vert_sub, vert_sub + i, c);
		}

		*out_hull_faces = 0;

		// alloc faces only if we're going to create them
		if (i < points || !colinear)
		{
			I hull_faces = 2 * points - 4;
			*out_hull_faces = hull_faces;

			if (max_faces < hull_faces)
			{
				if (max_faces)
				{
					free(face_alloc);
					//delete[] face_alloc;
				}
				max_faces = 0;

				/*
				try 
				{
					face_alloc = new Face[(size_t)hull_faces];
				}
				catch (...)
				{
					face_alloc = 0;
				}
				*/

				face_alloc = (Face*)malloc(sizeof(Face) * hull_faces);

				if (face_alloc)
					max_faces = hull_faces;
				else
				{
					if (errlog_proc)
						errlog_proc(errlog_file, "[ERR] Not enough memory, shop for some more RAM. See you!\n");
					return 0;
				}
			}

			face_alloc[0].flags = 0;
			for (I i = 1; i < hull_faces; i++)
			{
				face_alloc[i - 1].next = face_alloc + i;
				face_alloc[i].flags = 0;
			}
			face_alloc[hull_faces - 1].next = 0;

			*cache = face_alloc;
		}

		if (i == points)
		{
			if (errlog_proc)
			{
				// TODO:
				// this should be moved to end of Triangulate()
				// after printing final "[100] ... \n"
				if (colinear)
					errlog_proc(errlog_file, "[WRN] all input points are colinear\n");
				else
					errlog_proc(errlog_file, "[WRN] all input points are cocircular\n");
			}

			if (colinear)
			{
				// link verts into open list
				first_boundary_vert = vert_alloc + vert_sub[0];
				first_internal_vert = 0;
				out_boundary_verts = points;

				for (I j = 1; j < points; j++)
					vert_alloc[j - 1].next = vert_alloc + vert_sub[j];
				vert_alloc[points - 1].next = 0;

				free(vert_sub);
				return -points;
			}

			// we're almost done!
			// time to build final flat hull
			Face* next_p = Face::Alloc(cache);
			Face* next_q = Face::Alloc(cache);
			Face* prev_p = 0;
			Face* prev_q = 0;

			for (I j = 2; j < i; j++)
			{
				Face* p = next_p;
				p->v[0] = vert_alloc + vert_sub[0];
				p->v[1] = vert_alloc + vert_sub[j-1];
				p->v[2] = vert_alloc + vert_sub[j];

				// mirrored
				Face* q = next_q;
				q->v[0] = vert_alloc + vert_sub[0];
				q->v[1] = vert_alloc + vert_sub[j];
				q->v[2] = vert_alloc + vert_sub[j-1];

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

			f.v[0] = vert_alloc + vert_sub[0];
			f.v[1] = vert_alloc + vert_sub[1];
			f.v[2] = vert_alloc + vert_sub[2];

			int one = 1, two = 2;
			if (!f.dotNP(vert_alloc[vert_sub[i]]))
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

			for (I j = 2; j < i; j++)
			{
				Face* p = next_p;
				p->v[0] = vert_alloc + vert_sub[0];
				p->v[one] = vert_alloc + vert_sub[j - 1];
				p->v[two] = vert_alloc + vert_sub[j];

				Face* q;

				if (j == 2)
				{
					// first base triangle also build extra tip face
					q = Face::Alloc(cache);
					q->v[0] = vert_alloc + vert_sub[i];
					q->v[one] = vert_alloc + vert_sub[1];
					q->v[two] = vert_alloc + vert_sub[0];

					q->f[0] = p;
					q->f[one] = 0; // LAST_Q;
					q->f[two] = next_q;

					first_q = q;
					prev_q = q;
				}

				q = next_q;
				q->v[0] = vert_alloc + vert_sub[i];
				q->v[one] = vert_alloc + vert_sub[j];
				q->v[two] = vert_alloc + vert_sub[j-1];

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
					q->v[0] = vert_alloc + vert_sub[i];
					q->v[one] = vert_alloc + vert_sub[0];
					q->v[two] = vert_alloc + vert_sub[i-1];

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
		}

		free(vert_sub);

		*start = i;
		return points;
	}

	I Triangulate(I* other_faces)
	{
		I i = 0;
		Face* hull = 0;
		I hull_faces = 0;
		Face* cache = 0;

		uint64_t sort_stamp;
		I points = Prepare(&i, &hull, &hull_faces, &cache, &sort_stamp);
		unique_points = points < 0 ? -points : points;
		if (points <= 0)
		{
			return points;
		}

		/////////////////////////////////////////////////////////////////////////
		// ACTUAL ALGORITHM

		int pro = 0;
		for (; i < points; i++)
		{
			if (i >= pro)
			{
				uint64_t p = (int)((uint64_t)100 * i / points);
				pro = (int)((p+1) * points / 100);
				if (pro >= points)
					pro = (int)points - 1;
				if (i == points - 1)
					p = 100;
				if (errlog_proc)
					errlog_proc(errlog_file, "\r[%2d%s] convex hull triangulation ", p, p>=100 ? "":"%");
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
					f = face_alloc + (intptr_t)2 * i - 4 - 1;
					//while (f->dot(*q) <= 0)
					while (f->dotNP(*q))
					//while (f->dotN(*q)) // we want to consume coplanar faces
					{
						#ifdef DELABELLA_AUTOTEST
						assert(f != face_alloc); // no face is visible? you must be kidding!
						#endif
						f--;
					}
				}
			}

			// 2. DELETE VISIBLE FACES & ADD NEW ONES
			//    (we also build silhouette (vertex loop) between visible & invisible faces)

			I del = 0;
			I add = 0;

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
							#ifdef DELABELLA_AUTOTEST
							else
								assert(0);
							#endif

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
							#ifdef DELABELLA_AUTOTEST
							else
								assert(0);
							#endif

							// push neighbor face, it's visible and requires processing
							n->next = stack ? stack : n;
							stack = n;
						}
					}
				}
			}

			#ifdef DELABELLA_AUTOTEST
			// if add<del+2 hungry hull has consumed some point
			// that means we can't do delaunay for some under precission reasons
			// althought convex hull would be fine with it
			assert(add == del + 2);
			#endif

			// 3. SEW SIDES OF CONE BUILT ON SLIHOUTTE SEGMENTS

			hull = face_alloc + (intptr_t)2 * i - 4 + 1; // last added face

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

		#ifdef DELABELLA_AUTOTEST
		assert(2 * i - 4 == hull_faces);
		#endif

		for (I j = 0; j < points; j++)
		{
			vert_alloc[j].next = 0;
			vert_alloc[j].sew = 0;
		}

		I others = 0;

		i = 0;
		Face** prev_dela = &first_dela_face;
		Face** prev_hull = &first_hull_face;
		for (I j = 0; j < hull_faces; j++)
		{
			Face* f = face_alloc + j;

			// back-link all verts to some_face
			// yea, ~6x times, sorry
			((Vert*)f->v[0])->sew = f;
			((Vert*)f->v[1])->sew = f;
			((Vert*)f->v[2])->sew = f;

			if (f->signN())
			{
				f->index = i;  // store index in dela list
				f->flags = 0b01000000; // interior
				*prev_dela = f;
				prev_dela = (Face**)&f->next;
				i++;
			}
			else
			{
				f->index = others; // store index in hull list (~ to mark it is not dela)
				f->flags = 0b10000000; // hull
				*prev_hull = f;
				prev_hull = (Face**)&f->next;
				others++;
			}
		}

		if (other_faces)
			*other_faces = others;

		*prev_dela = 0;
		*prev_hull = 0;

		// let's trace boudary contour, at least one vertex of first_hull_face
		// must be shared with dela face, find that dela face
		Iter it;
		Vert* v = (Vert*)first_hull_face->v[0];
		Face* t = (Face*)v->StartIterator(&it);
		Face* e = t; // end
		
		first_boundary_vert = (Vert*)v;
		out_boundary_verts = 1;

		while (1)
		{
			if (t->IsDelaFace())
			{
				int pr = it.around-1; if (pr<0) pr = 2;
				int nx = it.around+1; if (nx>2) nx = 0;

				Face* fpr = (Face*)t->f[pr];
				if (!fpr->IsDelaFace())
				{
					// let's move from: v to t->v[nx]
					v->next = t->v[nx];
					v = (Vert*)v->next;
					if (v == first_boundary_vert)
						break; // lap finished
					out_boundary_verts++;
					t = (Face*)t->StartIterator(&it,nx);
					e = t;
					continue;
				}
			}
			t = (Face*)it.Next();

			#ifdef DELABELLA_AUTOTEST
			assert(t!=e);
			#endif
		}

		// link all other verts into internal list
		first_internal_vert = 0;
		Vert** prev_inter = &first_internal_vert;
		for (I j=0; j<points; j++)
		{
			if (!vert_alloc[j].next)
			{
				Vert* next = vert_alloc+j;
				*prev_inter = next;
				prev_inter = (Vert**)&next->next;
			}
		}

		if (errlog_proc)
			errlog_proc(errlog_file, "\r[100] convex hull triangulation (%lld ms)\n", (uSec() - sort_stamp) / 1000);

		return 3*i;
	}

	bool ReallocVerts(I points)
	{
		inp_verts = points;
		out_verts = 0;
		polygons = 0;

		first_dela_face = 0;
		first_hull_face = 0;
		first_boundary_vert = 0;

		if (max_verts < points)
		{
			if (max_verts)
			{
				free(vert_map);
				vert_map = 0;

				free(vert_alloc);
				//delete [] vert_alloc;
				vert_alloc = 0;
				max_verts = 0;
			}

			/*
			try
			{
				vert_alloc = new Vert[(size_t)points];
			}
			catch (...)
			{
				vert_alloc = 0;
			}
			*/

			vert_alloc = (Vert*)malloc(sizeof(Vert) * points);

			if (vert_alloc)
				vert_map = (I*)malloc(sizeof(I) * (size_t)points);

			if (vert_alloc && vert_map)
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

		Iter it;

		Face* first = (Face*)va->StartIterator(&it);
		Face* face = first;

		// find first face around va, containing offending edge
		Vert* v0;
		Vert* v1;
		Face* N = 0;
		int a, b, c;

		while (1)
		{
			if (!face->IsDelaFace())
			{
				face = (Face*)it.Next();
				#ifdef DELABELLA_AUTOTEST
				assert(face != first);
				#endif
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

			T a0b = predicates::adaptive::orient2d(va->x,va->y, v0->x,v0->y, vb->x,vb->y);
			T a1b = predicates::adaptive::orient2d(va->x,va->y, v1->x,v1->y, vb->x,vb->y);
			
			if (a0b <= 0 && a1b >= 0)
			{
				// note: 
				// check co-linearity only if v0,v1 are pointing
				// to the right direction (from va to vb)

				if (a0b == 0)
				{
					*restart = v0;
					*tail = 0;
					*ptail = list ? tail : 0;
					return list;
				}

				if (a1b == 0)
				{
					*restart = v1;
					*tail = 0;
					*ptail = list ? tail : 0;
					return list;
				}

				// offending edge!
				N = (Face*)face;
				break;
			}

			face = (Face*)it.Next();

			#ifdef DELABELLA_AUTOTEST
			assert(face != first);
			#endif
		}

		while (1)
		{
			if (a)
			{
				// rotate N->v[] and N->f 'a' times 'backward' such offending edge appears opposite to v[0]
				const int* r = rotate[a];

				Vert* v[3] = { (Vert*)N->v[0], (Vert*)N->v[1], (Vert*)N->v[2] };
				N->v[0] = v[r[0]];
				N->v[1] = v[r[1]];
				N->v[2] = v[r[2]];

				Face* f[3] = { (Face*)N->f[0], (Face*)N->f[1], (Face*)N->f[2] };
				N->f[0] = f[r[0]];
				N->f[1] = f[r[1]];
				N->f[2] = f[r[2]];

				if (classify)
				{
					// rotate face_bits!
					if (a == 1)
						N->RotateEdgeFlagsCW();
					else // a==2
						N->RotateEdgeFlagsCCW();
				}
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
			T abr = predicates::adaptive::orient2d(va->x, va->y, vb->x, vb->y, vr->x, vr->y);
			
			if (abr == 0)
			{
				*restart = vr;
				*tail = 0;
				*ptail = list ? tail : 0;
				return list;
			}

			if (abr > 0)
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

		#ifdef DELABELLA_AUTOTEST
		assert(0);
		#endif
		*restart = 0;
		*ptail = 0;
		return 0;
	}

	// private:
	bool FindEdgeFaces(const Vert* v1, const Vert* v2, Face* twin[2], int opposed[2])
	{
		static const int other[3] = { 2,0,1 };

		Iter it;
		Face* f = (Face*)v1->StartIterator(&it);
		Face* e = f;
		do
		{
			if (f->v[other[it.around]] == v2)
			{
				// at least one of them must be dela, it will appear as twin[0]
				if (f->IsDelaFace())
				{
					opposed[0] = other[it.around];
					twin[0] = f;
					twin[1] = (Face*)it.Next();
					opposed[1] = other[it.around];
					return true;
				}

				int tmp_op = other[it.around];
				Face* tmp_f = f;
				f = (Face*)it.Next();

				if (f->IsDelaFace())
				{
					opposed[0] = other[it.around];
					twin[0] = f;
					opposed[1] = tmp_op;
					twin[1] = tmp_f;
					return true;
				}
			}
			else
				f = (Face*)it.Next();
		} while (f != e);

		return false;
	}

	virtual I ConstrainEdges(I edges, const I* pa, const I* pb, size_t advance_bytes)
	{
		if (advance_bytes == 0)
			advance_bytes = 2*sizeof(I);

		bool unfixed = false; // wrn?

		uint64_t time0 = uSec();

		int pro = 0;
		for (I con = 0; con < edges; con++)
		{
			if (con >= pro)
			{
				uint64_t p = (int)((uint64_t)100 * con / edges);
				pro = (int)((p + 1) * edges / 100);
				if (pro >= edges)
					pro = (int)edges - 1;
				if (con == edges - 1)
					p = 100;
				if (errlog_proc)
					errlog_proc(errlog_file, "\r[%2d%s] constraining ", p, p >= 100 ? "" : "%");
			}

			I a = *(const I*)((const char*)pa + (intptr_t)con * advance_bytes);
			I b = *(const I*)((const char*)pb + (intptr_t)con * advance_bytes);

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

			if (!list)
			{
				if (restart)
				{
					if (classify)
					{
						Face* nf[2];
						int op[2];

						#ifdef DELABELLA_AUTOTEST
						assert(FindEdgeFaces(va, restart, nf, op));
						#else
						FindEdgeFaces(va, restart, nf, op);
						#endif

						nf[0]->ToggleEdgeFixed(op[0]);
						nf[1]->ToggleEdgeFixed(op[1]);
					}

					va = restart;
					continue;
				}

				if (classify)
				{
					Face* nf[2];
					int op[2];

					#ifdef DELABELLA_AUTOTEST
					assert(FindEdgeFaces(va, vc, nf, op));
					#else
					FindEdgeFaces(va, vc, nf, op);
					#endif

					nf[0]->ToggleEdgeFixed(op[0]);
					nf[1]->ToggleEdgeFixed(op[1]);
				}
			}

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

				// is v,v0,vr,v1 a convex quad?
				T v0r = predicates::adaptive::orient2d(v->x, v->y, v0->x, v0->y, vr->x, vr->y);
				T v1r = predicates::adaptive::orient2d(v->x, v->y, v1->x, v1->y, vr->x, vr->y);
				
				if (v0r >= 0 || v1r <= 0)
				{
					// CONCAVE CUNT!
					*tail = N;
					tail = (Face**)&N->next;
					continue;
				}

				#ifdef DELABELLA_AUTOTEST
				assert(v0r < 0 && v1r > 0);
				#endif

				if (classify)
				{
					if (N->IsEdgeFixed(a))
					{
						unfixed = true;
						#ifdef DELABELLA_AUTOTEST
						assert(F->IsEdgeFixed(d));
					}
					else
					{
						assert(!F->IsEdgeFixed(d));
						#endif
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

						if (classify)
						{
							// TRANSFER bits:
							// N->b -> F->d
							// N->c -> N->b
							// F->e -> N->c

							F->SetEdgeBits(d, N->GetEdgeBits(b));
							N->SetEdgeBits(b, N->GetEdgeBits(c));
							N->SetEdgeBits(c, F->GetEdgeBits(e));

							if (v == va && vr == vb)
							{
								// set N->a & F->e
								N->SetEdgeBits(a, 0b00001001);
								F->SetEdgeBits(e, 0b00001001);
							}
							else
							{
								// clear N->a & F->e
								N->SetEdgeBits(a, 0);
								F->SetEdgeBits(e, 0);
							}
						}
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

						if (classify)
						{
							// TRANSFER bits:
							// N->c -> F->d
							// N->b -> N->c
							// F->f -> N->b
							// F->e -> F->e

							F->SetEdgeBits(d, N->GetEdgeBits(c));
							N->SetEdgeBits(c, N->GetEdgeBits(b));
							N->SetEdgeBits(b, F->GetEdgeBits(f));

							if (v == va && vr == vb)
							{
								// set N->a & F->f
								N->SetEdgeBits(a, 0b00001001);
								F->SetEdgeBits(f, 0b00001001);
							}
							else
							{
								N->SetEdgeBits(a, 0);
								F->SetEdgeBits(f, 0);
							}
						}
					}

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
						// check if v and vr are on the same side of a--b
						T abv = predicates::adaptive::orient2d(va->x, va->y, vb->x, vb->y, v->x, v->y);
						T abr = predicates::adaptive::orient2d(va->x, va->y, vb->x, vb->y, vr->x, vr->y);
						
						if (abv >= 0 && abr >= 0 || abv <= 0 && abr <= 0)
						{
							// resolved
							N->next = flipped;
							flipped = N;
						}
						else
						{
							// unresolved
							*tail = N;
							tail = (Face**)&N->next;
						}
					}
				}
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

					// can we check if it was flipped last time?
					// if ((N->index & 0x40000000) == 0)
					if (N->dotP(*vr) /* || F->dotP(*v)*/) 
					{
						if (classify)
						{
							// check if N->a or F->d is marked as fixed
							// it should not ever happen -- i think
							
							// ... but if it does happen
							// remember we've broken fixed edge
							// then display warning right before we return

							if (N->IsEdgeFixed(a))
							{
								printf("IT DOES HAPPEN!\n");
								unfixed = true;
								#ifdef DELABELLA_AUTOTEST
								assert(F->IsEdgeFixed(d));
							}
							else
							{
								assert(!F->IsEdgeFixed(d));
								#endif
							}
						}

						no_flips = false;

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

							if (classify)
							{
								// TRANSFER bits:
								// N->b -> F->d
								// N->c -> N->b
								// F->f -> F->f
								// F->e -> N->c

								F->SetEdgeBits(d, N->GetEdgeBits(b));
								N->SetEdgeBits(b, N->GetEdgeBits(c));
								N->SetEdgeBits(c, F->GetEdgeBits(e));

								// clear N->a & F->e
								N->SetEdgeBits(a, 0);
								F->SetEdgeBits(e, 0);
							}
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

							if (classify)
							{
								// TRANSFER bits:
								// N->c -> F->d
								// N->b -> N->c
								// F->f -> N->b
								// F->e -> F->e

								F->SetEdgeBits(d, N->GetEdgeBits(c));
								N->SetEdgeBits(c, N->GetEdgeBits(b));
								N->SetEdgeBits(b, F->GetEdgeBits(f));

								// clear N->a & F->f
								N->SetEdgeBits(a, 0);
								F->SetEdgeBits(f, 0);
							}
						}

						// can we un-mark not flipped somehow?
						// N->index &= 0x3fffffff;
						// F->index &= 0x3fffffff;
					}
					else
					{
						// can we mark it as not flipped somehow?
						// N->index |= 0x40000000;
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
		I hull_faces = 2 * unique_points - 4;
		Face** tail = &first_dela_face;
		I index = 0;
		for (I i = 0; i < hull_faces; i++)
		{
			if (face_alloc[i].IsDelaFace())
			{
				*tail = face_alloc + i;
				tail = (Face**)&face_alloc[i].next;
				face_alloc[i].index = index++;
			}
		}
		*tail = 0;

		polygons = index;

		if (errlog_proc)
			errlog_proc(errlog_file, "\r[100] constraining (%lld ms)\n", (uSec() - time0) / 1000);

		return polygons;
	}

	virtual I Polygonize(const typename IDelaBella2<T,I>::Simplex* poly[])
	{
		const I marker = (I)(-1);
		

		uint64_t time0 = uSec();
		Face** buf = 0;
		if (!poly)
		{
			buf = (Face**)malloc(sizeof(Face*) * (size_t)out_verts / 3);
			poly = (const typename IDelaBella2<T,I>::Simplex**)buf;
			if (!poly)
				return -1;
		}

		// clear poly indices;
		Face* f = first_dela_face;
		while (f)
		{
			f->index = marker;
			f = (Face*)f->next;
		}

		I num = 0;
		f = first_dela_face;
		int pro = 0;
		I faces = out_verts / 3;
		I i = 0;
		while (f)
		{
			if (i >= pro)
			{
				uint64_t p = (int)((uint64_t)100 * i / faces);
				pro = (int)((p + 1) * faces / 100);
				if (pro >= faces)
					pro = (int)faces - 1;
				if (i == faces - 1)
					p = 100;
				if (errlog_proc)
					errlog_proc(errlog_file, "\r[%2d%s] polygonizing ", p, p >= 100 ? "" : "%");
			}
			i++;

			bool new_poly = true;
			Face* next = (Face*)f->next;
			for (int i = 0; i < 3; i++)
			{
				if (classify && f->IsEdgeFixed(i))
					continue;

				Face* a = (Face*)f->f[i];

				I index = a->index;

				if (index != marker && a->IsDelaFace())
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
						I dest = f->index;
						if (dest != marker)
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

		polygons = num;

		#if 1
		// ALTER POST PROC:
		// re-order triangles in every polygon and re-order indices in faces:
		// - first face defines first 3 verts in contour with v[0],v[1],v[2]
		// - every next face define just 1 additional vertex at its v[0]
		// thay can form fan or strip or arbitraty mix of both
		// without changing existing edges!

		for (I p = num-1; p >= 0; p--)
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
			}

			#ifdef DELABELLA_AUTOTEST
			assert(first);
			#endif

			Face* list = first; // here list is empty, first is used as sentinel

			f = first; // current face

			Face* last = 0; // will be one inserted right after first
			
			bool step_on = false; // is current vertex inserted

			Iter it;
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

							if (classify)
								f->RotateEdgeFlagsCW();
						}
						else // it.around == 2
						{
							f->f[0] = f->f[2];
							f->f[2] = f->f[1];
							f->f[1] = fr;
							f->v[0] = f->v[2];
							f->v[2] = f->v[1];
							f->v[1] = vr;

							if (classify)
								f->RotateEdgeFlagsCCW();
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

		if (errlog_proc)
			errlog_proc(errlog_file, "\r[100] polygonizing (%lld ms)\n", (uSec() - time0) / 1000);

		return num;
	}

	virtual I Triangulate(I points, const T* x, const T* y, size_t advance_bytes)
	{
		// const size_t max_triangulate_indices = (size_t)points * 6 - 15;
		// const size_t max_voronoi_edge_indices = (size_t)points * 6 - 12;
		const size_t max_voronoi_poly_indices = (size_t)points * 7 - 9; // winner of shame!

		// upcoming x16 - 40 (every triangle with 3 edge flags)

		if (max_voronoi_poly_indices > (size_t)std::numeric_limits<I>::max())
		{
			if (errlog_proc)
				errlog_proc(errlog_file, "[ERR] index type too small for provided number of points!\n");
			return 0;
		}

		if (!x)
			return 0;
		
		if (!y)
			y = x + 1;
		
		if (advance_bytes < sizeof(T) * 2)
			advance_bytes = sizeof(T) * 2;

		if (!ReallocVerts(points))
			return 0;

		for (I i = 0; i < points; i++)
		{
			Vert* v = vert_alloc + i;
			v->i = i;
			v->x = *(const T*)((const char*)x + i*advance_bytes);
			v->y = *(const T*)((const char*)y + i*advance_bytes);
		}

		out_hull_faces = 0;
		unique_points = 0;
		out_verts = Triangulate(&out_hull_faces);
		polygons = out_verts / 3;
		return out_verts;
	}

	virtual void Destroy()
	{
		if (vert_map)
			free(vert_map);

		if (face_alloc)
		{
			//delete [] face_alloc;
			free(face_alloc);
		}

		if (vert_alloc)
		{
			//delete [] vert_alloc;
			free(vert_alloc);
		}

		delete this;
	}

	// num of points passed to last call to Triangulate()
	virtual I GetNumInputPoints() const
	{
		return inp_verts;
	}

	// num of verts returned from last call to Triangulate()
	virtual I GetNumOutputIndices() const
	{
		return out_verts;
	}

	virtual I GetNumOutputHullFaces() const
	{
		return out_hull_faces;
	}

	virtual I GetNumBoundaryVerts() const
	{
		return out_verts < 0 ? -out_verts : out_boundary_verts;
	}

	virtual I GetNumInternalVerts() const
	{
		return out_verts < 0 ? 0 : unique_points - out_boundary_verts;
	}

	// num of polygons
	virtual I GetNumPolygons() const
	{
		return polygons;
	}

	virtual const typename IDelaBella2<T,I>::Simplex* GetFirstDelaunaySimplex() const
	{
		return first_dela_face;
	}

	virtual const typename IDelaBella2<T,I>::Simplex* GetFirstHullSimplex() const
	{
		return first_hull_face;
	}

	virtual const typename IDelaBella2<T,I>::Vertex* GetFirstBoundaryVertex() const
	{
		return first_boundary_vert;
	}

	virtual const typename IDelaBella2<T,I>::Vertex* GetFirstInternalVertex() const
	{
		return first_internal_vert;
	}

	virtual const typename IDelaBella2<T,I>::Vertex* GetVertexByIndex(I i) const
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

	virtual I GenVoronoiDiagramVerts(T* x, T* y, size_t advance_bytes) const
	{
		if (!first_dela_face)
			return 0;

		const I polys = polygons;
		const I contour = out_boundary_verts;
		I ret = polys + contour;

		if (!x || !y)
			return ret;

		if (advance_bytes < sizeof(T) * 2)
			advance_bytes = sizeof(T) * 2;

		const Face* f = first_dela_face;
		while (f)
		{
			const T v1x = f->v[1]->x - f->v[0]->x;
			const T v1y = f->v[1]->y - f->v[0]->y;
			const T v2x = f->v[2]->x - f->v[0]->x;
			const T v2y = f->v[2]->y - f->v[0]->y;

			const T v11 = v1x * v1x + v1y * v1y;
			const T v22 = v2x * v2x + v2y * v2y;
			const T v12 = (v1x * v2y - v1y * v2x) * 2;

			T cx = f->v[0]->x + (v2y * v11 - v1y * v22) / v12;
			T cy = f->v[0]->y + (v1x * v22 - v2x * v11) / v12;

			// yes, for polys < tris
			// we possibly calc it multiple times
			// and overwrite already calculated centers
			size_t offs = advance_bytes * (size_t)f->index;
			*(T*)((char*)x + offs) = cx;
			*(T*)((char*)y + offs) = cy;

			f = (Face*)f->next;
		}

		{
			size_t offs = advance_bytes * (size_t)polys;
			x = (T*)((char*)x + offs);
			y = (T*)((char*)y + offs);
		}

		Vert* prev = first_boundary_vert;
		Vert* vert = (Vert*)prev->next;
		for (I i = 0; i < contour; i++)
		{
			T nx = prev->y - vert->y;
			T ny = vert->x - prev->x;

			T nn = 1/sqrt(nx * nx + ny * ny);
			nx *= nn;
			ny *= nn;

			size_t offs = advance_bytes * (size_t)i;
			*(T*)((char*)x + offs) = nx;
			*(T*)((char*)y + offs) = ny;

			prev = vert;
			vert = (Vert*)vert->next;
		}

		return ret;
	}

	virtual I GenVoronoiDiagramEdges(I* indices, size_t advance_bytes) const
	{
		if (!first_dela_face)
			return 0;

		const I polys = polygons;
		const I verts = unique_points;
		I ret = 2 * (verts + polys - 1);

		if (!indices)
			return ret;

		if (advance_bytes < sizeof(I))
			advance_bytes = sizeof(I);

		const I contour = out_boundary_verts;
		const I inter = verts - contour;

		I* idx = indices;

		Vert* vert = first_internal_vert;
		for (I i = 0; i < inter; i++)
		{
			Iter it;
			Face* t = (Face*)vert->StartIterator(&it);
			Face* e = t;

			I a = t->index; // begin
			do
			{
				I b = t->index;
				if (a < b)
				{
					*idx = a;
					idx = (I*)((char*)idx + advance_bytes);
					*idx = b;
					idx = (I*)((char*)idx + advance_bytes);
				}
				a = b;
				t = (Face*)it.Next();
			} while (t != e);

			// loop closing seg
			I b = t->index;
			if (a < b)
			{
				*idx = a;
				idx = (I*)((char*)idx + advance_bytes);
				*idx = b;
				idx = (I*)((char*)idx + advance_bytes);
			}
			a = b;

			vert = (Vert*)vert->next;
		}

		Vert* prev = first_boundary_vert;
		vert = (Vert*)prev->next;
		for (I i = 0; i < contour; i++)
		{
			I a = i + polys; // begin

			// iterate all dela faces around prev
			// add their voro-vert index == dela face index
			Iter it;
			Face* t = (Face*)prev->StartIterator(&it);

			// it starts at random face, so lookup the prev->vert edge
			while (1)
			{
				if (t->IsDelaFace())
				{
					if (t->v[0] == prev && t->v[1] == vert ||
						t->v[1] == prev && t->v[2] == vert ||
						t->v[2] == prev && t->v[0] == vert)
						break;
				}
				t = (Face*)it.Next();
			}

			// now iterate around, till we're inside the boundary
			while (t->IsDelaFace())
			{
				I b = t->index;

				if (a<b)
				{
					*idx = a;
					idx = (I*)((char*)idx + advance_bytes);
					*idx = b;
					idx = (I*)((char*)idx + advance_bytes);
				}
				a = b;

				t = (Face*)it.Next();
			}

			I b = (i == 0 ? contour - 1 : i - 1) + polys; // loop-wrapping!
			if (a < b)
			{
				*idx = a;
				idx = (I*)((char*)idx + advance_bytes);
				*idx = b;
				idx = (I*)((char*)idx + advance_bytes);
			}
			a = b;

			prev = vert;
			vert = (Vert*)vert->next;
		}

		#ifdef DELABELLA_AUTOTEST
		assert(((char*)idx-(char*)indices) / advance_bytes == ret);
		#endif

		return ret;
	}

	virtual I GenVoronoiDiagramPolys(I* indices, size_t advance_bytes, I* closed_indices) const
	{
		if (!first_dela_face)
			return 0;

		const I polys = polygons;
		const I contour = out_boundary_verts;
		const I verts = unique_points;
		I ret = 3 * verts + 2 * (polys - 1) + contour;

		if (!indices)
			return ret;

		if (advance_bytes < sizeof(I))
			advance_bytes = sizeof(I);

		const I inter = verts - contour;
		const I poly_ending = ~0;

		I* idx = indices;
		
		I closed = 0;
		Vert* vert = first_internal_vert;
		for (I i = 0; i < inter; i++)
		{
			Iter it;
			Face* t = (Face*)vert->StartIterator(&it);
			Face* e = t;

			do
			{
				I a = t->index;
				*idx = a;
				idx = (I*)((char*)idx + advance_bytes);
				closed++;

				t = (Face*)it.Next();
			} while (t != e);

			// loop closing seg
			I a = poly_ending;
			*idx = a;
			idx = (I*)((char*)idx + advance_bytes);
			closed++;

			vert = (Vert*)vert->next;
		}

		if (closed_indices)
			*closed_indices = closed;

		Vert* prev = first_boundary_vert;
		vert = (Vert*)prev->next;
		for (I i = 0; i < contour; i++)
		{
			I a = i + polys; // begin
			*idx = a;
			idx = (I*)((char*)idx + advance_bytes);

			// iterate all dela faces around prev
			// add their voro-vert index == dela face index
			Iter it;
			Face* t = (Face*)prev->StartIterator(&it);

			// it starts at random face, so lookup the prev->vert edge
			while (1)
			{
				if (t->IsDelaFace())
				{
					if (t->v[0] == prev && t->v[1] == vert ||
						t->v[1] == prev && t->v[2] == vert ||
						t->v[2] == prev && t->v[0] == vert)
						break;
				}
				t = (Face*)it.Next();
			}

			// now iterate around, till we're inside the boundary
			while (t->IsDelaFace())
			{
				I a = t->index;
				*idx = a;
				idx = (I*)((char*)idx + advance_bytes);

				t = (Face*)it.Next();
			}

			a = (i == 0 ? contour - 1 : i - 1) + polys; // loop-wrapping!
			*idx = a;
			idx = (I*)((char*)idx + advance_bytes);

			a = poly_ending;
			*idx = a;
			idx = (I*)((char*)idx + advance_bytes);

			prev = vert;
			vert = (Vert*)vert->next;
		}

		#ifdef DELABELLA_AUTOTEST
		assert(((char*)idx - (char*)indices) / advance_bytes == ret);
		#endif

		return verts;
	}
};

template <typename T, typename I>
IDelaBella2<T,I>* IDelaBella2<T,I>::Create()
{
	IDelaBella2<T, I>* ret = 0;
	try
	{ 
		ret = new CDelaBella2<T,I>;
	}
	catch (...)
	{
		ret = 0;
	}
	return ret;
}

// this should cover all malcontents
template IDelaBella2<float, int8_t>* IDelaBella2<float, int8_t>::Create();
template IDelaBella2<double, int8_t>* IDelaBella2<double, int8_t>::Create();
template IDelaBella2<long double, int8_t>* IDelaBella2<long double, int8_t>::Create();

template IDelaBella2<float, int16_t>* IDelaBella2<float, int16_t>::Create();
template IDelaBella2<double, int16_t>* IDelaBella2<double, int16_t>::Create();
template IDelaBella2<long double, int16_t>* IDelaBella2<long double, int16_t>::Create();

template IDelaBella2<float, int32_t>* IDelaBella2<float, int32_t>::Create();
template IDelaBella2<double, int32_t>* IDelaBella2<double, int32_t>::Create();
template IDelaBella2<long double, int32_t>* IDelaBella2<long double, int32_t>::Create();

template IDelaBella2<float, int64_t>* IDelaBella2<float, int64_t>::Create();
template IDelaBella2<double, int64_t>* IDelaBella2<double, int64_t>::Create();
template IDelaBella2<long double, int64_t>* IDelaBella2<long double, int64_t>::Create();

