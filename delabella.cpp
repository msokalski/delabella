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

struct CDelaBella : IDelaBella
{
	struct Face;

	struct Vert : DelaBella_Vertex
	{
		Unsigned28 z;
		//Face* sew; now it is inherited!

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
			return u28cmp(this, &v) < 0;
		}

		static int u28cmp(const void* a, const void* b)
		{
			const Vert* va = (const Vert*)a;
			const Vert* vb = (const Vert*)b;
			if (va->z < vb->z)
				return -1;
			if (va->z > vb->z)
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
		// Norm n;

		static Face* Alloc(Face** from)
		{
			Face* f = *from;
			*from = (Face*)f->next;
			f->next = 0;
			// init all arithmetics to 0
			return f;
		}

		void Free(Face** to)
		{
			// free all arithmetics
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

		bool dot0(const Vert& p)
		{
			return 
				(Signed62)n.x * (Signed62)((Signed15)p.x - (Signed15)v[0]->x) + 
				(Signed62)n.y * (Signed62)((Signed15)p.y - (Signed15)v[0]->y) == 
				(Signed62)n.z * (Signed62)((Signed29)((Vert*)v[0])->z - (Signed29)p.z);
		}

		bool dotP(const Vert& p)
		{
			return 
				(Signed62)n.x * (Signed62)((Signed15)p.x - (Signed15)v[0]->x) + 
				(Signed62)n.y * (Signed62)((Signed15)p.y - (Signed15)v[0]->y) > 
				(Signed62)n.z * (Signed62)((Signed29)((Vert*)v[0])->z - (Signed29)p.z);
		}

		bool dotN(const Vert& p)
		{
			return
				(Signed62)n.x * (Signed62)((Signed15)p.x - (Signed15)v[0]->x) +
				(Signed62)n.y * (Signed62)((Signed15)p.y - (Signed15)v[0]->y) <
				(Signed62)n.z * (Signed62)((Signed29)((Vert*)v[0])->z - (Signed29)p.z);
		}

		bool dotNP(const Vert& p)
		{
			return 
				(Signed62)n.x * (Signed62)((Signed15)p.x - (Signed15)v[0]->x) + 
				(Signed62)n.y * (Signed62)((Signed15)p.y - (Signed15)v[0]->y) <= 
				(Signed62)n.z * (Signed62)((Signed29)((Vert*)v[0])->z - (Signed29)p.z);
		}

		Signed62 dot(const Vert& p) const // dot
		{
			Vect d = p - *(Vert*)v[0];
			return (Signed62)n.x * (Signed62)d.x + (Signed62)n.y * (Signed62)d.y + (Signed62)n.z * (Signed62)d.z;
		}

		Norm cross() const // cross of diffs
		{
			return (*(Vert*)v[1] - *(Vert*)v[0]).cross(*(Vert*)v[2] - *(Vert*)v[0]);
		}
	};

	Vert* vert_alloc;
	Face* face_alloc;
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
		int points = inp_verts;
		std::sort(vert_alloc, vert_alloc + points);

		// rmove dups
		{
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

				// copy next no-dups block
				while (r < points && !Vert::overlap(vert_alloc + r, vert_alloc + r - 1))
					vert_alloc[w++] = vert_alloc[r++];
			}

			if (points - w)
			{
				if (errlog_proc)
					errlog_proc(errlog_file, "[WRN] detected %d dups in xy array!\n", points - w);
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

			return -points;
		}

		int hull_faces = 2 * points - 4;
		*out_hull_faces = hull_faces;

		if (max_faces < hull_faces)
		{
			if (max_faces)
			{
				//free(face_alloc);
				delete [] face_alloc;
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

		Face f; // tmp
		f.v[0] = vert_alloc + 0;
		f.v[1] = vert_alloc + 1;
		f.v[2] = vert_alloc + 2;
		f.n = f.cross();

		bool colinear = f.n.z == 0;
		int i = 3;

		/////////////////////////////////////////////////////////////////////////
		// UNTIL INPUT IS COPLANAR, GROW IT IN FORM OF A 2D CONTOUR
		/*
		. |                |         after adding     . |        ________* L
		. \ Last points to / Head     next point      . \ ______/        /
		.  *____          |             ----->        .H *____          |
		.  |\_  \_____    |                           .  |\_  \_____    |
		.   \ \_      \__* - Tail points to Last      .   \ \_      \__* T
		.    \  \_      /                             .    \  \_      /
		.     \__ \_ __/                              .     \__ \_ __/
		.        \__* - Head points to Tail           .        \__/
		*/

		Vert* head = (Vert*)f.v[0];
		Vert* tail = (Vert*)f.v[1];
		Vert* last = (Vert*)f.v[2];

		head->next = tail;
		tail->next = last;
		last->next = head;

		//while (i < points && f.dot(vert_alloc[i]) == 0)
		while (i < points && f.dot0(vert_alloc[i]))
		{
			Vert* v = vert_alloc + i;

			// it is enough to test just 1 non-zero coord
			// but we want also to test stability (assert)
			// so we calc all signs...

			// why not testing sign of dot prod of 2 normals?
			// that way we'd fall into precission problems

			Norm LvH = (*v - *last).cross(*head - *last);
			bool lvh =
				(f.n.x > 0 && LvH.x > 0) || (f.n.x < 0 && LvH.x < 0) ||
				(f.n.y > 0 && LvH.y > 0) || (f.n.y < 0 && LvH.y < 0) ||
				(f.n.z > 0 && LvH.z > 0) || (f.n.z < 0 && LvH.z < 0);

			Norm TvL = (*v - *tail).cross(*last - *tail);
			bool tvl =
				(f.n.x > 0 && TvL.x > 0) || (f.n.x < 0 && TvL.x < 0) ||
				(f.n.y > 0 && TvL.y > 0) || (f.n.y < 0 && TvL.y < 0) ||
				(f.n.z > 0 && TvL.z > 0) || (f.n.z < 0 && TvL.z < 0);

			if (lvh && !tvl) // insert new f on top of e(2,0) = (last,head)
			{
				// f.v[0] = head;
				f.v[1] = last;
				f.v[2] = v;

				last->next = v;
				v->next = head;
				tail = last;
			}
			else
				if (tvl && !lvh) // insert new f on top of e(1,2) = (tail,last)
				{
					f.v[0] = last;
					//f.v[1] = tail;
					f.v[2] = v;

					tail->next = v;
					v->next = last;
					head = last;
				}
				else
				{
					// wtf? dilithium crystals are fucked.
					assert(0);
				}

			last = v;
			i++;
		}

		if (i == points)
		{
			if (colinear)
			{
				if (errlog_proc)
					errlog_proc(errlog_file, "[WRN] all input points are colinear, returning segment list!\n");
				first_boundary_vert = head;
				first_internal_vert = 0;
				last->next = 0; // break contour, make it a list
				return -points;
			}
			else
			{
				if (points > 3)
				{
					if (errlog_proc)
						errlog_proc(errlog_file, "[NFO] all input points are cocircular.\n");
				}
				else
				{
					if (errlog_proc)
						errlog_proc(errlog_file, "[NFO] trivial case of 3 points, thank you.\n");

					first_dela_face = Face::Alloc(cache);
					first_dela_face->next = 0;
					first_hull_face = Face::Alloc(cache);
					first_hull_face->next = 0;

					first_dela_face->f[0] = first_dela_face->f[1] = first_dela_face->f[2] = first_hull_face;
					first_hull_face->f[0] = first_hull_face->f[1] = first_hull_face->f[2] = first_dela_face;

					head->sew = tail->sew = last->sew = first_hull_face;

					if (f.n.z < 0)
					{
						first_dela_face->v[0] = head;
						first_dela_face->v[1] = tail;
						first_dela_face->v[2] = last;
						first_hull_face->v[0] = last;
						first_hull_face->v[1] = tail;
						first_hull_face->v[2] = head;

						// reverse silhouette
						head->next = last;
						last->next = tail;
						tail->next = head;

						first_boundary_vert = last;
						first_internal_vert = 0;
					}
					else
					{
						first_dela_face->v[0] = last;
						first_dela_face->v[1] = tail;
						first_dela_face->v[2] = head;
						first_hull_face->v[0] = head;
						first_hull_face->v[1] = tail;
						first_hull_face->v[2] = last;

						first_boundary_vert = head;
						first_internal_vert = 0;
					}

					first_dela_face->n = first_dela_face->cross();
					first_hull_face->n = first_hull_face->cross();

					return 3;
				}

				// retract last point it will be added as a cone's top later
				last = head;
				head = (Vert*)head->next;
				i--;
			}
		}

		/////////////////////////////////////////////////////////////////////////
		// CREATE CONE HULL WITH TOP AT cloud[i] AND BASE MADE OF CONTOUR LIST
		// in 2 ways :( - depending on at which side of the contour a top vertex appears

		{
			Vert* q = vert_alloc + i;

			//if (f.dot(*q) > 0)
			if (f.dotP(*q))
			{
				Vert* p = last;
				Vert* n = (Vert*)p->next;

				Face* first_side = Face::Alloc(cache);
				first_side->v[0] = p;
				first_side->v[1] = n;
				first_side->v[2] = q;
				first_side->n = first_side->cross();
				*hull = first_side;

				p = n;
				n = (Vert*)n->next;

				Face* prev_side = first_side;
				Face* prev_base = 0;
				Face* first_base = 0;

				do
				{
					Face* base = Face::Alloc(cache);
					base->v[0] = n;
					base->v[1] = p;
					base->v[2] = last;
					base->n = base->cross();

					Face* side = Face::Alloc(cache);
					side->v[0] = p;
					side->v[1] = n;
					side->v[2] = q;
					side->n = side->cross();

					side->f[2] = base;
					base->f[2] = side;

					side->f[1] = prev_side;
					prev_side->f[0] = side;

					base->f[0] = prev_base;
					if (prev_base)
						prev_base->f[1] = base;
					else
						first_base = base;

					prev_base = base;
					prev_side = side;

					p = n;
					n = (Vert*)n->next;
				} while (n != last);

				Face* last_side = Face::Alloc(cache);
				last_side->v[0] = p;
				last_side->v[1] = n;
				last_side->v[2] = q;
				last_side->n = last_side->cross();

				last_side->f[1] = prev_side;
				prev_side->f[0] = last_side;

				last_side->f[0] = first_side;
				first_side->f[1] = last_side;

				first_base->f[0] = first_side;
				first_side->f[2] = first_base;

				last_side->f[2] = prev_base;
				prev_base->f[1] = last_side;

				i++;
			}
			else
			{
				Vert* p = last;
				Vert* n = (Vert*)p->next;

				Face* first_side = Face::Alloc(cache);
				first_side->v[0] = n;
				first_side->v[1] = p;
				first_side->v[2] = q;
				first_side->n = first_side->cross();
				*hull = first_side;

				p = n;
				n = (Vert*)n->next;

				Face* prev_side = first_side;
				Face* prev_base = 0;
				Face* first_base = 0;

				do
				{
					Face* base = Face::Alloc(cache);
					base->v[0] = p;
					base->v[1] = n;
					base->v[2] = last;
					base->n = base->cross();

					Face* side = Face::Alloc(cache);
					side->v[0] = n;
					side->v[1] = p;
					side->v[2] = q;
					side->n = side->cross();

					side->f[2] = base;
					base->f[2] = side;

					side->f[0] = prev_side;
					prev_side->f[1] = side;

					base->f[1] = prev_base;
					if (prev_base)
						prev_base->f[0] = base;
					else
						first_base = base;

					prev_base = base;
					prev_side = side;

					p = n;
					n = (Vert*)n->next;
				} while (n != last);

				Face* last_side = Face::Alloc(cache);
				last_side->v[0] = n;
				last_side->v[1] = p;
				last_side->v[2] = q;
				last_side->n = last_side->cross();

				last_side->f[0] = prev_side;
				prev_side->f[1] = last_side;

				last_side->f[1] = first_side;
				first_side->f[0] = last_side;

				first_base->f[1] = first_side;
				first_side->f[2] = first_base;

				last_side->f[2] = prev_base;
				prev_base->f[0] = last_side;

				i++;
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
		if (points <= 0 || points == 3)
			return points;

		/////////////////////////////////////////////////////////////////////////
		// ACTUAL ALGORITHM
		int proc = -1;
		for (; i < points; i++)
		{
			//ValidateHull(alloc, 2 * i - 4);

			int _proc = i * 100 / points;
			if (_proc != proc)
			{
				proc = _proc;
				printf("%d%%\n", proc);
			}

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

							s->n = s->cross();
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

		int others = 0;

		i = 0;
		Face** prev_dela = &first_dela_face;
		Face** prev_hull = &first_hull_face;
		for (int j = 0; j < hull_faces; j++)
		{
			Face* f = face_alloc + j;

			// back-link all verts to some_face
			// yea, ~6x times, sorry
			((Vert*)f->v[0])->sew = f;
			((Vert*)f->v[1])->sew = f;
			((Vert*)f->v[2])->sew = f;

			if (f->n.z < 0)
			{
				f->index = i;  // store index in dela list
				*prev_dela = f;
				prev_dela = (Face**)&f->next;
				i++;
			}
			else
			{
				f->index = others; // store index in hull list
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
		DelaBella_Iterator it;
		DelaBella_Vertex* v = first_hull_face->v[0];
		const DelaBella_Triangle* t = v->StartIterator(&it);
		const DelaBella_Triangle* e = t; // end
		
		first_boundary_vert = (Vert*)v;
		out_boundary_verts = 1;

		while (1)
		{
			if (t->n.z < 0)
			{
				int pr = it.around-1; if (pr<0) pr = 2;
				int nx = it.around+1; if (nx>2) nx = 0;

				if (t->f[pr]->n.z >= 0)
				{
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
		}

		out_hull_faces = 0;
		unique_points = 0;
		out_verts = Triangulate(&out_hull_faces);
		return out_verts;
	}	

	virtual void Destroy()
	{
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
