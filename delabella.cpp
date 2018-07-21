#include <assert.h>
#include <stdio.h>
#include <search.h>
#include <malloc.h>
#include <algorithm>
#include "delabella.h" // we just need LOG() macro

// assuming BITS is max(X_BITS,Y_BITS)

typedef double Unsigned28; // 2xBITS
typedef double Signed14;   // BITS
typedef double Signed15;   // BITS + 1
typedef double Signed29;   // 2xBITS + 1
typedef double Signed31;   // 2xBITS + 3
typedef double Signed45;   // 3xBITS + 3
typedef double Signed62;   // 4xBITS + 6

static Unsigned28 s14sqr(const Signed14& s)
{
	return (Unsigned28)((Signed29)s*s);
}

struct Norm
{
	Signed45 x;
	Signed45 y;
	Signed31 z;
};

struct Vect
{
	Signed15 x, y;
	Signed29 z;

	Norm cross (const Vect& v) const // cross prod
	{
		Norm n;
		n.x = (Signed45)y*v.z - (Signed45)z*v.y;
		n.y = (Signed45)z*v.x - (Signed45)x*v.z;
		n.z = (Signed29)x*v.y - (Signed29)y*v.x;
		return n;
	}
};

struct Face;

struct Vert
{
	int i; // user's array idx
	Signed14 x, y;
	Unsigned28 z;
	Vert* next;
	Face* sew;

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
		Vert* va = (Vert*)a;
		Vert* vb = (Vert*)b;
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
			return -1;
		return 0;
	}
};

struct Face
{
	Vert* v[3];
	Face* f[3];
	Norm n;
	Face* next;

	static Face* Alloc(Face** from)
	{
		Face* f = *from;
		*from = f->next;
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
			return f[1];
		if (v[1] == p)
			return f[2];
		if (v[2] == p)
			return f[0];
		return 0;
	}

	Signed62 dot (const Vert& p) const // dot
	{
		Vect d = p - *v[0];
		return (Signed62)n.x * d.x + (Signed62)n.y * d.y + (Signed62)n.z * d.z;
	}

	Norm cross () const // cross of diffs
	{
		return (*v[1] - *v[0]).cross(*v[2] - *v[0]);
	}
};

void ValidateHull(Face* arr, int n)
{
	for (int i = 0; i < n; i++)
	{
		Face* f = arr + i;

		assert(f->f[0] && f->f[1] && f->f[2]);
		assert(f->v[0] && f->v[1] && f->v[2] );
		assert(f->v[0] != f->v[1] && f->v[1] != f->v[2] && f->v[2] != f->v[0]);

		for (int k = 0; k < 3; k++)
		{
			int fv[2] = { (k + 1) % 3, (k + 2) % 3 };
			Vert* fe[2] = { f->v[fv[0]], f->v[fv[1]] };

			Face* a = f->f[k];
			assert(a->f[0] && a->f[1] && a->f[2]);
			assert(a->v[0] && a->v[1] && a->v[2]);
			assert(a->v[0] != a->v[1] && a->v[1] != a->v[2] && a->v[2] != a->v[0]);

			int av[2] = { -1,-1 };
			for (int j = 0; j < 3; j++)
			{
				if (a->v[j] == fe[1])
					av[0] = j;
				if (a->v[j] == fe[0])
					av[1] = j;
			}

			assert(av[0] >= 0 && av[1] >= 0);
			assert(av[0]+1 == av[1] || av[0] == 2 && av[1]==0);

			int j = (av[1] + 1) % 3;

			assert(a->f[j] == f);
		}
	}
}

int DelaBella(int points, const double* xy/*[points][2]*/, int* abc/*[2*points-5][3]*/, int(*errlog)(const char* fmt, ...))
{
	if (points <= 0)
	{
		if (errlog)
			errlog("[WRN] no points given, nutting to return!\n");
		return 0;
	}

	if (points == 1)
	{
		if (errlog)
			errlog("[WRN] given single point, returning single point!\n");
		abc[0] = 0;
		return -1;
	}

	Vert* cloud = (Vert*)malloc(sizeof(Vert)*points); // new Vert[points];
	if (!cloud)
	{
		if (errlog)
			errlog("[ERR] Not enough memory, shop for some more RAM. See you!\n");
		return 0;
	}

	for (int i = 0; i < points; i++)
	{
		cloud[i].i = i;
		cloud[i].x = xy[2 * i];
		cloud[i].y = xy[2 * i + 1];
		cloud[i].z = s14sqr(cloud[i].x) + s14sqr(cloud[i].y); // x*x + y*y
	}

	//qsort(cloud, points, sizeof(Vert), Vert::u28cmp);
	std::sort(cloud, cloud + points);

	// rmove dups
	{
		int dups = 0;

		int w = 0, r = 1; // skip initial no-dups block
		while (r < points && !Vert::overlap(cloud + r,cloud + w))
		{
			w++;
			r++;
		}

		w++;

		while (r < points)
		{
			r++;

			// skip dups
			while (r < points && Vert::overlap(cloud + r, cloud + r-1))
				r++;

			// copy next no-dups block
			while (r < points && !Vert::overlap(cloud + r, cloud + r - 1))
				cloud[w++] = cloud[r++];
		}

		if (points - w)
		{
			if (errlog)
				errlog("[WRN] detected %d dups in xy array!\n", points - w);
			points = w;
		}
	}

	if (points < 3)
	{
		if (points == 2)
		{
			if (errlog)
				errlog("[WRN] all input points are colinear, returning single segment!\n");
			abc[0] = cloud[0].i;
			abc[1] = cloud[1].i;
		}
		else
		{
			if (errlog)
				errlog("[WRN] all input points are identical, returning signle point!\n");
			abc[0] = cloud[0].i;
		}

		free(cloud); // delete[] cloud;
		return -points;
	}

	Face f;
	f.v[0] = cloud + 0;
	f.v[1] = cloud + 1;
	f.v[2] = cloud + 2;
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

	Vert* head = f.v[0];
	Vert* tail = f.v[1];
	Vert* last = f.v[2];

	head->next = tail;
	tail->next = last;
	last->next = head;

	while (i < points && f.dot(cloud[i]) == 0)
	{
		Vert* v = cloud + i;

		// it is enough to test just 1 non-zero coord
		// but we want also to test stability (assert) 
		// so we calc all signs...

		// why not testing sign of dot prod of 2 normals?
		// that way we'd fall into precission problems

		Norm LvH = (*v - *last).cross(*head - *last);
		bool lvh =
			f.n.x > 0 && LvH.x > 0 || f.n.x < 0 && LvH.x < 0 ||
			f.n.y > 0 && LvH.y > 0 || f.n.y < 0 && LvH.y < 0 ||
			f.n.z > 0 && LvH.z > 0 || f.n.z < 0 && LvH.z < 0;

		Norm TvL = (*v - *tail).cross(*last - *tail);
		bool tvl =
			f.n.x > 0 && TvL.x > 0 || f.n.x < 0 && TvL.x < 0 ||
			f.n.y > 0 && TvL.y > 0 || f.n.y < 0 && TvL.y < 0 ||
			f.n.z > 0 && TvL.z > 0 || f.n.z < 0 && TvL.z < 0;

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
			if (errlog)
				errlog("[WRN] all input points are colinear, returning segment list!\n");
			Vert* v = head;
			int i = 0;
			do
			{
				abc[i++] = v->i;
				v = v->next;
			} while (v != head);
			
			free(cloud); // delete[] cloud;

			return -i;
		}
		else
		{
			if (errlog)
			{
				if (points > 3)
					errlog("[NFO] all input points are cocircular.\n");
				else
					errlog("[NFO] trivial case of 3 points, thank you.\n");
			}

			Vert* p = head->next;
			Vert* v = p->next;

			int i = 0;
			if (f.n.z < 0)
			{
				do
				{
					abc[i++] = head->i;
					abc[i++] = p->i;
					abc[i++] = v->i;
					p = v;
					v = v->next;
				} while (v != head);
			}
			else
			{
				do
				{
					abc[i++] = head->i;
					abc[i++] = v->i;
					abc[i++] = p->i;
					p = v;
					v = v->next;
				} while (v != head);
			}

			free(cloud); // delete[] cloud;

			return i;
		}
	}

	// prepare cache list
	int max_faces = 2 * points - 4; // +1 max delaunay faces!
	Face* alloc = (Face*)malloc(sizeof(Face)*max_faces); // new Face[max_faces];
	if (!alloc)
	{
		if (errlog)
			errlog("[ERR] Not enough memory, shop for some more RAM. See you!\n");
		free(cloud);
		return 0;
	}

	for (int i = 1; i < max_faces; i++)
		alloc[i - 1].next = alloc + i;
	alloc[max_faces - 1].next = 0;

	Face* cache = alloc;
	Face* hull = 0; 

	/////////////////////////////////////////////////////////////////////////
	// CREATE CONE HULL WITH TOP AT cloud[i] AND BASE MADE OF CONTOUR LIST
	// in 2 ways :( - depending on at which side of the contour a top vertex appears

	Vert* q = cloud + i;

	if (f.dot(*q) > 0)
	{
		Vert* p = last;
		Vert* n = p->next;

		Face* first_side = Face::Alloc(&cache);
		first_side->v[0] = p;
		first_side->v[1] = n;
		first_side->v[2] = q;
		first_side->n = first_side->cross();
		hull = first_side;

		p = n;
		n = n->next;

		// Vert* p = head;
		// Vert* n = p->next;

		Face* prev_side = first_side;
		Face* prev_base = 0;
		Face* first_base = 0;

		do
		{
			Face* base = Face::Alloc(&cache);
			base->v[0] = n;
			base->v[1] = p;
			base->v[2] = last;
			base->n = base->cross();

			Face* side = Face::Alloc(&cache);
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
			n = n->next;
		} while (n != last);

		Face* last_side = Face::Alloc(&cache);
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
		Vert* n = p->next;

		Face* first_side = Face::Alloc(&cache);
		first_side->v[0] = n;
		first_side->v[1] = p;
		first_side->v[2] = q;
		first_side->n = first_side->cross();
		hull = first_side;

		p = n;
		n = n->next;

		Face* prev_side = first_side;
		Face* prev_base = 0;
		Face* first_base = 0;

		// Vert* p = head;
		// Vert* n = p->next;

		do
		{
			Face* base = Face::Alloc(&cache);
			base->v[0] = p;
			base->v[1] = n;
			base->v[2] = last;
			base->n = base->cross();

			Face* side = Face::Alloc(&cache);
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
			n = n->next;
		} while (n != last);

		Face* last_side = Face::Alloc(&cache);
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

	/////////////////////////////////////////////////////////////////////////
	// ACTUAL ALGORITHM

	for (; i < points; i++)
	{
		//ValidateHull(alloc, 2 * i - 4);

		Vert* q = cloud + i;
		Vert* p = cloud + i - 1;
		Face* f = hull;

		// 1. FIND FIRST VISIBLE FACE 
		//    simply iterate around last vertex using last added triange adjecency info
		while (f->dot(*q) <= 0)
		{
			f = f->Next(p);
			if (f == hull)
			{
				// if no visible face can be located at last vertex,
				// let's run through all faces (approximately last to first), 
				// yes this is emergency fallback and should not ever happen.
				f = alloc + 2 * i - 4 - 1;
				while (f->dot(*q) <= 0)
				{
					assert(f != alloc); // no face is visible? you must be kidding!
					f--;
				}
			}
		}

		// 2. DELETE VISIBLE FACES & ADD NEW ONES
		//    (we also build silhouette (vertex loop) between visible & invisible faces)

		// push first visible face onto stack (of visible faces)
		Face* stack = f;
		f->next = f; // old trick to use list pointers as 'on-stack' markers
		while (stack)
		{
			// pop, take care of last item ptr (it's not null!)
			f = stack;
			stack = f->next;
			if (stack == f)
				stack = 0;
			f->next = 0;

			// copy parts of old face that we still need after removal
			Vert* fv[3] = { f->v[0],f->v[1],f->v[2] };
			Face* ff[3] = { f->f[0],f->f[1],f->f[2] };

			// delete visible face
			f->Free(&cache);

			// check all 3 neighbors
			for (int e = 0; e < 3; e++)
			{
				Face* n = ff[e];
				if (n && !n->next) // ensure neighbor is not processed yet & isn't on stack 
				{
					if (n->dot(*q) <= 0) // if neighbor is not visible we have slihouette edge
					{
						// build face

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

		// 3. SEW SIDES OF CONE BUILT ON SLIHOUTTE SEGMENTS

		hull = alloc + 2 * i - 4 + 1; // last added face

		// last face must contain part of the silhouette
		// (edge between its v[0] and v[1])
		Vert* entry = hull->v[0]; 

		Vert* pr = entry; 
		do
		{
			// sew pr<->nx
			Vert* nx = pr->next;
			pr->sew->f[0] = nx->sew;
			nx->sew->f[1] = pr->sew;
			pr = nx;
		} while (pr != entry);
	}

	assert(2 * i - 4 == max_faces);
	//ValidateHull(alloc, max_faces);

	i = 0;
	for (int j = 0; j < max_faces; j++)
	{
		if (alloc[j].n.z < 0)
		{
			Face* h = alloc + j;
			abc[i++] = h->v[0]->i;
			abc[i++] = h->v[1]->i;
			abc[i++] = h->v[2]->i;
		}
	}

	free(cloud); // delete[] cloud;
	free(alloc); // delete[] alloc;

	return i;
}

