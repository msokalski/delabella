/*
DELABELLA - Delaunay triangulation library
Copyright (C) 2018 GUMIX - Marcin Sokalski
*/

#define _CRT_SECURE_NO_WARNINGS

//#define DELAUNATOR
#define Cdt

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <assert.h>

#include <vector>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include <random>

#include "delabella.h"

#include <GL/gl.h>


#ifndef _WIN32
#include <GL/glext.h>
#else
#pragma comment(lib,"OpenGL32.lib")
#pragma comment(lib,"SDL2.lib")
#undef main // on windows SDL does this weird thing
#endif

// competitors
#ifdef DELAUNATOR
#include "delaunator/delaunator-header-only.hpp"
#endif

#ifdef Cdt
#include "CDT/include/CDT.h"
#endif

PFNGLBINDBUFFERPROC glBindBuffer = 0;
PFNGLDELETEBUFFERSPROC glDeleteBuffers = 0;
PFNGLGENBUFFERSPROC glGenBuffers = 0;
PFNGLBUFFERDATAPROC glBufferData = 0;
PFNGLBUFFERSUBDATAPROC glBufferSubData = 0;
PFNGLMAPBUFFERPROC glMapBuffer = 0;
PFNGLUNMAPBUFFERPROC glUnmapBuffer = 0;
PFNGLPRIMITIVERESTARTINDEXPROC glPrimitiveRestartIndex = 0;
bool BindGL(bool prim_restart)
{
	#define BINDGL(proc) if ((*(void**)&proc = SDL_GL_GetProcAddress(#proc)) == 0) return false;
	BINDGL(glBindBuffer);
	BINDGL(glDeleteBuffers);
	BINDGL(glGenBuffers);
	BINDGL(glBufferData);
	BINDGL(glBufferSubData);
	BINDGL(glMapBuffer);
	BINDGL(glUnmapBuffer);
    if (prim_restart)
	    BINDGL(glPrimitiveRestartIndex);
	#undef BINDGL
	return true;
}


int errlog(void* stream, const char* fmt, ...)
{
	va_list arg;
	va_start(arg,fmt);
	int ret = vfprintf((FILE*)stream, fmt, arg);
	va_end(arg);
	return ret;
}

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

// wrap GL buffer, 
// so mapping works even if it doesnt
struct Buf
{
    GLuint buf;
    GLenum target;
    GLsizei size;
    void* map;
    bool mapped;

    Buf() : buf(0),target(0),size(0),map(0),mapped(false) {}

    GLuint Gen(GLenum t, GLsizei s)
    {
        target = t;
        size = s;

        GLint push = 0;
        if (target == GL_ARRAY_BUFFER)
            glGetIntegerv(GL_ARRAY_BUFFER_BINDING, &push);
        else
        if (target == GL_ELEMENT_ARRAY_BUFFER)
            glGetIntegerv(GL_ELEMENT_ARRAY_BUFFER_BINDING, &push);

        glGenBuffers(1, &buf);
        glBindBuffer(target, buf);
        glBufferData(target, size, 0, GL_STATIC_DRAW);

        glBindBuffer(target, push);
        return buf;
    }

    void Del()
    {
        Unmap();
        glDeleteBuffers(1,&buf);
    }

    void* Map()
    {
        if (map)
            return 0;

        GLint push = 0;
        if (target == GL_ARRAY_BUFFER)
            glGetIntegerv(GL_ARRAY_BUFFER_BINDING, &push);
        else
        if (target == GL_ELEMENT_ARRAY_BUFFER)
            glGetIntegerv(GL_ELEMENT_ARRAY_BUFFER_BINDING, &push);

        glBindBuffer(target, buf);
        map = glMapBuffer(target, GL_WRITE_ONLY);
        mapped = true;
        if (!map)
        {
            map = malloc(size);
            mapped = false;
        }

        glBindBuffer(target, push);
        return map;
    }

    void Unmap()
    {
        if (!map)
            return;

        GLint push = 0;
        if (target == GL_ARRAY_BUFFER)
            glGetIntegerv(GL_ARRAY_BUFFER_BINDING, &push);
        else
        if (target == GL_ELEMENT_ARRAY_BUFFER)
            glGetIntegerv(GL_ELEMENT_ARRAY_BUFFER_BINDING, &push);

        glBindBuffer(target, buf);
        if (mapped)
            glUnmapBuffer(target);
        else
        {
            glBufferSubData(target, 0, size, map);
            free(map);
        }

        map = 0;
        mapped = false;

        glBindBuffer(target, push);
    }

    void Bind()
    {
        glBindBuffer(target,buf);
    }
};

// big circle r = 2503515625, (first octant)
static const unsigned int big_circle[] =
{
    2503515625,0,
    2503514500,2373375,
    2503486825,12008400,
    2503462500,16309375,
    2503312500,31890625,
    2503051625,48198000,
    2502983032,51637449,
    2502890625,55937500,
    2502122175,83517100,
    2501975000,87815625,
    2501890625,90187500,
    2501525000,99815625,
    2501429249,102187068,
    2501349825,104113100,
    2501250000,106484375,
    2500653375,119679500,
    2499924600,134042825,
    2499820648,135967689,
    2499690625,138337500,
    2498521545,158052440,
    2497489431,173600192,
    2497187500,177890625,
    2496484375,187500000,
    2496305500,189866625,
    2495975625,194155000,
    2495015625,206125000,
    2494657820,210411135,
    2493684375,221650000,
    2493299900,225933825,
    2492187500,237890625,
    2491440000,245595625,
    2488109375,277312500,
    2487845360,279671145,
    2487629280,281586665,
    2486750591,289243812,
    2486250000,293515625,
    2485970625,295872500,
    2484375000,308984375,
    2482309375,325162500,
    2482000000,327515625,
    2481433665,331779220,
    2480138400,341326825,
    2476504375,366765000,
    2475074500,376293375,
    2474424375,380545000,
    2474062500,382890625,
    2472197460,394753345,
    2471515625,399000000,
    2469692628,410132671,
    2469376100,412034175,
    2468984375,414375000,
    2466968375,426213000,
    2466232500,430450625,
    2460549175,461831400,
    2459912240,465212055,
    2459109375,469437500,
    2456829375,481227500,
    2455999000,485447625,
    2453786625,496509500,
    2452930000,500724375,
    2450500000,512484375,
    2450013055,514807260,
    2449615935,516693580,
    2449125000,519015625,
    2447109375,528437500,
    2438437500,567109375,
    2437898775,569420800,
    2437459607,571297776,
    2434137500,585290625,
    2431015625,598125000,
    2430447500,600429375,
    2429412375,604604000,
    2427067500,613949375,
    2426484375,616250000,
    2425422076,620417793,
    2418437500,647109375,
    2416425000,654584375,
    2412109375,670312500,
    2407890625,685312500,
    2406709704,689448247,
    2404575745,696854340,
    2404038465,698705620,
    2403375000,700984375,
    2402709375,703262500,
    2399985000,712504375,
    2398757375,716626500,
    2394250625,731542500,
    2389434120,747124535,
    2388405225,750407200,
    2387112500,754509375,
    2382738743,768209976,
    2378652500,780770625,
    2377307625,784856000,
    2376562500,787109375,
    2373515625,796250000,
    2372759700,798499775,
    2372144164,800326527,
    2371384375,802575000,
    2367116980,815075265,
    2362395625,828660000,
    2360968500,832717375,
    2359822017,835960844,
    2354326000,851316375,
    2348981800,865953225,
    2347490625,869987500,
    2344125000,879015625,
    2343290625,881237500,
    2342454144,883458583,
    2341773200,885261975,
    2337500000,896484375,
    2331875000,911015625,
    2330306433,915020444,
    2325890625,926187500,
    2323015625,933375000,
    2310937500,962890625,
    2310023625,965081000,
    2306292300,973964225,
    2304615625,977925000,
    2303687500,980109375,
    2298484375,992250000,
    2291971500,1007202625,
    2291015625,1009375000,
    2285361353,1022112504,
    2274750000,1045515625,
    2270709375,1054262500,
    2269708896,1056414697,
    2268894800,1058162025,
    2267890625,1060312500,
    2262778625,1071178500,
    2260935000,1075064375,
    2256067775,1085241300,
    2255231487,1086978116,
    2254200000,1089115625,
    2248950000,1099915625,
    2247057025,1103777700,
    2245538180,1106864385,
    2232814416,1132311913,
    2231256375,1135379000,
    2229302500,1139210625,
    2223812500,1149890625,
    2222721385,1151998320,
    2221833705,1153709440,
    2216612500,1163709375,
    2214609975,1167515800,
    2208984375,1178125000,
    2207866500,1180218625,
    2205835625,1184010000,
    2201262500,1192490625,
    2188906815,1215021580,
    2184214464,1223436823,
    2182109375,1227187500,
    2180945000,1229255625,
    2172890625,1243437500,
    2166300000,1254884375,
    2165109375,1256937500,
    2162946760,1260655305,
    2159055423,1267308236,
    2158078975,1268970300,
    2156875000,1271015625,
    2140509375,1298387500,
    2136484375,1305000000,
    2130200215,1315232880,
    2127937500,1318890625,
    2121972937,1328465784,
    2119687500,1332109375,
    2113273500,1342261375,
    2110964375,1345890000,
    2109687500,1347890625,
    2104484375,1356000000,
    2102151660,1359613505,
    2093648700,1372670975,
    2082755000,1389144375,
    2080365375,1392720500,
    2072330400,1404648425,
    2064890625,1415562500,
    2063547720,1417519415,
    2062455640,1419107895,
    2061109375,1421062500,
    2055625000,1428984375,
    2054269375,1430932500,
    2051808000,1434459625,
    2035875000,1456984375,
    2033368895,1460479860,
    2031360100,1463272575,
    2021784375,1476475000,
    2012555625,1489030000,
    2009994500,1492485375,
    2004235625,1500210000,
    2002812500,1502109375,
    2001387575,1504007400,
    2000228919,1505547992,
    1995584420,1511698815,
    1992984375,1515125000,
    1983515625,1527500000,
    1973522500,1540389375,
    1968750000,1546484375,
    1951587480,1568086985,
    1948890625,1571437500,
    1947400000,1573284375,
    1941330625,1580767500,
    1938612000,1584100375,
    1937109375,1585937500,
    1928715000,1596135625,
    1919513985,1607189020,
    1918275905,1608666540,
    1916750000,1610484375,
    1903985772,1625554879,
    1891015625,1640625000,
    1884687500,1647890625,
    1883124425,1649676600,
    1881853641,1651126088,
    1880287500,1652909375,
    1872337500,1661909375,
    1869479575,1665123600,
    1861957500,1673530625,
    1859079625,1676727000,
    1851015625,1685625000,
    1848116988,1688802559,
    1845794625,1691340500,
    1824100000,1714715625,
    1822473601,1716444132,
    1821151425,1717846900,
    1812890625,1726562500,
    1811253000,1728280375,
    1802109375,1737812500,
    1799121152,1740905961,
    1790750000,1749515625,
    1789090625,1751212500,
    1786079400,1754283575,
    1779314625,1761144500,
//  1761144500,1779314625,
//  1754283575,1786079400,
//  ...
};

int main(int argc, char* argv[])
{
	#ifdef _WIN32
	SetProcessDPIAware();
	#endif

	if (argc<2)
	{
        printf("usage: %s <input.xy> [output.abc]\n", argc > 0 ? argv[0] : "<executable>");
		printf("required argument (.xy file name / number of random points) missing, terminating!\n");
		return -1;
	}

	struct MyPoint
	{
		/*long*/ double x;
		/*long*/ double y;

        bool operator == (const MyPoint& p) const
        {
            return x == p.x && y == p.y;
        }

        /*
        bool operator < (const MyPoint& p)
        {
            return x * x + y * y > p.x* p.x + p.y * p.y;
        }
        */
	};
	
	std::vector<MyPoint> cloud;
    std::vector<int> constr;

	FILE* f = fopen(argv[1],"r");
    int n = atoi(argv[1]);
    if (!f)
    {
        if (n <= 2)
        {
            printf("can't open %s file, terminating!\n", argv[1]);
            return -1;
        }
        printf("generating random %d points\n", n);
        std::random_device rd{};
        std::mt19937_64 gen{ rd() };

        std::uniform_real_distribution</*long*/ double> d(-2.503515625, +2.503515625);
        //std::normal_distribution</*long*/ double> d{0.0,2.0};
        //std::gamma_distribution</*long*/ double> d(0.1,2.0);

        for (int i = 0; i < n; i++)
        {
            //MyPoint p = { (d(gen) - 5.0)*m, (d(gen) - 5.0)*m };
            MyPoint p = { d(gen), d(gen) };
            assert(isfinite(p.x) && isfinite(p.y));
            cloud.push_back(p);
        }

        for (int i = 0; i < 100; i++)
        {
            int f = rand() % n;
            int t = (f + rand() % (n-1)) % n;
            constr.push_back(f);
            constr.push_back(t);
        }

        
        //std::sort(cloud.begin(), cloud.end());

        /*
        for (int i = 0; i < n; i++)
        {
            MyPoint p = { d(gen), 0 };
            //p.y = p.x;
            cloud.push_back(p);
        }
        MyPoint p = { 0, 1};
        cloud.push_back(p);
        */
        
	}
    else
    {
        int r,n,c;
        r = fscanf(f, "%d %d", &n, &c);

        for (int i=0; i<n; i++)
        {
            double x,y;
            // allow variety of separators and extra fields till end of the line
            int n = fscanf(f,"%lf%*[,; \v\t]%lf%*[^\n]", &x, &y);

            MyPoint p = {x,y};
            cloud.push_back(p);
        }

        for (int i = 0; i < c; i++)
        {
            int a, b;
            r = fscanf(f, "%d %d", &a, &b);
            constr.push_back(a);
            constr.push_back(b);
        }

        fclose(f);
    }

    int points = (int)cloud.size();

    #ifdef Cdt

        std::vector<CDT::V2d<double>> nodups;
        for (int i = 0; i < cloud.size(); i++)
        {
            CDT::V2d<double> v;
            v.x = cloud[i].x;
            v.y = cloud[i].y;
            nodups.push_back(v);
        }

        std::vector<CDT::Edge> edges;
        for (int c = 0; c < constr.size()/2/*CONSTRAINTS*/; c++)
        {
            CDT::Edge e
            { 
                (size_t)constr[2 * c + 0],//constraint_from[c], 
                (size_t)constr[2 * c + 1] //constraint_to[c] 
            };
            edges.push_back(e);
        }

        printf("cdt removing dups... ");
        uint64_t t0 = uSec();

        CDT::DuplicatesInfo dups = CDT::RemoveDuplicatesAndRemapEdges(nodups,edges);
        /*
        for (int c = 0; c < CONSTRAINTS; c++)
        {
            CDT::Edge e = edges[c];
            CDT::Edge f{ e.v1() - account_super_triangle,e.v2() - account_super_triangle };
            edges[c] = f;
        }
        */

        uint64_t t1 = uSec();
        printf("%d ms\n", (int)((t1 - t0) / 1000));

        printf("cdt triangulation... ");
        CDT::Triangulation<double> cdt;
        cdt.insertVertices(nodups);

        uint64_t t2 = uSec();
        printf("%d ms\n", (int)((t2 - t1) / 1000));

        if (constr.size()>=2)
        {
            printf("cdt constraining... ");
            cdt.insertEdges(edges);
            uint64_t t3 = uSec();
            printf("%d ms\n", (int)((t3 - t2) / 1000));
        }

        uint64_t t4 = uSec();
        printf("cdt erasing super tri... ");
        cdt.eraseSuperTriangle();
        uint64_t t5 = uSec();
        printf("%d ms\n", (int)((t5 - t4) / 1000));
        printf("CDT TOTAL: %d\n", (int)((t5 - t0) / 1000));

    #endif


    #ifdef DELAUNATOR
    {
        std::vector<double> coords;
        for (int i=0; i<points; i++)
        {
            coords.push_back(cloud[i].x);
            coords.push_back(cloud[i].y);
        }
        uint64_t t0 = uSec();
        printf("running delaunator...\n");
        delaunator::Delaunator* d = 0;
        try
        {
            d = new delaunator::Delaunator(coords);
        }
        catch (...)
        {
            printf("delaunator threw an exception!\n");
            d = 0;
        }
        int tris_delaunator = d ? (int)d->triangles.size() / 3 : 0;
        uint64_t t1 = uSec();
        printf("elapsed %d ms\n", (int)((t1-t0)/1000));
        printf("delaunator triangles: %d\n", tris_delaunator);
        /*
        for(std::size_t i = 0; i < d.triangles.size(); i+=3) 
        {
            printf("%d %d %d\n", (int)d.triangles[i], (int)d.triangles[i+1], (int)d.triangles[i+2]);
        }
        */
    }
    #endif

	#ifdef CRUDE_XA
	xa_pool_alloc(1000);
	#endif

	IDelaBella* idb = IDelaBella::Create();
	idb->SetErrLog(errlog, stdout);
	
    printf("running delabella...\n");
    uint64_t t6 = uSec();
	int verts = idb->Triangulate(points, &cloud.data()->x, &cloud.data()->y, sizeof(MyPoint));
	int tris_delabella = verts > 0 ? verts / 3 : 0;
    int contour = idb->GetNumBoundaryVerts();
    int non_contour = idb->GetNumInternalVerts();
	int vert_num = contour + non_contour;

    /*
    {
        printf("VERTICES DUMP:\n");
        for (int i = 0; i < cloud.size(); i++)
        {
            printf("%d: %f %f\n", i, cloud[i].x, cloud[i].y);
        }

        printf("TRIANGLES DUMP:\n");
        const DelaBella_Triangle* dela = idb->GetFirstDelaunayTriangle();
        for (int i = 0; i < tris_delabella; i++)
        {
            printf("%d: %d %d %d\n",
                dela->index,
                dela->v[0]->i,
                dela->v[1]->i,
                dela->v[2]->i);
            dela = dela->next;
        }
    }
    */

    uint64_t t7 = uSec();
    printf("elapsed %d ms\n", (int)((t7-t6)/1000));
    printf("delabella triangles: %d\n", tris_delabella);
    printf("delabella contour: %d\n", contour);

    if (constr.size()>=2)
    {
        uint64_t c0 = uSec();
        const int* data = constr.data();
        //int flips = idb->Constrain(CONSTRAINTS, constraint_from, constraint_to, sizeof(int));
        int flips = idb->Constrain(constr.size()/2, data, data + 1, 2 * sizeof(int));
        uint64_t c1 = uSec();
        printf("%d flips in %d ms\n", flips, (int)((c1 - c0) / 1000));
    }

    uint64_t t8 = uSec();
    printf("Delabella TOTAL: %d\n", (int)((t8 - t6) / 1000));

    //return 0;

	// if positive, all ok 
	if (verts<=0)
	{
        printf("nothing interesting to show, exiting!\n");
		// no points given or all points are colinear
		// make emergency call ...
		idb->Destroy();
		return -2;
	}

    // COMPARE CDT with Dela
    // if (0)
    {
        #ifdef Cdt

        int tris_cdt = (int)cdt.triangles.size();
        assert(tris_delabella == tris_cdt);

        struct MyTri
        {
            MyTri(const MyTri& tri)
            {
                p[0] = tri.p[0];
                p[1] = tri.p[1];
                p[2] = tri.p[2];
            }

            bool dot0(const MyPoint* q) const
            {
                XA_REF px = q->x;
                XA_REF py = q->y;

                XA_REF adx = (XA_REF)p[0].x - px;
                XA_REF ady = (XA_REF)p[0].y - py;
                XA_REF bdx = (XA_REF)p[1].x - px;
                XA_REF bdy = (XA_REF)p[1].y - py;
                XA_REF cdx = (XA_REF)p[2].x - px;
                XA_REF cdy = (XA_REF)p[2].y - py;

                return
                    (adx * adx + ady * ady) * (cdx * bdy - bdx * cdy) +
                    (bdx * bdx + bdy * bdy) * (adx * cdy - cdx * ady) ==
                    (cdx * cdx + cdy * cdy) * (adx * bdy - bdx * ady);
            }

            MyTri& operator = (const MyTri& tri)
            {
                p[0] = tri.p[0];
                p[1] = tri.p[1];
                p[2] = tri.p[2];
                return *this;
            }

            bool operator == (const MyTri& tri) const
            {
                if (memcmp(this, &tri, sizeof(MyTri)) == 0)
                    return true;
                return false;
            }

            bool operator < (const MyTri& tri) const
            {
                for (int i = 0; i < 3; i++)
                {
                    if (p[i].x < tri.p[i].x)
                        return true;
                    if (p[i].x > tri.p[i].x)
                        return false;
                    if (p[i].y < tri.p[i].y)
                        return true;
                    if (p[i].y > tri.p[i].y)
                        return false;
                }

                return false; // equal
            }

            MyTri(const double xy3[6])
            {
                p[0].x = xy3[0]; p[0].y = xy3[1];
                p[1].x = xy3[2]; p[1].y = xy3[3];
                p[2].x = xy3[4]; p[2].y = xy3[5];

                struct C
                {
                    bool operator () (const MyPoint& p1, const MyPoint& p2) const
                    {
                        if (p1.x < p2.x)
                            return true;
                        if (p1.x > p2.x)
                            return false;
                        if (p1.y < p2.y)
                            return true;
                        return false;
                    }
                };

                C c;
                std::sort(p, p + 3, c);
            }
            MyPoint p[3];
        };

        std::vector<MyTri> cdt_set;

        int pro = 0;
        for (int i = 0; i < tris_cdt; i++)
        {
            if (i >= pro)
            {
                int p = (int)((uint64_t)100 * i / tris_cdt);
                pro = (int)((uint64_t)(p + 1) * tris_cdt / 100);
                if (pro >= tris_cdt)
                    pro = tris_cdt - 1;
                if (i == tris_cdt - 1)
                    p = 100;
                printf("\r[%2d%s] analysing cdt %s", p, p >= 100 ? "" : "%", p >= 100 ? "\n" : "");
            }

            int a = cdt.triangles[i].vertices[0];
            int b = cdt.triangles[i].vertices[1];
            int c = cdt.triangles[i].vertices[2];

            double abc[6] =
            {
                cdt.vertices[a].x,
                cdt.vertices[a].y,
                cdt.vertices[b].x,
                cdt.vertices[b].y,
                cdt.vertices[c].x,
                cdt.vertices[c].y
            };

            MyTri tri{abc};
            cdt_set.push_back(tri);
        }
        std::sort(cdt_set.begin(), cdt_set.end());

        std::vector<MyTri> dela_set;

        const DelaBella_Triangle* dela = idb->GetFirstDelaunayTriangle();
        pro = 0;
        for (int i = 0; i < tris_delabella; i++)
        {
            if (i >= pro)
            {
                int p = (int)((uint64_t)100 * i / tris_delabella);
                pro = (int)((uint64_t)(p + 1) * tris_delabella / 100);
                if (pro >= tris_delabella)
                    pro = tris_delabella - 1;
                if (i == tris_delabella - 1)
                    p = 100;
                printf("\r[%2d%s] analysing dela %s", p, p >= 100 ? "" : "%", p >= 100 ? "\n" : "");
            }

            int a = dela->v[0]->i;
            int b = dela->v[1]->i;
            int c = dela->v[2]->i;

            double abc[6] =
            {
                dela->v[0]->x,
                dela->v[0]->y,
                dela->v[1]->x,
                dela->v[1]->y,
                dela->v[2]->x,
                dela->v[2]->y
            };

            MyTri tri{ abc };
            dela_set.push_back(tri);

            dela = dela->next;
        }
        std::sort(dela_set.begin(), dela_set.end());

        printf("COMPARING...\n");
        int diffs = -1;
        if (tris_delabella == tris_cdt)
        {
            diffs = 0;
            for (int i = 0; i < tris_delabella; i++)
            {
                if (dela_set[i] == cdt_set[i])
                    // ok, exact match
                    continue;

                // try advancing one of sets
                if (dela_set[i] < cdt_set[i])
                {
                    bool found = false;
                    for (int j = i + 1; j < tris_delabella; j++)
                    {
                        if (dela_set[j].p[0].x > cdt_set[i].p[2].x)
                            break;
                        if (dela_set[j] == cdt_set[i])
                        {
                            found = true;
                            break;
                        }
                    }

                    if (found)
                        continue;
                }
                else
                {
                    bool found = false;
                    for (int j = i + 1; j < tris_delabella; j++)
                    {
                        if (dela_set[i].p[2].x < cdt_set[j].p[0].x)
                            break;
                        if (dela_set[i] == cdt_set[j])
                        {
                            found = true;
                            break;
                        }
                    }

                    if (found)
                        continue;
                }

                /*
                if (dela_set[i].dot0(cdt_set[i].p + 0) &&
                    dela_set[i].dot0(cdt_set[i].p + 1) &&
                    dela_set[i].dot0(cdt_set[i].p + 2))
                    // probably ok, (degenerated triangulation)
                    // checking it requires polygonization on both data!
                    continue;
                */

                diffs++;
            }
        }

        if (!diffs)
            printf("COMPARE OK!\n");
        else
            printf("COMPARE FAIL %d DIFFS!\n", diffs);


        /*
        // may return false negative for degenerated triangulations
        if (dela_set != cdt_set)
            printf("COMPARE FAIL!\n");
        else
            printf("COMPARE OK!\n");
        */

        #endif
    }

    if (argc>=3)
    {
        f = fopen(argv[2],"w");
        if (f)
        {
            for (int i=0; i<tris_delabella; i++)
            {
                const DelaBella_Triangle* dela = idb->GetFirstDelaunayTriangle();
                fprintf(f,"%d %d %d\n", 
                    dela->v[0]->i,
                    dela->v[1]->i,
                    dela->v[2]->i);
                dela = dela->next;
            }
            fclose(f);
        }
    }

    SDL_Init(SDL_INIT_VIDEO);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_GL_SetAttribute(SDL_GL_ACCELERATED_VISUAL, 1);
    SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 8);
    SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE, 8);

    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 2);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1);
    //SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);

	// create viewer wnd
    int width = 800, height = 600;
    #ifdef CRUDE_XA
    const char* title = "delablella-xa";
    #else
    const char* title = "delablella";
    #endif
    SDL_Window * window = SDL_CreateWindow( title, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, width, height, SDL_WINDOW_ALLOW_HIGHDPI | SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE);
    if (!window)
    {
        printf("SDL_CreateWindow failed, terminating!\n");
        idb->Destroy();
        return -1;
    }

    SDL_GLContext context = SDL_GL_CreateContext( window );
    if (!context)
    {
        printf("SDL_GL_CreateContext failed, terminating!\n");
        idb->Destroy();
        return -1;
    }

    bool prim_restart = false;
    {
        const char* ext = (const char*)glGetString(GL_EXTENSIONS);
        while (ext)
        {
            ext = strstr(ext,"GL_NV_primitive_restart");
            if (ext)
            {
                if (ext[23]==0 || ext[23]==' ')
                {
                    prim_restart = true;
                    break;
                }
                ext += 23;
            }
        }
    }

	if (!BindGL(prim_restart))
	{
		printf("Can't bind to necessary GL functions, terminating!\n");
		idb->Destroy();
		return -1;
	}

	printf("preparing graphics...\n");

	typedef GLfloat gl_t;
	GLenum gl_e = GL_FLOAT;
	#define glLoadMatrix(m) glLoadMatrixf(m)

	// useless until glVertexAttribLPointer and accompanying shader using dmat/dvec/double arithmetic is being used
	// typedef GLdouble gl_t;
	// GLenum gl_e = GL_DOUBLE;
	// #define glLoadMatrix(m) glLoadMatrixd(m)


    // create vbo and ibo
    Buf vbo, ibo_delabella;

	vbo.Gen(GL_ARRAY_BUFFER, sizeof(gl_t[3]) * points);
    gl_t* vbo_ptr = (gl_t*)vbo.Map();

    // let's give a hand to gpu by centering vertices around 0,0
    double box[4]={(double)cloud[0].x, (double)cloud[0].y, (double)cloud[0].x, (double)cloud[0].y};
	for (int i = 0; i<points; i++)
    {
        box[0] = fmin(box[0], (double)cloud[i].x);
        box[1] = fmin(box[1], (double)cloud[i].y);
        box[2] = fmax(box[2], (double)cloud[i].x);
        box[3] = fmax(box[3], (double)cloud[i].y);
    }

    double vbo_x = 0.5 * (box[0]+box[2]);
    double vbo_y = 0.5 * (box[1]+box[3]);

	printf("box: x0:%f y0:%f x1:%f y1:%f\n", box[0], box[1], box[2], box[3]);

    box[0] -= vbo_x;
    box[1] -= vbo_y;
    box[2] -= vbo_x;
    box[3] -= vbo_y;

	for (int i = 0; i<points; i++)
    {
        vbo_ptr[3*i+0] = (gl_t)(cloud[i].x - vbo_x);
        vbo_ptr[3*i+1] = (gl_t)(cloud[i].y - vbo_y);
        vbo_ptr[3*i+2] = (gl_t)(1.0); // i%5; // color
    }
    vbo.Unmap();


    // pure indices, without: center points, restarts, loop closing
    // points may be a bit too much (cuza duplicates)
    int voronoi_indices = 2 * (vert_num + tris_delabella - 1) + contour;
    int voronoi_vertices = tris_delabella + contour;

    if (prim_restart)
        voronoi_indices += vert_num; // add primitive restarts
    else
        voronoi_indices = 4 * (vert_num + tris_delabella - 1); // almost 2x bigger ehh

    int ibo_voronoi_idx = 0;
    Buf vbo_voronoi, ibo_voronoi;
	vbo_voronoi.Gen(GL_ARRAY_BUFFER, sizeof(gl_t[3]) * voronoi_vertices);
	ibo_voronoi.Gen(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * voronoi_indices);
    gl_t* vbo_voronoi_ptr = (gl_t*)vbo_voronoi.Map();
    GLuint* ibo_voronoi_ptr = (GLuint*)ibo_voronoi.Map();

    ibo_delabella.Gen(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint[3]) * tris_delabella + sizeof(GLuint) * contour);	
    GLuint* ibo_ptr = (GLuint*)ibo_delabella.Map();
    const DelaBella_Triangle* dela = idb->GetFirstDelaunayTriangle();
	for (int i = 0; i<tris_delabella; i++)
	{
        int v0 = dela->v[0]->i;
        int v1 = dela->v[1]->i;
        int v2 = dela->v[2]->i;

        ibo_ptr[3*i+0] = (GLuint)v0;
        ibo_ptr[3*i+1] = (GLuint)v1;
        ibo_ptr[3*i+2] = (GLuint)v2;

        // put it into vbo_voronoi at 'i'

		// almost exact in hybrid
		/*
		vbo_voronoi_ptr[3 * i + 0] = (gl_t)dela->n.x;
		vbo_voronoi_ptr[3 * i + 1] = (gl_t)dela->n.y;
		vbo_voronoi_ptr[3 * i + 2] = (gl_t)(-2.0*dela->n.z);
		*/

		// less jumping on extreme zooming
        vbo_voronoi_ptr[3*i+0] = (gl_t)(-0.5 * (double)dela->n.x / (double)dela->n.z - vbo_x);
        vbo_voronoi_ptr[3*i+1] = (gl_t)(-0.5 * (double)dela->n.y / (double)dela->n.z - vbo_y);
        vbo_voronoi_ptr[3*i+2] = (gl_t)(1.0);

		dela = dela->next;
	}

	const DelaBella_Vertex* prev = idb->GetFirstBoundaryVertex();
    const DelaBella_Vertex* vert = prev->next;
    int contour_min = points-1;
    int contour_max = 0;
    for (int i = 0; i<contour; i++)    
    {
        ibo_ptr[i + 3*tris_delabella] = (GLuint)vert->i;
        contour_min = vert->i < contour_min ? vert->i : contour_min;
        contour_max = vert->i > contour_max ? vert->i : contour_max;

        // put infinite edge normal to vbo_voronoi at tris_delabella + 'i'

		// coords are overwritten with edge normals in edge postproc!
        /*
        vbo_voronoi_ptr[3*(tris_delabella+i)+0] = (gl_t)(prev->x);
        vbo_voronoi_ptr[3*(tris_delabella+i)+1] = (gl_t)(prev->y);
        vbo_voronoi_ptr[3*(tris_delabella+i)+2] = (gl_t)(0.0);
        */

        double nx = prev->y - vert->y;
        double ny = vert->x - prev->x;
        vbo_voronoi_ptr[3 * (tris_delabella + i) + 0] = (gl_t)(nx);
        vbo_voronoi_ptr[3 * (tris_delabella + i) + 1] = (gl_t)(ny);
        vbo_voronoi_ptr[3 * (tris_delabella + i) + 2] = (gl_t)(0.0);

        // create special-fan / line_strip in ibo_voronoi around this boundary vertex
        ibo_voronoi_ptr[ibo_voronoi_idx++] = (GLuint)i + tris_delabella; // begin

        // iterate all dela faces around prev
        // add their voro-vert index == dela face index
        DelaBella_Iterator it;
        const DelaBella_Triangle* t = prev->StartIterator(&it); 

        // it starts at random face, so lookup the prev->vert edge
        while (1)
        {
            if (t->index >= 0)
            {
                if (t->v[0] == prev && t->v[1] == vert ||
                    t->v[1] == prev && t->v[2] == vert ||
                    t->v[2] == prev && t->v[0] == vert)
                    break;
            }
            t = it.Next();
        }

        // now iterate around, till we're inside the boundary
        while (t->index >= 0)
        {
            ibo_voronoi_ptr[ibo_voronoi_idx++] = t->index; // end
            if (!prim_restart)
                ibo_voronoi_ptr[ibo_voronoi_idx++] = t->index; // begin of next line segment
            t = it.Next();
        }

        ibo_voronoi_ptr[ibo_voronoi_idx++] = (GLuint)( i==0 ? contour-1 : i-1 ) + tris_delabella; // loop-wrapping!
        
        if (prim_restart)
            ibo_voronoi_ptr[ibo_voronoi_idx++] = (GLuint)~0; // primitive restart

        prev = vert;
        vert = vert->next;
    }

    int voronoi_strip_indices = ibo_voronoi_idx;

    // finally, for all internal vertices
    vert = idb->GetFirstInternalVertex();
    for (int i = 0; i<non_contour; i++)
    {
        // create regular-fan / line_loop in ibo_voronoi around this internal vertex
        DelaBella_Iterator it;
        const DelaBella_Triangle* t = vert->StartIterator(&it);
        const DelaBella_Triangle* e = t;
        do
        {
            assert(t->index>=0);
            ibo_voronoi_ptr[ibo_voronoi_idx++] = t->index; // begin
            t = it.Next();
            if (!prim_restart)
                ibo_voronoi_ptr[ibo_voronoi_idx++] = t->index; // end
        } while (t!=e);
        
        if (prim_restart)
            ibo_voronoi_ptr[ibo_voronoi_idx++] = (GLuint)~0; // primitive restart

        vert = vert->next;
    }

    int voronoi_loop_indices = ibo_voronoi_idx - voronoi_strip_indices;

    assert(ibo_voronoi_idx == voronoi_indices);

    ibo_delabella.Unmap();
    vbo_voronoi.Unmap();
    ibo_voronoi.Unmap();

    /*
    #ifdef DELAUNATOR
    Buf ibo_delaunator;
    ibo_delaunator.Gen(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint[3]) * tris_delaunator);
    if (d)
    {
        ibo_ptr = (GLuint*)ibo_delaunator.Map();
        for (int i = 0; i < tris_delaunator; i++)
        {
            ibo_ptr[3 * i + 0] = (GLuint)d->triangles[3 * i + 0];
            ibo_ptr[3 * i + 1] = (GLuint)d->triangles[3 * i + 1];
            ibo_ptr[3 * i + 2] = (GLuint)d->triangles[3 * i + 2];
        }
        ibo_delaunator.Unmap();
    }
    #endif
    */

    #ifdef Cdt
    int tris_cdt = (int)cdt.triangles.size();
    Buf ibo_delaunator;
    ibo_delaunator.Gen(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint[3]) * tris_cdt);
    {
        ibo_ptr = (GLuint*)ibo_delaunator.Map();

        if (dups.duplicates.size())
        {
            // let's make inverse mapping first
            int* invmap = 0;
            int invmap_size = dups.mapping.size() - dups.duplicates.size();
            invmap = (int*)malloc(sizeof(int) * invmap_size);
            for (int i = 0; i < dups.mapping.size(); i++)
                invmap[dups.mapping[i]] = i;

            for (int i = 0; i < tris_cdt; i++)
            {
                int a = cdt.triangles[i].vertices[0];
                int b = cdt.triangles[i].vertices[1];
                int c = cdt.triangles[i].vertices[2];
                ibo_ptr[3 * i + 0] = (GLuint)invmap[c];
                ibo_ptr[3 * i + 1] = (GLuint)invmap[b];
                ibo_ptr[3 * i + 2] = (GLuint)invmap[a];
            }

            free(invmap);
        }
        else
        {
            // 1:1 mapping (no dups)
            for (int i = 0; i < tris_cdt; i++)
            {
                int a = cdt.triangles[i].vertices[0];
                int b = cdt.triangles[i].vertices[1];
                int c = cdt.triangles[i].vertices[2];
                ibo_ptr[3 * i + 0] = (GLuint)c;
                ibo_ptr[3 * i + 1] = (GLuint)b;
                ibo_ptr[3 * i + 2] = (GLuint)a;
            }
        }
        ibo_delaunator.Unmap();
    }
    #endif

    // now, everything is copied to gl, free delabella
	#ifdef CRUDE_XA
	// close pool before destroy so we won't move 
	// everything unneccessarily to the pool
	xa_pool_free(); 
	#endif

    //idb->Destroy();

    int vpw, vph;
    SDL_GL_GetDrawableSize(window, &vpw, &vph);

    double cx = 0.5 * (box[0]+box[2]);
    double cy = 0.5 * (box[1]+box[3]);
    double scale = 2.0 * fmin((double)vpw/(box[2]-box[0]),(double)vph/(box[3]-box[1]));
    int zoom = -3+(int)round(log(scale) / log(1.01));

    int drag_x, drag_y, drag_zoom;
    double drag_cx, drag_cy;
    int drag = 0;

    if (prim_restart)
    {
        glPrimitiveRestartIndex((GLuint)~0);
        glEnable(GL_PRIMITIVE_RESTART);
    }

	printf("going interactive.\n");

    for( ;; )
    {
        bool x = false;
        SDL_Event event;
        while( SDL_PollEvent( &event ) )
        {
            switch( event.type )
            {
                case SDL_MOUSEWHEEL:
                {
                    if (drag==0)
                    {
                        scale = pow(1.01, zoom);
                        int vpw, vph;
                        SDL_GL_GetDrawableSize(window, &vpw, &vph);

                        SDL_GetMouseState(&drag_x,&drag_y);

                        drag_cx = cx + 2.0 * (drag_x - vpw*0.5) / scale;
                        drag_cy = cy + 2.0 * (vph*0.5 - drag_y) / scale;

                        zoom += event.wheel.y * 10;
                        scale = pow(1.01, zoom);

                        cx = (drag_cx - 2.0 * (drag_x - vpw*0.5) / scale);
                        cy = (drag_cy - 2.0 * (vph*0.5 - drag_y) / scale);
                    }

                    break;
                }

                case SDL_MOUSEMOTION:
                {
                    if (drag == 1)
                    {
                        int dx = event.motion.x - drag_x;
                        int dy = event.motion.y - drag_y;

                        double scale = pow(1.01, zoom);
                        cx = drag_cx - 2*dx/scale;
                        cy = drag_cy + 2*dy/scale;
                    }

                    if (drag == 2)
                    {
                        int dx = event.motion.x - drag_x;
                        int dy = event.motion.y - drag_y;

                        zoom = drag_zoom - dy;

                        scale = pow(1.01, zoom);
                        int vpw, vph;
                        SDL_GL_GetDrawableSize(window, &vpw, &vph);

                        cx = (drag_cx - 2.0 * (drag_x - vpw*0.5) / scale);
                        cy = (drag_cy - 2.0 * (vph*0.5 - drag_y) / scale);                        
                    }
                    break;
                }

                case SDL_MOUSEBUTTONDOWN:
                {
                    if (!drag)
                    {
                        if (event.button.button == SDL_BUTTON_LEFT)
                        {
                            drag = 1;
                            drag_x = event.button.x;
                            drag_y = event.button.y;
                            drag_cx = cx;
                            drag_cy = cy;
                            drag_zoom = zoom;
                            SDL_CaptureMouse(SDL_TRUE);
                        }
                        if (event.button.button == SDL_BUTTON_RIGHT)
                        {
                            scale = pow(1.01, zoom);
                            int vpw, vph;
                            SDL_GL_GetDrawableSize(window, &vpw, &vph);
                            double dx = cx + 2.0 * (event.button.x - vpw*0.5) / scale;
                            double dy = cy + 2.0 * (vph*0.5 - event.button.y) / scale;

                            drag = 2;
                            drag_x = event.button.x;
                            drag_y = event.button.y;
                            drag_cx = dx;
                            drag_cy = dy;
                            drag_zoom = zoom;
                            SDL_CaptureMouse(SDL_TRUE);
                        }
                    }
                    break;
                }

                case SDL_MOUSEBUTTONUP:
                {
                    if ((event.button.button == SDL_BUTTON_LEFT && drag == 1) ||
                        (event.button.button == SDL_BUTTON_RIGHT && drag == 2))
                    {
                        drag = 0;
                        SDL_CaptureMouse(SDL_FALSE);
                    }
                    break;
                }

                case SDL_WINDOWEVENT:
                {
                    switch (event.window.event) {

                        case SDL_WINDOWEVENT_CLOSE:   // exit game
                            x = true;
                            break;

                        default:
                            break;
                    }
                    break;
                }

                case SDL_KEYUP:
                    if (event.key.keysym.sym == SDLK_ESCAPE)
                    {
                        /*
                        printf("CONSTRAINTS DUMP:\n");
                        for (int i = 0; i < CONSTRAINTS; i++)
                            printf("%d - %d\n", constraint_from[i], constraint_to[i]);
                        printf("VERTICES DUMP:\n");
                        for (int i = 0; i < cloud.size(); i++)
                        {
                            printf("%f %f\n", cloud[i].x, cloud[i].y);
                        }
                        */
                    }
                    break;
            }

            if (x)
                break;
        }

        if (x)
            break;

        int vpw, vph;
        SDL_GL_GetDrawableSize(window, &vpw, &vph);
        glViewport(0,0,vpw,vph);

        double scale = pow(1.01, zoom);

        glClearColor(0,0,0,0);
        glClear(GL_COLOR_BUFFER_BIT);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(cx - vpw/scale, cx + vpw/scale, cy - vph/scale, cy + vph/scale, -1, +1);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glScalef(1,1,0); // flatten z when not shaded

        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
        glFrontFace(GL_CW);

        // hello ancient world
        vbo.Bind();
        //glInterleavedArrays(GL_V3F,0,0); // x,y, palette_index(not yet)
		glVertexPointer(3, gl_e, 0, 0);
		glEnableClientState(GL_VERTEX_ARRAY);

        ibo_delabella.Bind();

        glColor4f(0.2f,0.2f,0.2f,1.0f);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glDrawElements(GL_TRIANGLES, /*0,points-1,*/ tris_delabella*3, GL_UNSIGNED_INT, 0);

        // put verts over fill
        glColor4f(1.0f,1.0f,0.0f,1.0f);
        glPointSize(3.0f);
        glDrawArrays(GL_POINTS, 0, points);
        glPointSize(1.0f);

        //glColor4f(0.5f,0.5f,0.5f,1.0f);
        glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glDrawElements(GL_TRIANGLES, /*0,points-1,*/ tris_delabella * 3, GL_UNSIGNED_INT, 0);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        glColor4f(0.0f,0.0f,1.0f,1.0f);
        glLineWidth(3.0f);
        glDrawElements(GL_LINE_LOOP, /*contour_min, contour_max,*/ contour, GL_UNSIGNED_INT, (GLuint*)0 + tris_delabella*3);
        glLineWidth(1.0f);

        // compare with CDT
		#if 1
        #ifdef Cdt
        glEnable(GL_BLEND);
        glBlendFunc(GL_ONE, GL_ONE);
        ibo_delaunator.Bind();
        glColor4f(0.0f,1.0f,0.0f,1.0f);
        glLineWidth(1.0f);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glDrawElements(GL_TRIANGLES, /*0,points-1,*/ tris_cdt * 3, GL_UNSIGNED_INT, 0);
        glDisable(GL_BLEND);
        #endif
		#endif

        // voronoi!
        if (0)
        {
        vbo_voronoi.Bind();
        //glInterleavedArrays(GL_V3F,0,0); // x,y, palette_index(not yet)
		glVertexPointer(3, gl_e, 0, 0);
		glEnableClientState(GL_VERTEX_ARRAY);

        const static gl_t z2w[16]=
        {
            1,0,0,0,
            0,1,0,0,
            0,0,0,1,
            0,0,0,0
        };

        glLoadMatrix(z2w);

		// voro-verts in back
		glColor4f(1.0f, 1.0f, 0.0f, 1.0f);
		glPointSize(3.0f);
		glDrawArrays(GL_POINTS, 0, tris_delabella);
		glPointSize(1.0f);
        ibo_voronoi.Bind();

        glColor4f(0.0f,0.75f,0.0f,1.0f);

        if (prim_restart)
        {
            // first, draw open cells, maybe split infinite lines to separate call? (so we could color them differently)
            glDrawElements(GL_LINE_STRIP, voronoi_strip_indices, GL_UNSIGNED_INT, (GLuint*)0);

            // then closed cells
            glDrawElements(GL_LINE_LOOP, voronoi_loop_indices, GL_UNSIGNED_INT, (GLuint*)0 + voronoi_strip_indices);
        }
        else
        {
			// first, draw open cells, maybe split infinite lines to separate call? (so we could color them differently)
			glDrawElements(GL_LINES, voronoi_strip_indices, GL_UNSIGNED_INT, (GLuint*)0);

            // then closed cells
            glDrawElements(GL_LINES, voronoi_loop_indices, GL_UNSIGNED_INT, (GLuint*)0 + voronoi_strip_indices);
        }
        }

        SDL_GL_SwapWindow(window);
        SDL_Delay(15);
    }

    vbo.Del();
    ibo_delabella.Del();

    #ifdef DELAUNATOR
    ibo_delaunator.Del();
    #endif

    vbo_voronoi.Del();
    ibo_voronoi.Del();

    SDL_GL_DeleteContext( context );
    SDL_DestroyWindow( window );
    SDL_Quit();

    #ifdef XA_VAL_LEAKS
    printf("LEAKED %d allocs\n", xa_leaks(0));
    #endif

	printf("exiting!\n");

	return 0;
}

