/*
DELABELLA - Delaunay triangulation library
Copyright (C) 2018 GUMIX - Marcin Sokalski
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <vector>

#define GL_GLEXT_PROTOTYPES

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include <random>

#include "delabella.h"

#include <GL/gl.h>
#include <GL/glext.h>

// competitor
#include <assert.h>
#include "delaunator/delaunator-header-only.hpp"


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

int main(int argc, char* argv[])
{
	if (argc<2)
	{
        printf("usage: delabella[-xa]-sdl2 <input.xy> [output.abc]\n");
		printf("required argument (.xy file name) missing, terminating!\n");
		return -1;
	}

	struct MyPoint
	{
		/*long*/ double x;
		/*long*/ double y;
	};
	
	std::vector<MyPoint> cloud;

	FILE* f = fopen(argv[1],"r");
    int n = atoi(argv[1]);
	if (!f)
	{
        if (n<=2)
        {
		    printf("can't open %s file, terminating!\n", argv[1]);
            return -1;
        }
        printf("generating random %d points\n", n);
       	std::random_device rd{};
    	std::mt19937_64 gen{rd()};

        //std::uniform_real_distribution</*long*/ double> d(-1.L,1.L);
        std::normal_distribution</*long*/ double> d{0.0, 2.0};
        //std::gamma_distribution</*long*/ double> d(0.1L,2.0L);

        for (int i=0; i<n; i++)
        {
            MyPoint p = { d(gen), d(gen) };
            //p.y=p.x*0.3;
            cloud.push_back(p);
        }
	}
    else
    {
        while (1)
        {
            /*
            long double x,y;
            if (fscanf(f,"%Lf %Lf", &x, &y) != 2)
                break;
            */
            double x,y;
            if (fscanf(f,"%lf %lf", &x, &y) != 2)
                break;
            MyPoint p = {x,y};
            cloud.push_back(p);
        }
        fclose(f);
    }

    int points = cloud.size();

    std::vector<double> coords;
    for (int i=0; i<points; i++)
    {
        coords.push_back(cloud[i].x);
        coords.push_back(cloud[i].y);
    }

    #if 1
    uint64_t t0 = uSec();
    printf("running delaunator...\n");
    delaunator::Delaunator d(coords);
    int tris_delaunator = d.triangles.size() / 3;
    uint64_t t1 = uSec();
    printf("elapsed %d ms\n", (int)((t1-t0)/1000));
    printf("delaunator triangles: %d\n", tris_delaunator);
    /*
    for(std::size_t i = 0; i < d.triangles.size(); i+=3) 
    {
        printf("%d %d %d\n", (int)d.triangles[i], (int)d.triangles[i+1], (int)d.triangles[i+2]);
    }
    */
    #endif
    
	IDelaBella* idb = IDelaBella::Create();
	idb->SetErrLog(errlog, stdout);
	
    printf("running delabella...\n");
    uint64_t t2 = uSec();
	int verts = idb->Triangulate(points, &cloud.data()->x, &cloud.data()->y, sizeof(MyPoint));
	int tris_delabella = verts / 3;
    int contour = idb->GetNumBoundaryVerts();

    uint64_t t3 = uSec();
    printf("elapsed %d ms\n", (int)((t3-t2)/1000));
    printf("delabella triangles: %d\n", tris_delabella);
    printf("delabella contour: %d\n", contour);
    printf("delabella is %s\n", tris_delabella == 2*points - 4 - (contour-2) ? "CORRECT" : "WRONG");


	// if positive, all ok 
	if (verts<=0)
	{
		// no points given or all points are colinear
		// make emergency call ...
		idb->Destroy();
		return -2;
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
    //SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_COMPATIBILITY);

	// create viewer wnd
    int width = 800, height = 600;
    #ifdef CRUDE_XA
    const char* title = "delablella-xa";
    #else
    const char* title = "delablella";
    #endif
    SDL_Window * window = SDL_CreateWindow( title, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, width, height, SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE);
    SDL_GLContext context = SDL_GL_CreateContext( window );

	// create vbo and ibo
    GLuint vbo, ibo_delabella, ibo_delaunator;
	glGenBuffers( 1, &vbo );
	glGenBuffers( 1, &ibo_delabella );
	glGenBuffers( 1, &ibo_delaunator );

    glBindBuffer( GL_ARRAY_BUFFER, vbo );
    glBufferData( GL_ARRAY_BUFFER, sizeof(GLfloat[3]) * points, 0, GL_STATIC_DRAW );
    GLfloat* vbo_ptr = (GLfloat*)glMapBuffer( GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    float box[4]={(float)cloud[0].x, (float)cloud[0].y, (float)cloud[0].x, (float)cloud[0].y};
	for (int i = 0; i<points; i++)
    {
        vbo_ptr[3*i+0] = (GLfloat)cloud[i].x;
        vbo_ptr[3*i+1] = (GLfloat)cloud[i].y;
        vbo_ptr[3*i+2] = i%5; // color

        box[0] = fmin(box[0], (float)cloud[i].x);
        box[1] = fmin(box[1], (float)cloud[i].y);
        box[2] = fmax(box[2], (float)cloud[i].x);
        box[3] = fmax(box[3], (float)cloud[i].y);
    }
    glUnmapBuffer(GL_ARRAY_BUFFER);

    int vert_num, vert_adv;
    const char* base = (char*)idb->GetVertexArray(&vert_num,&vert_adv);

    // pure indices, without: center points, restarts, loop closing
    // points may be a bit too much (cuza duplicates)
    int voronoi_indices = 2 * (vert_num + tris_delabella - 1) + contour;
    int voronoi_vertices = tris_delabella + contour;

    // add primitive restarts
    voronoi_indices += vert_num; // num of all verts (no dups)

    int ibo_voronoi_idx = 0;
    GLuint vbo_voronoi, ibo_voronoi;
	glGenBuffers( 1, &vbo_voronoi );
	glGenBuffers( 1, &ibo_voronoi );
    glBindBuffer(GL_ARRAY_BUFFER, vbo_voronoi);
    glBufferData( GL_ARRAY_BUFFER, sizeof(GLfloat[3]) * voronoi_vertices, 0, GL_STATIC_DRAW );
    GLfloat* vbo_voronoi_ptr = (GLfloat*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_voronoi);
    glBufferData( GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * voronoi_indices, 0, GL_STATIC_DRAW );
    GLuint* ibo_voronoi_ptr = (GLuint*)glMapBuffer(GL_ELEMENT_ARRAY_BUFFER, GL_WRITE_ONLY);

    glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, ibo_delabella );
    glBufferData( GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint[3]) * tris_delabella + sizeof(GLuint) * contour, 0, GL_STATIC_DRAW );
    GLuint* ibo_ptr = (GLuint*)glMapBuffer( GL_ELEMENT_ARRAY_BUFFER, GL_WRITE_ONLY);
	const DelaBella_Triangle* dela = idb->GetFirstDelaunayTriangle();
	for (int i = 0; i<tris_delabella; i++)
	{
        int v0 = dela->v[0]->i;
        int v1 = dela->v[1]->i;
        int v2 = dela->v[2]->i;

        ibo_ptr[3*i+0] = (GLuint)v0;
        ibo_ptr[3*i+1] = (GLuint)v1;
        ibo_ptr[3*i+2] = (GLuint)v2;

        // calc voronoi cell boundary vertex
        double x1 = cloud[v0].x, y1 = cloud[v0].y;
        double x2 = cloud[v1].x, y2 = cloud[v1].y;
        double x3 = cloud[v2].x, y3 = cloud[v2].y;

        double x12 = x1-x2, x23 = x2-x3, x31 = x3-x1;
        double y12 = y1-y2, y23 = y2-y3, y31 = y3-y1;

        double cx = (x1*x1 * y23 + x2*x2 * y31 + x3*x3 * y12 - y12 * y23 * y31) /
                    (2 * (x1 * y23 + x2 * y31 + x3 * y12));

        double cy = (y1*y1 * x23 + y2*y2 * x31 + y3*y3 * x12 - x12 * x23 * x31) / 
                    (2 * (y1 * x23 + y2 * x31 + y3 * x12));

        // put it into vbo_voronoi at 'i'
        vbo_voronoi_ptr[3*i+0] = cx;
        vbo_voronoi_ptr[3*i+1] = cy;
        vbo_voronoi_ptr[3*i+2] = 1.0;

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

        double nx = cloud[prev->i].y - cloud[vert->i].y;
        double ny = cloud[vert->i].x - cloud[prev->i].x;
        double nd = 1.0/sqrt(nx*nx + ny*ny);
        nx *= nd;
        ny *= nd;


        // put infinite edge normal to vbo_voronoi at tris_delabella + 'i'
        //double mx = 0.5 * (cloud[prev->i].x + cloud[vert->i].x);
        //double my = 0.5 * (cloud[prev->i].y + cloud[vert->i].y);
        vbo_voronoi_ptr[3*(tris_delabella+i)+0] = nx; 
        vbo_voronoi_ptr[3*(tris_delabella+i)+1] = ny; 
        vbo_voronoi_ptr[3*(tris_delabella+i)+2] = 0.0;             

        // create special-fan / line_strip in ibo_voronoi around this boundary vertex
        ibo_voronoi_ptr[ibo_voronoi_idx++] = (GLuint)i + tris_delabella;

        // iterate all dela faces around prev
        // add their voro-vert index == dela face index
        DelaBella_Iterator it;
        const DelaBella_Triangle* t = prev->StartIterator(&it); 

        // it starts at random face, so lookup the prev->vert edge
        while (1)
        {
            if (t->GetSignature()<0)
            {
                if (t->v[0] == prev && t->v[1] == vert ||
                    t->v[1] == prev && t->v[2] == vert ||
                    t->v[2] == prev && t->v[0] == vert)
                    break;
            }
            t = it.Next();
        }

        // now iterate around, till we're inside the boundary
        while (t->GetSignature()<0)
        {
            ibo_voronoi_ptr[ibo_voronoi_idx++] = t->GetListIndex();
            t = it.Next();
        }

        ibo_voronoi_ptr[ibo_voronoi_idx++] = (GLuint)( i==0 ? contour-1 : i-1 ) + tris_delabella; // loop-wrapping!
        ibo_voronoi_ptr[ibo_voronoi_idx++] = (GLuint)~0; // primitive restart

        prev = vert;
        vert = vert->next;
    }

    int voronoi_strip_indices = ibo_voronoi_idx;

    // finally, for all internal vertices
    for (int i = 0; i<vert_num; i++) // we ommit duplicated points!
    {
        const DelaBella_Vertex* vert = (const DelaBella_Vertex*)(base + i*vert_adv);

        // note: all and only bondary verts have a next pointer != 0
        // (they form an endless loop)
        if (vert->next)
            continue;

        // create regular-fan / line_loop in ibo_voronoi around this internal vertex
        DelaBella_Iterator it;
        const DelaBella_Triangle* t = vert->StartIterator(&it);
        const DelaBella_Triangle* e = t;
        do
        {
            assert(t->GetSignature()<0);
            ibo_voronoi_ptr[ibo_voronoi_idx++] = t->GetListIndex();
            t = it.Next();
        } while (t!=e);
        
        ibo_voronoi_ptr[ibo_voronoi_idx++] = (GLuint)~0; // primitive restart
    }

    int voronoi_loop_indices = ibo_voronoi_idx - voronoi_strip_indices;

    // TODO: FIX CONTOUR BUILDING FOR TRIVIAL 3-POINTS CASE !!!
    assert(ibo_voronoi_idx == voronoi_indices);

    printf("contour min:%d max:%d\n",contour_min,contour_max);

    glUnmapBuffer(GL_ELEMENT_ARRAY_BUFFER); // ibo_delabella

    glBindBuffer(GL_ARRAY_BUFFER, vbo_voronoi);
    glUnmapBuffer(GL_ARRAY_BUFFER);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_voronoi);
    glUnmapBuffer(GL_ELEMENT_ARRAY_BUFFER);

    glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, ibo_delaunator );
    glBufferData( GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint[3]) * tris_delaunator, 0, GL_STATIC_DRAW );
    ibo_ptr = (GLuint*)glMapBuffer( GL_ELEMENT_ARRAY_BUFFER, GL_WRITE_ONLY);
    for (int i = 0; i<tris_delaunator; i++)    
    {
        ibo_ptr[3*i+0] = (GLuint)d.triangles[3*i+0];
        ibo_ptr[3*i+1] = (GLuint)d.triangles[3*i+1];
        ibo_ptr[3*i+2] = (GLuint)d.triangles[3*i+2];
    }
    glUnmapBuffer(GL_ELEMENT_ARRAY_BUFFER);

    // now, everything is copied to gl, free delabella
    idb->Destroy();

    int vpw, vph;
    SDL_GL_GetDrawableSize(window, &vpw, &vph);

    double cx = 0.5 * (box[0]+box[2]);
    double cy = 0.5 * (box[1]+box[3]);
    double scale = 2.0 * fmin((double)vpw/(box[2]-box[0]),(double)vph/(box[3]-box[1]));
    int zoom = -3+(int)round(log(scale) / log(1.01));

    int drag_x, drag_y, drag_zoom;
    double drag_cx, drag_cy;
    int drag = 0;


    glPrimitiveRestartIndex((GLuint)~0);
    glEnable(GL_PRIMITIVE_RESTART);

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
                    if( event.key.keysym.sym == SDLK_ESCAPE )
                        return 0;
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
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glInterleavedArrays(GL_V3F,0,0); // x,y, palette_index(not yet)

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_delabella);

        glColor4f(0.2f,0.2f,0.2f,1.0f);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glDrawRangeElements(GL_TRIANGLES, 0,points-1, tris_delabella*3, GL_UNSIGNED_INT, 0);

        // put verts over fill
        glColor4f(1.0f,1.0f,0.0f,1.0f);
        glPointSize(3.0f);
        glDrawArrays(GL_POINTS, 0, points);
        glPointSize(1.0f);

        glColor4f(0.5f,0.5f,0.5f,1.0f);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glDrawRangeElements(GL_TRIANGLES, 0,points-1, tris_delabella * 3, GL_UNSIGNED_INT, 0);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        glColor4f(1.0f,0.0f,0.0f,1.0f);
        //glLineWidth(3.0f);
        glDrawRangeElements(GL_LINE_LOOP, contour_min, contour_max, contour, GL_UNSIGNED_INT, (GLuint*)0 + tris_delabella*3);
        //glLineWidth(1.0f);

        // delaunator
        /*
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_delaunator);
        glColor4f(1.0f,1.0f,1.0f,1.0f);
        glLineWidth(1.0f);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glDrawRangeElements(GL_TRIANGLES, 0,points-1, tris_delaunator * 3, GL_UNSIGNED_INT, 0);
        */

        // voronoi!
        glBindBuffer(GL_ARRAY_BUFFER, vbo_voronoi);
        glInterleavedArrays(GL_V3F,0,0); // x,y, palette_index(not yet)

        // voro-verts in back
        glColor4f(1.0f,1.0f,0.0f,1.0f);
        glPointSize(3.0f);
        glDrawArrays(GL_POINTS, 0, tris_delabella);
        glPointSize(1.0f);        

        const static GLfloat z2w[16]=
        {
            1,0,0,0,
            0,1,0,0,
            0,0,0,1,
            0,0,0,0
        };

        glLoadMatrixf(z2w);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_voronoi);

        // first, draw open cells
        glColor4f(0.0f,0.75f,0.0f,1.0f);
        glDrawElements(GL_LINE_STRIP/*GL_POLYGON*/, voronoi_strip_indices, GL_UNSIGNED_INT, (GLuint*)0);

        // then closed cells
        glDrawElements(GL_LINE_LOOP/*GL_POLYGON*/, voronoi_loop_indices, GL_UNSIGNED_INT, (GLuint*)0 + voronoi_strip_indices);

        SDL_GL_SwapWindow(window);
        SDL_Delay(15);
    }

    glDeleteBuffers(1,&vbo);
    glDeleteBuffers(1,&ibo_delabella);
    glDeleteBuffers(1,&ibo_delaunator);

    SDL_GL_DeleteContext( context );
    SDL_DestroyWindow( window );
    SDL_Quit();

	return 0;
}