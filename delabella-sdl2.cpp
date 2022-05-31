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

#include "delabella.h"

#include <GL/gl.h>
#include <GL/glext.h>

int errlog(void* stream, const char* fmt, ...)
{
	va_list arg;
	va_start(arg,fmt);
	int ret = vfprintf((FILE*)stream, fmt, arg);
	va_end(arg);
	return ret;
}

int main(int argc, char* argv[])
{
	if (argc<2)
	{
        printf("usage: delabella[-xa]-sdl2 <input.xy> [output.abc]\n");
		printf("required argument (.xy file name) missing, terminating!\n");
		return -1;
	}

	FILE* f = fopen(argv[1],"r");
	if (!f)
	{
		printf("can't open %s file, terminating!\n", argv[1]);
		return -1;
	}

	struct MyPoint
	{
		long double x;
		long double y;
	};
	
	std::vector<MyPoint> cloud;

	while (1)
	{
		long double x,y;
		if ( fscanf(f,"%Lf %Lf", &x, &y) != 2)
			break;
		MyPoint p = {x,y};
		cloud.push_back(p);
	}
	fclose(f);

	IDelaBella* idb = IDelaBella::Create();
	idb->SetErrLog(errlog, stdout);
	
	int verts = idb->Triangulate(cloud.size(), &cloud.data()->x, &cloud.data()->y, sizeof(MyPoint));
	int tris = verts / 3;
    int contour = idb->GetNumOutputHullVerts();

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
            for (int i=0; i<tris; i++)
            {
                const DelaBella_Triangle* dela = idb->GetFirstDelaunayTriangle();
                for (int i = 0; i<tris; i++)
                {
                    fprintf(f,"%d %d %d\n", 
                        dela->v[0]->i,
                        dela->v[1]->i,
                        dela->v[2]->i);
                    dela = dela->next;
                }
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
    GLuint vbo, ibo;
	glGenBuffers( 1, &vbo );
	glGenBuffers( 1, &ibo );
    glBindBuffer( GL_ARRAY_BUFFER, vbo );
    glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, ibo );
    glBufferData( GL_ARRAY_BUFFER, sizeof(GLfloat[2]) * cloud.size(), 0, GL_STATIC_DRAW );
    glBufferData( GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint[3]) * tris + sizeof(GLuint) * contour, 0, GL_STATIC_DRAW );

    GLfloat* vbo_ptr = (GLfloat*)glMapBuffer( GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    int points = cloud.size();
    float box[4]={(float)cloud[0].x, (float)cloud[0].y, (float)cloud[0].x, (float)cloud[0].y};
	for (int i = 0; i<points; i++)
    {
        vbo_ptr[2*i+0] = (GLfloat)cloud[i].x;
        vbo_ptr[2*i+1] = (GLfloat)cloud[i].y;

        box[0] = fmin(box[0], (float)cloud[i].x);
        box[1] = fmin(box[1], (float)cloud[i].y);
        box[2] = fmax(box[2], (float)cloud[i].x);
        box[3] = fmax(box[3], (float)cloud[i].y);
    }
    glUnmapBuffer(GL_ARRAY_BUFFER);

    GLuint* ibo_ptr = (GLuint*)glMapBuffer( GL_ELEMENT_ARRAY_BUFFER, GL_WRITE_ONLY);
	const DelaBella_Triangle* dela = idb->GetFirstDelaunayTriangle();
	for (int i = 0; i<tris; i++)
	{
        ibo_ptr[3*i+0] = (GLuint)dela->v[0]->i;
        ibo_ptr[3*i+1] = (GLuint)dela->v[1]->i;
        ibo_ptr[3*i+2] = (GLuint)dela->v[2]->i;
		dela = dela->next;
	}

	const DelaBella_Vertex* vert = idb->GetFirstHullVertex();
    for (int i = 3*tris; i<3*tris+contour; i++)    
    {
        ibo_ptr[i] = (GLuint)vert->i;
        vert = vert->next;
    }
    glUnmapBuffer(GL_ELEMENT_ARRAY_BUFFER);

    // now, everything is copied to gl, free delabella
    idb->Destroy();

    // hello ancient world
    glInterleavedArrays(GL_V2F,0,0);

    double cx = 0.5 * (box[0]+box[2]);
    double cy = 0.5 * (box[1]+box[3]);
    double dd = fmax(box[2]-box[0], box[3]-box[1]);
    int zoom = (int)round(log(dd) / log(1.01));

    int drag_x, drag_y;
    bool drag = false;

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
                    zoom += event.wheel.y * 10;
                    break;
                }

                case SDL_MOUSEMOTION:
                {
                    if (drag)
                    {
                        int dx = event.motion.x - drag_x;
                        int dy = event.motion.y - drag_y;

                        double scale = pow(1.01, zoom);
                        cx -= 2*dx/scale;
                        cy += 2*dy/scale;

                        drag_x = event.motion.x;
                        drag_y = event.motion.y;
                    }
                    break;
                }

                case SDL_MOUSEBUTTONDOWN:
                {
                    if (event.button.button == SDL_BUTTON_LEFT)
                    {
                        drag = true;
                        drag_x = event.button.x;
                        drag_y = event.button.y;
                        SDL_CaptureMouse(SDL_TRUE);
                    }
                    break;
                }

                case SDL_MOUSEBUTTONUP:
                {
                    if (event.button.button == SDL_BUTTON_LEFT && drag)
                    {
                        drag = false;
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

        glColor4f(0.2f,0.2f,0.2f,1.0f);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glDrawElements(GL_TRIANGLES, tris*3, GL_UNSIGNED_INT, 0);

        glColor4f(1.0f,0.0f,0.0f,1.0f);
        glLineWidth(3.0f);
        glDrawElements(GL_LINE_LOOP, contour, GL_UNSIGNED_INT, (GLuint*)0 + tris*3);

        glColor4f(1.0f,1.0f,1.0f,1.0f);
        glLineWidth(1.0f);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glDrawElements(GL_TRIANGLES, tris*3, GL_UNSIGNED_INT, 0);

        SDL_GL_SwapWindow(window);
        SDL_Delay(15);
    }

    glDeleteBuffers(1,&vbo);
    glDeleteBuffers(1,&ibo);

    SDL_GL_DeleteContext( context );
    SDL_DestroyWindow( window );
    SDL_Quit();

	return 0;
}