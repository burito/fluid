/*
Copyright (c) 2012-2016 Daniel Burke

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

   1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.

   2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

   3. This notice may not be removed or altered from any source
   distribution.
*/
#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#include <OpenGL/gl3.h>
#else
#include <GL/glew.h>
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "global.h"
#include "3dmaths.h"
#include "log.h"
#include "vr.h"
#include "fps_movement.h"
#include "shader.h"
#include "mesh_gl.h"
#include "spacemouse.h"

//#include "fluid.h"
#include "fluidtest.h"

long long time_start = 0;
float time = 0;

float step = 0.0f;

struct MESH_OPENGL *bunny;
struct GLSLSHADER *mesh_shader;

void gfx_init(void);
void gfx_end(void);
void gfx_swap(void);


int main_init(int argc, char *argv[])
{
	gfx_init();
	log_info("GL Vendor   : %s", glGetString(GL_VENDOR) );
	log_info("GL Renderer : %s", glGetString(GL_RENDERER) );
	log_info("GL Version  : %s", glGetString(GL_VERSION) );
	log_info("SL Version  : %s", glGetString(GL_SHADING_LANGUAGE_VERSION) );

	int gl_major_version = 0;
	int gl_minor_version = 0;
	glGetIntegerv(GL_MAJOR_VERSION, &gl_major_version);
	glGetIntegerv(GL_MAJOR_VERSION, &gl_minor_version);
	log_info("glGetIntVer : %d.%d", gl_major_version, gl_minor_version);
//	glEnable(GL_FRAMEBUFFER_SRGB);

	fluidtest_init();

	mesh_shader = shader_load(
		"data/shaders/mesh.vert",
		"data/shaders/mesh.frag" );
	shader_uniform(mesh_shader, "modelview");
	shader_uniform(mesh_shader, "projection");

	bunny = mesh_load("data/models/bunny/bunny.obj");
//	bunny = mesh_load("data/models/powerplant/powerplant.obj");
//	bunny = mesh_load("data/models/lpshead/heads.obj");
//	bunny = mesh_load("data/models/lpshead/head.OBJ");
//	bunny = mesh_load("data/models/buddha/buddha.obj");
//	bunny = mesh_load("data/models/hairball/hairball.obj");
//	bunny = mesh_load("data/models/sponza/sponza.obj");
//	bunny = mesh_load("data/models/San_Miguel/san-miguel.obj");
//	bunny = mesh_load("data/models/stargate/stargate.obj");
//	bunny = mesh_load("data/models/xyzrgb_dragon/xyzrgb_dragon.obj");
	if(!bunny)
	{
		return 1;
	}
//	vr_init();


//	glGenBuffers( 1, &ea_fluid_line_elements );
//	glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, ea_fluid_line_elements );
//	glBufferData( GL_ELEMENT_ARRAY_BUFFER, (sim->cells+1)*(sim->cells+1)*6 * sizeof(int), line_elements, GL_STATIC_DRAW );

	spacemouse_init();

	time_start = sys_time();
	log_info("Initialised : OK");
	return 0;   // it worked!
}


void main_end(void)
{
	mesh_free(bunny);
	spacemouse_shutdown();
	if(vr_using)
	{
		vr_end();
	}
	gfx_end();
	log_info("Shutdown    : OK");
}


// last digit of angle is x-fov, in radians
vec4 position = {{1.636183, -0.490000, 1.073543, 0.0}};
vec4 angle = {{0.021000, 0.993001, 0.000000, M_PI*0.5}};

void render(mat4x4 view, mat4x4 projection)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);
	glDepthRangef( 0.1f, 30.0f);

	mat4x4 model = mat4x4_identity();
	model = mul( model, mat4x4_rot_y(step) );		// rotate the bunny
	model = mul( model, mat4x4_translate_float(-0.5, 0, -0.5) ); // around it's own origin
	model = mul( mat4x4_translate_float( 0, 0, -2), model );	// move it 2 metres infront of the origin

	model = mul(mat4x4_translate_vec3( position.xyz ), model);	// move to player position
	model = mul(mat4x4_rot_y(-angle.y ), model);
	model = mul(mat4x4_rot_x(-angle.x ), model);

	mat4x4 modelview = mul( view, model);

	glUseProgram(mesh_shader->program);
	glUniformMatrix4fv(mesh_shader->uniforms[0], 1, GL_TRUE, modelview.f);
	glUniformMatrix4fv(mesh_shader->uniforms[1], 1, GL_TRUE, projection.f);
	glColor4f( 1., 1., 1., 1.);
	mesh_draw(bunny);

	fluidtest_draw(modelview, projection);

//	glDrawElements( GL_LINES, (sim->cells+1)*(sim->cells+1), GL_UNSIGNED_INT, 0 );
	glBindVertexArray( 0 );

	glUseProgram(0);
}


void main_loop(void)
{
	spacemouse_tick();
	fluidtest_tick();

	// do normal render loop stuff
	// if(step > 2*M_PI)
	// 	step -= 2*M_PI;
	// else
	// 	step += 0.01;

	glClearColor(0.3f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	glClearColor(0.2f, 0.2f, 0.3f, 1.0f);

	if(!vr_using)
	{
		mat4x4 projection;
		projection = mat4x4_perspective(0.1, 300, (float)vid_width / (float)vid_height*0.1, 0.1);
		float ratio = (float)vid_height / (float)vid_width;
//		projection = mat4x4_glfrustum(-0.5, 0.5, -(ratio*0.5), ratio*0.5, 1, 30);
//		projection = mat4x4_orthographic(0.1, 30, 1, (float)vid_height / (float)vid_width);
		mat4x4 view = mat4x4_translate_float(0, 0, 0); // move the camera 1m above ground
		render(view, projection);
	}
	else
	{
		vr_loop(render);
	}

	if(keys[KEY_ESCAPE])
	{
		log_info("Shutdown on : Button press (ESC)");
		killme=1;
	}

	if(keys[KEY_F9])
	{
		keys[KEY_F9] = 0;
		log_info("VR %s", (vr_using?"Shutdown":"Startup") );
		if(!vr_using)vr_init();
		else vr_end();
	}

	fps_movement(&position, &angle, 0.007);

	time = (float)(sys_time() - time_start)/(float)sys_ticksecond;
	gfx_swap();
}
