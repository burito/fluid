/*
Copyright (c) 2011,2020 Daniel Burke

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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "shader.h"
#include "log.h"
#include "fluid.h"


static void fluidtest_build_lines(void);


fluid_sim * sim;
struct particle *particles;
int n_part;

GLSLSHADER *particle_shader;
GLSLSHADER *line_shader;

GLuint va_fluid;
GLuint b_fluid;
GLuint va_fluid_line_vecs;
GLuint b_fluid_line_vecs;
vec3 *line_vecs;


void fluidtest_init(void)
{
	float s;
	vorton_list * list;
	vorton *vort;

	n_part = 30*30*30;
	particles = malloc(sizeof(struct particle)*n_part + 30*30);

	float scale = 1.0 / 30.;
	int i=0;
	for(int x=0; x<30; x++)
	for(int y=0; y<30; y++)
	for(int z=0; z<30; z++)
	{
		particles[i].p.x = (float)x*scale;
		particles[i].p.y = (float)y*scale;
		particles[i].p.z = (float)z*scale;
		particles[i].r = ((float)x*(255.0/30.0));
		particles[i].g = ((float)y*(255.0/30.0));
		particles[i].b = ((float)z*(255.0/30.0));
		particles[i].a = 255;
		i++;
	}

	s = 1.0f;

	sim = fluid_init(s,s,s, 3);

	list = malloc(sizeof(vorton_list));
	vort = malloc(sizeof(vorton));

	list->prev = list->next = NULL;
	list->vort = vort;

	memset(vort, 0, sizeof(vorton));
	vort->p.x = vort->p.y = 0.5f;
	vort->p.z = 0.5f;
	vort->w.x = 0.5f;

	sim->vortons = list;

	glGenVertexArrays(1, &va_fluid);
	glBindVertexArray(va_fluid);
	glGenBuffers(1, &b_fluid);
	glBindBuffer(GL_ARRAY_BUFFER, b_fluid);
	glBufferData(GL_ARRAY_BUFFER, n_part * sizeof(struct particle), particles, GL_DYNAMIC_DRAW);
	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(3);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 16, (void*)0);
	glVertexAttribPointer(3, 4, GL_UNSIGNED_BYTE, GL_TRUE, 16, (void*)12);
	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	particle_shader = shader_load(
		"data/shaders/particle_vertex.shader",
		"data/shaders/particle_frag.shader" );
	shader_uniform(particle_shader, "modelview");
	shader_uniform(particle_shader, "projection");

	line_shader = shader_load(
		"data/shaders/line_vertex.shader",
		"data/shaders/line_frag.shader" );
	shader_uniform(line_shader, "modelview");
	shader_uniform(line_shader, "projection");

	line_vecs = malloc((sim->cells+1)*(sim->cells+1)*6 * sizeof(vec3));
	fluidtest_build_lines();
	glGenVertexArrays(1, &va_fluid_line_vecs);
	glBindVertexArray(va_fluid_line_vecs);
	glGenBuffers(1, &b_fluid_line_vecs);
	glBindBuffer(GL_ARRAY_BUFFER, b_fluid_line_vecs);
	glBufferData(GL_ARRAY_BUFFER, (sim->cells+1)*(sim->cells+1)*6 * sizeof(vec3), line_vecs, GL_DYNAMIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 12, (void*)0);
	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void fluidtest_tick(void)
{
	// advec3 the fluid
	fluid_tick(sim);
	fluid_advect_tracers(sim, particles, n_part);
	for(int i=0; i< n_part; i++)
	{
		fluid_bound(sim, &particles[i]);
	}
	fluid_update_box(sim);
}

/*
vec3 t;
void f_diag(void)
{
	log_debug("VortonP x=%f, y=%f, z=%f\n", sim->vortons->vort->p.x, sim->vortons->vort->p.y, sim->vortons->vort->p.z);
	log_debug("VortonW x=%f, y=%f, z=%f\n", sim->vortons->vort->w.x, sim->vortons->vort->w.y, sim->vortons->vort->w.z);
	log_debug("Last x=%f, y=%f, z=%f\n", particles[n_part-1].p.x, particles[n_part-1].p.y, particles[n_part-1].p.z);
	t.x = particles[n_part-1].p.x;
	t.y = particles[n_part-1].p.y;
	t.z = particles[n_part-1].p.z;
}

void f_diag2(void)
{
	int off;
//	off = pos_to_offset(sim, sim->depth, &t);

	log_debug("Off = %d\n", off);
	log_debug("CellV x=%f, y=%f, z=%f\n", 
	sim->tree[off].v.x,
	sim->tree[off].v.y,
	sim->tree[off].v.z);
}
*/

static void fluidtest_build_lines(void)
{
	for(int x=0; x<= sim->cells; x++)
	for(int y=0; y<= sim->cells; y++)
	{
		int i = y*(sim->cells+1) + x;
		line_vecs[ (i*6)+0 ] = (vec3){{
				sim->origin.x + sim->step.x * (float)x,
				sim->origin.y + sim->step.y * (float)y,
				sim->origin.z }};
		line_vecs[ (i*6)+1 ] = (vec3){{
				sim->origin.x + sim->step.x * (float)x,
				sim->origin.y + sim->step.y * (float)y,
				sim->origin.z + sim->size.z}};
//		line_elements[ (i*3)+0 ] = (int2){(i*6)+0, (i*6)+1};

		line_vecs[ (i*6)+2 ] = (vec3){{
				sim->origin.x,
				sim->origin.y + sim->step.y * (float)y,
				sim->origin.z + sim->step.z * (float)x}};
		line_vecs[ (i*6)+3 ] = (vec3){{
				sim->origin.x + sim->size.x,
				sim->origin.y + sim->step.y * (float)y,
				sim->origin.z + sim->step.z * (float)x}};
//		line_elements[ (i*3)+1 ] = (int2){(i*6)+2, (i*6)+3};

		line_vecs[ (i*6)+4 ] = (vec3){{
				sim->origin.x + sim->step.x * (float)x,
				sim->origin.y, 
				sim->origin.z + sim->step.z * (float)y}};
		line_vecs[ (i*6)+5 ] = (vec3){{
				sim->origin.x + sim->step.x * (float)x,
				sim->origin.y + sim->size.y,
				sim->origin.z + sim->step.z * (float)y}};
//		line_elements[ (i*3)+2 ] = (int2){(i*6)+4, (i*6)+5};
	}
	glBindBuffer(GL_ARRAY_BUFFER, b_fluid_line_vecs);
	glBufferData(GL_ARRAY_BUFFER, (sim->cells+1)*(sim->cells+1)*6 * sizeof(vec3), line_vecs, GL_DYNAMIC_DRAW);
}


void fluidtest_draw(mat4x4 modelview, mat4x4 projection)
{
	// draw the tracers
	glPointSize( 3.0 );
	glUseProgram(particle_shader->prog);
	glUniformMatrix4fv(particle_shader->unif[0], 1, GL_FALSE, modelview.f);
	glUniformMatrix4fv(particle_shader->unif[1], 1, GL_FALSE, projection.f);
	glBindVertexArray(va_fluid);
	glBindBuffer(GL_ARRAY_BUFFER, b_fluid);
	glBufferData(GL_ARRAY_BUFFER, n_part * sizeof(struct particle), particles, GL_DYNAMIC_DRAW);
	glDrawArrays( GL_POINTS, 0, 30*30*30);

	// draw a bounding volume
	glUseProgram(line_shader->prog);
	glUniformMatrix4fv(line_shader->unif[0], 1, GL_FALSE, modelview.f);
	glUniformMatrix4fv(line_shader->unif[1], 1, GL_FALSE, projection.f);
	glBindVertexArray( va_fluid_line_vecs );
	fluidtest_build_lines();

	glDrawArrays( GL_LINES, 0, (sim->cells+1)*(sim->cells+1)*6);

}
