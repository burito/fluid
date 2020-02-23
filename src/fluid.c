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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <errno.h>

#include "fluid.h"
#include "log.h"


struct fluid_sim* fluid_init(float x, float y, float z, int max_depth)
{
	struct fluid_sim *sim;

	sim = malloc(sizeof(struct fluid_sim));
	if(sim == NULL)
	{
		log_fatal("malloc(fluid_sim) %s", strerror(errno));
		return NULL;
	}
	memset(sim, 0, sizeof(struct fluid_sim));
	sim->octtree = octree_init(20);
	if(sim->octtree == NULL)
	{
		log_fatal("octtree_init() failed");
		free(sim);
		return NULL;
	}
	sim->max_depth = 5; // chosen by fair dice roll
	sim->max_vortons = 10;	// why not
	sim->vortons = malloc(sim->max_vortons * sizeof(struct vorton));
	if(sim->vortons == NULL)
	{
		log_fatal("malloc(sim->vortons) %s", strerror(errno));
		free(sim);
		return NULL;
	}
	memset(sim->vortons, 0, sim->max_vortons * sizeof(struct vorton));

	// allocated all of our stuff
	return sim;
}

void fluid_end(struct fluid_sim *sim)
{
	octtree_free(sim->octtree);
	free(sim->vortons);
	free(sim);
}


void fluid_tree_update(fluid_sim *sim)
{
	int poff;
	float weight;
	float magnitude;
	vorton *parent, *child;
	vorton_list* this = sim->vortons;

	// delete the vorton list on each leaf node, and each vorton
	memset(sim->tree, 0, fluid_tree_size(sim->depth)*sizeof(vorton));
	while(this)
	{
		this->vort->next = NULL;
		this = this->next;
	}

	// walk all vortons, assigning them to a leaf node in the octtree
	this = sim->vortons;
	while(this)
	{
		parent = &sim->tree[pos_to_offset(sim, sim->depth, &this->vort->p)];
		this->vort->next = parent->next;
		parent->next = this->vort;
		this = this->next;
	}

	// weighted average the list of vortons on each leaf node
	for(int i = fluid_tree_size(sim->depth-1); i< fluid_tree_size(sim->depth); i++)
	{
		parent = &sim->tree[i];
		child = parent->next;
		weight = 0.000000001;
		while(child)
		{
			magnitude = sqrt(mag(child->w));//sqrt() may be optional
			child->magnitude = magnitude;
			weight += magnitude;
			parent->w = add(mul(child->w, magnitude), parent->w);
			parent->p = add(mul(child->p, magnitude), parent->p);
			child = child->next;
		}
		parent->w = div(parent->w, weight);
		parent->p = div(parent->p, weight);
		parent->p = sub(parent->p, sim->origin);
	}

	// weighted average each branch of the oct-tree to its parent node
	poff = fluid_tree_size(sim->depth-1);
	while(poff)
	{
		poff--;
		weight = 0.0000000001;
		parent = &sim->tree[poff];
		memset(parent, 0, sizeof(vorton));
		for(int i=0; i<8; i++)
		{
			child = &sim->tree[(poff<<3)+i];
			magnitude = sqrt(mag(child->w));//sqrt() may be optional
			child->magnitude = magnitude;
			weight += magnitude;
			parent->w = add(mul(child->w, magnitude), parent->w);
			parent->p = add(mul(child->p, magnitude), parent->p);

		}
		parent->weight = weight;
		parent->w = div(parent->w, weight);
		parent->p = div(parent->p, weight);
	}

}


int particle_inside_bound(vec3 particle, vec3 origin, vec3 volume)
{
	if(particle.x < origin.x)return 0;
	if(particle.x > origin.x+volume.x)return 0;

	if(particle.y < origin.y)return 0;
	if(particle.y > origin.y+volume.y)return 0;

	if(particle.z < origin.z)return 0;
	if(particle.z > origin.z+volume.z)return 0;
	return 1;
}



void fluid_accumulate_velocity(vec3 *v, vec3 *p, vorton *vort)
{
	float distmag;
	float oneOverDist;
	float distLaw;
	float radius = 50.0f;
	float rad2 = radius * radius;
	vec3 dist, w, result;
	dist = sub(*p, vort->p);

	distmag = mag(dist) + 0.001;
	oneOverDist = finvsqrt(distmag);
//	vect_smul(&dir, &dist, oneOverDist);
	distLaw = (distmag < rad2)
		? (oneOverDist / rad2) : (oneOverDist / distmag);

	dist = mul(dist, distLaw);
	w = mul(vort->w,
//		(1.0f / (4.0f * 3.1415926535f)) * (8.0f * rad2 * radius));
		0.636619772367f * rad2 * radius);
	result = vec3_cross(w, dist);
	*v = add(*v, result);

}



void fluid_vorton_exchange(vorton *left, vorton *right)
{
	float viscosity = 0.01f;
	float deltatime = 1.0f / 60.0f;
	vec3 delta, exchange;
	delta = sub(left->w, right->w);
	exchange = mul(delta, viscosity * deltatime);
	left->w = sub(left->w, exchange);
	right->w = add(right->w, exchange);
}


void fluid_diffuse(fluid_sim *sim)
{
	vorton* this;
	vorton* other;
	vec3 temp;
	int layer = fluid_tree_size(sim->depth - 1);
	int step;

	int cells = sim->cells - 1;
	for(int x = 0; x < cells; x++)
	for(int y = 0; y < cells; y++)
	for(int z = 0; z < cells; z++)
	{
		step = 0;
		this = sim->tree[layer+cell_offset(x,y,z)].next;
		while(this)
		{
			switch(step)
			{
			case 0:		// loop over local cell
				other = this->next; // local exchange should be mul by 2
				break;
			case 1:		// loop over x+1 cell
				other = sim->tree[layer+cell_offset(x+1,y,z)].next;
				break;
			case 2:		// loop over y+1 cell
				other = sim->tree[layer+cell_offset(x,y+1,z)].next;
				break;
			case 3:		// loop over z+1 cell
				other = sim->tree[layer+cell_offset(x,y,z+1)].next;
				break;
			}

			while(other)
			{
				fluid_vorton_exchange(this, other);
				other = other->next;
			}

			if(step < 3)
			{
				step++;
			}
			else
			{ // viscosity
				temp = mul(this->w, (1.0f/60.0f) * 0.01f);
				this->w = sub(this->w, temp);
				step = 0;
				this = this->next;
			}
		}
	}
}


void fluid_tree_velocity(fluid_sim *sim, vec3 *result, vec3 *pos)
{
	vec3 tpos;
	int layer = fluid_tree_size(sim->depth-1);
	int offset = pos_to_offset(sim, sim->depth, pos)-layer;

	tpos = sub(*pos, sim->origin);

	*result = (vec3){{0.0, 0.0, 0.0}};
	fluid_accumulate_velocity(result, &tpos, &sim->tree[offset] );
	while(offset)
	{
		int parent = offset & 0xFFFFFFF8;
		int child = offset & 7;
		for(int i=0; i<8; i++)
		if(i != child)
		{
			fluid_accumulate_velocity(result, &tpos, &sim->tree[layer+parent+i]);
		}
		offset = offset >> 3;
		layer = layer >> 3;
	}
}


void fluid_velocity_grid(fluid_sim *sim)
{
	int layer = fluid_tree_size(sim->depth - 1);

	int cells = sim->cells;
	for(int y = 1; y < cells; y++)
	for(int x = 1; x < cells; x++)
	for(int z = 1; z < cells; z++)
	{
		vec3 *v = &sim->tree[layer + cell_offset(x,y,z)].v;
		vec3 pos = {{x,y,z}};
		pos = add(mul(pos, sim->step), sim->origin);
		fluid_tree_velocity(sim, v, &pos);
	}
}

float mix(float x, float y, float a)
{
	return x * (1.0-a) + y*a;
}


void fluid_interpolate_velocity(fluid_sim *sim, vec3 *result, vec3 *pos)
{
	int layer = fluid_tree_size(sim->depth - 1);

	vec3 rpos = sub(*pos, sim->origin);
	if( rpos.x < 0.0f || rpos.x > sim->size.x )
	{
		log_warning("Flail x! x=%f, y=%f, z=%f", pos->x, pos->y, pos->z);
		return;
	}
	if( rpos.y < 0.0f || rpos.y > sim->size.y )
	{
		log_warning("Flail y! x=%f, y=%f, z=%f", pos->x, pos->y, pos->z);
		return;
	}
	if( rpos.z < 0.0f || rpos.z > sim->size.z )
	{
		log_warning("Flail z! x=%f, y=%f, z=%f", pos->x, pos->y, pos->z);
		return;
	}

	int xp, yp, zp;
	xp = (int)(rpos.x / sim->step.x);
	yp = (int)(rpos.y / sim->step.y);
	zp = (int)(rpos.z / sim->step.z);

	float xd, yd, zd;
	xd = fmod(rpos.x, sim->step.x) * sim->oneOverStep.x;
	yd = fmod(rpos.y, sim->step.y) * sim->oneOverStep.y;
	zd = fmod(rpos.z, sim->step.z) * sim->oneOverStep.z;

	float ixd, iyd, izd;
	ixd = 1.0f - xd;
	iyd = 1.0f - yd;
	izd = 1.0f - zd;

	// http://en.wikipedia.org/wiki/Trilinear_interpolation
	vec3 *a, *b, *c, *d, *e, *f, *g, *h;
	a = &sim->tree[layer+cell_offset(xp, yp, zp)].v;	// c000
	b = &sim->tree[layer+cell_offset(xp, yp, zp+1)].v;	// c001
	c = &sim->tree[layer+cell_offset(xp, yp+1, zp)].v;	// c010
	d = &sim->tree[layer+cell_offset(xp, yp+1, zp+1)].v;	// c011
	e = &sim->tree[layer+cell_offset(xp+1, yp, zp)].v;	// c100
	f = &sim->tree[layer+cell_offset(xp+1, yp, zp+1)].v;	// c101
	g = &sim->tree[layer+cell_offset(xp+1, yp+1, zp)].v;	// c110
	h = &sim->tree[layer+cell_offset(xp+1, yp+1, zp+1)].v;	// c111

	vec3 i1 = add(mul(*a, izd), mul(*b, zd));
	vec3 i2 = add(mul(*c, izd), mul(*d, zd));
	vec3 j1 = add(mul(*e, izd), mul(*f, zd));
	vec3 j2 = add(mul(*g, izd), mul(*h, zd));

	vec3 w1 = add(mul(i1, iyd), mul(i2, yd));
	vec3 w2 = add(mul(j1, iyd), mul(j2, yd));

	*result = add(mul(w1, ixd), mul(w2, xd));
}

void fluid_stretch_tilt(fluid_sim *sim)
{
	float deltatime = 1.0f / 60.0f;
	vorton_list * this;
	int layer = fluid_tree_size(sim->depth - 1);
	int offset;

	this = sim->vortons;
	while(this)
	{
		if( particle_inside_bound(this->vort->p, sim->origin, sim->size) )
		{

			vec3 p = sub(this->vort->p, sim->origin);

			int px = (int)(p.x / sim->step.x);
			int py = (int)(p.y / sim->step.y);
			int pz = (int)(p.z / sim->step.z);

			offset = layer+cell_offset(px,py,pz);
			// FIXME: this line causes a segfault
			// compute jacobian matrix
			vec3 diff;
			diff.x = sim->tree[offset].v.x - // TODO: segfault
				sim->tree[layer+cell_offset(px+1,py,pz)].v.x;
			diff.y = sim->tree[offset].v.y -
				sim->tree[layer+cell_offset(px,py+1,pz)].v.x;
			diff.z = sim->tree[offset].v.z -
				sim->tree[layer+cell_offset(px,py,pz+1)].v.x;

			mat3x3 jacobian = vec3_jacobian_vec3(diff, sim->step);
			// multiply jacobian by vorticity vector
			vec3 dw = mul(jacobian, this->vort->w);

#define TILT_FUDGE 0.2f
			// integrate with euler
			vec3 *w = &this->vort->w;
			*w = add(*w, mul(dw, deltatime * TILT_FUDGE));
#undef TILT_FUDGE
		}
		this = this->next;
	}
}

void fluid_advect_tracers(fluid_sim *sim, struct particle *particles,
		int count)
{
	float deltatime = 1.0f / 60.0f;
	vec3 velocity;

	for(int i=0; i<count; i++)
	{
		if( !particle_inside_bound(particles[i].p, sim->origin, sim->size) )
			continue;
		fluid_interpolate_velocity(sim, &velocity, &particles[i].p);
		velocity = mul(velocity, deltatime);
		particles[i].p = add(particles[i].p, velocity);
	}

}

void fluid_advect_vortons(fluid_sim *sim)
{
	float deltatime = 1.0f / 60.0f;
	vec3 velocity;
	vorton_list *this = sim->vortons;

	while(this)
	{
		if( particle_inside_bound(this->vort->p, sim->origin, sim->size) )
		{
			fluid_interpolate_velocity(sim, &velocity, &this->vort->p);
			velocity = mul(velocity, deltatime);
			this->vort->p = add(this->vort->p, velocity);
		}
		this = this->next;
	}
}


void fluid_tick(fluid_sim *sim)
{
	fluid_tree_update(sim);
	fluid_diffuse(sim);
	fluid_velocity_grid(sim);
	fluid_stretch_tilt(sim);
	fluid_advect_vortons(sim);
}

void fluid_bound(fluid_sim *sim, struct particle *pos)
{
	float step;

	step = sim->step.x * 2.0f;

	if(pos->p.x < sim->origin.x+step)
		sim->origin.x = pos->p.x - step;
	if(pos->p.x > sim->origin.x + sim->size.x - step)
		sim->size.x = pos->p.x - sim->origin.x + step;

	step = sim->step.y * 2.0f;

	if(pos->p.y < sim->origin.y+step)
		sim->origin.y = pos->p.y - step;
	if(pos->p.y > sim->origin.y + sim->size.y - step)
		sim->size.y = pos->p.y - sim->origin.y + step;

	step = sim->step.z * 2.0f;

	if(pos->p.z < sim->origin.z+step)
		sim->origin.z = pos->p.z - step;
	if(pos->p.z > sim->origin.z + sim->size.z - step)
		sim->size.z = pos->p.z - sim->origin.z + step;
}

void fluid_update_box(fluid_sim *sim)
{
	sim->step.x = sim->size.x / (float)sim->cells;
	sim->step.y = sim->size.y / (float)sim->cells;
	sim->step.z = sim->size.z / (float)sim->cells;

	sim->oneOverStep.x = 1.0f / sim->step.x;
	sim->oneOverStep.y = 1.0f / sim->step.y;
	sim->oneOverStep.z = 1.0f / sim->step.z;
}

