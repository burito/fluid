/*
 * Naugthy Assumptions I've made
 * Integers are at least 32-bits.
 * >> shifts towards the LSB, filling the MSB's with 0
 *
 *
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <errno.h>

#include "fluid.h"
#include "log.h"


int spread_bits(int x)
{
	x = (x | (x << 16)) & 0x030000FF;
	x = (x | (x <<  8)) & 0x0300F00F;
	x = (x | (x <<  4)) & 0x030C30C3;
	x = (x | (x <<  2)) & 0x09249249;
	return x;
}

int cell_offset(int x, int y, int z)
{
	return spread_bits(x)
		| spread_bits(y) << 1
		| spread_bits(z) << 2;
}

int fluid_tree_size(int depth)
{
	return 011111111111 & (0xFFFFFFFF >> (32-(depth*3 + 1)));
}

int pos_to_offset(fluid_sim *sim, int depth, vec3 *pos)
{
	int segments = pow(2,depth);
	float sx, sy, sz;
	sx = sim->size.x / (float)segments;
	sy = sim->size.y / (float)segments;
	sz = sim->size.z / (float)segments;

	return cell_offset((int)((pos->x-sim->origin.x) / sx),
			(int)((pos->y-sim->origin.y) / sy),
			(int)((pos->z-sim->origin.z) / sz) )
		+ fluid_tree_size(depth-1);
}



fluid_sim * fluid_init(float x, float y, float z, int depth)
{
	fluid_sim* sim;
	int tree_size;

	// if you don't know what you're doing, sane default
	if(depth > 10) depth = 2;
	if(depth < 2) depth = 2;

	tree_size = fluid_tree_size(depth);

	sim = malloc(sizeof(fluid_sim));
	if(sim == NULL)
	{
		log_fatal("malloc(fluid_sim)");
		return NULL;
	}
	memset(sim, 0, sizeof(fluid_sim));
	sim->depth = depth;
	sim->tree = malloc( tree_size * sizeof(vorton) );
	if(sim->tree == NULL)
	{
		log_fatal("malloc(sim->tree)");
		free(sim);
		return NULL;
	}
	memset(sim->tree, 0, tree_size * sizeof(vorton));

	sim->size.x = x;
	sim->size.y = y;
	sim->size.z = z;

	sim->cells = pow(2,depth);

	sim->step.x = x / (float)sim->cells;
	sim->step.y = y / (float)sim->cells;
	sim->step.z = z / (float)sim->cells;

	sim->oneOverStep.x = 1.0f / sim->step.x;
	sim->oneOverStep.y = 1.0f / sim->step.y;
	sim->oneOverStep.z = 1.0f / sim->step.z;

	return sim;
}


void fluid_end(fluid_sim *sim)
{
	free(sim->tree);
	free(sim);
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
//			(1.0f / (4.0f * 3.1415926535f)) * (8.0f * rad2 * radius));
			0.636619772367f * rad2 * radius);
	result = vec3_cross(w, dist);
	*v = add(*v, result);

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
	int mx, my, mz;
	int layer = fluid_tree_size(sim->depth - 1);
	int step;

	mx = my = mz = sim->cells - 1;

	for(int x = 0; x < mx; x++)
	for(int y = 0; y < my; y++)
	for(int z = 0; z < mz; z++)
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
			{								// viscosity
				temp = mul(this->w, (1.0f/60.0f) * 0.01f);
				this->w = sub(this->w, temp);
				step = 0;
				this = this->next;
			}
		}
	}
}


void fluid_tree_velocity(fluid_sim *sim, vec3 *result,
		vec3 *pos)
{
	vec3 tpos;
	int layer = fluid_tree_size(sim->depth-1);
	int offset = pos_to_offset(sim, sim->depth, pos)-layer;
	int child, parent;

	tpos = sub(*pos, sim->origin);

	result->x = result->y = result->z = 0.0f;
	fluid_accumulate_velocity(result, &tpos, &sim->tree[offset] );
	while(offset)
	{
		parent = offset & 0xFFFFFFF8;
		child = offset & 7;
		for(int i=0; i< 8; i++)
		if(i != child)
		{
			fluid_accumulate_velocity(result, &tpos, &sim->tree[layer+parent+i]);
		}
		offset = offset >> 3;
		layer = layer>>3;
	}
}


void fluid_velocity_grid(fluid_sim *sim)
{
	int mx, my, mz;
	int layer = fluid_tree_size(sim->depth - 1);
	vec3 *v;
	vec3 pos;

	mx = my = mz = sim->cells;

	for(int x = 1; x < mx; x++)
	for(int y = 1; y < my; y++)
	for(int z = 1; z < mz; z++)
	{
		v = &sim->tree[layer + cell_offset(x,y,z)].v;
		pos.x = (float)x * sim->step.x + sim->origin.x;
		pos.y = (float)y * sim->step.y + sim->origin.y;
		pos.z = (float)z * sim->step.z + sim->origin.z;
		fluid_tree_velocity(sim, v, &pos);
	}
}

void fluid_interpolate_velocity(fluid_sim *sim, vec3 *result,
		vec3 *pos)
{
	int xp, yp, zp;
	float xd, yd, zd, ixd, iyd, izd;
	float i1, i2, j1, j2;
	float w1, w2;
	int layer = fluid_tree_size(sim->depth - 1);
	vec3 *a, *b, *c, *d, *e, *f, *g, *h;
	vec3 rpos;

	rpos = sub(*pos, sim->origin);

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
	xp = (int)(rpos.x / sim->step.x);
	yp = (int)(rpos.y / sim->step.y);
	zp = (int)(rpos.z / sim->step.z);
	xd = fmod(rpos.x, sim->step.x) * sim->oneOverStep.x;
	yd = fmod(rpos.y, sim->step.y) * sim->oneOverStep.y;
	zd = fmod(rpos.z, sim->step.z) * sim->oneOverStep.z;
	ixd = 1.0f - xd;
	iyd = 1.0f - yd;
	izd = 1.0f - zd;

	// http://en.wikipedia.org/wiki/Trilinear_interpolation
	a = &sim->tree[layer+cell_offset(xp, yp, zp)].v;		// c000
	b = &sim->tree[layer+cell_offset(xp, yp, zp+1)].v;		// c001
	c = &sim->tree[layer+cell_offset(xp, yp+1, zp)].v;		// c010
	d = &sim->tree[layer+cell_offset(xp, yp+1, zp+1)].v;	// c011
	e = &sim->tree[layer+cell_offset(xp+1, yp, zp)].v;		// c100
	f = &sim->tree[layer+cell_offset(xp+1, yp, zp+1)].v;	// c101
	g = &sim->tree[layer+cell_offset(xp+1, yp+1, zp)].v;	// c110
	h = &sim->tree[layer+cell_offset(xp+1, yp+1, zp+1)].v;	// c111

	i1 = a->x * izd + b->x * zd;
	i2 = c->x * izd + d->x * zd;
	j1 = e->x * izd + f->x * zd;
	j2 = g->x * izd + h->x * zd;
	w1 = i1 * iyd + i2 * yd;
	w2 = j1 * iyd + j2 * yd;
	result->x = w1 * ixd + w2 * xd;

	i1 = a->y * izd + b->y * zd;
	i2 = c->y * izd + d->y * zd;
	j1 = e->y * izd + f->y * zd;
	j2 = g->y * izd + h->y * zd;
	w1 = i1 * iyd + i2 * yd;
	w2 = j1 * iyd + j2 * yd;
	result->y = w1 * ixd + w2 * xd;

	i1 = a->z * izd + b->z * zd;
	i2 = c->z * izd + d->z * zd;
	j1 = e->z * izd + f->z * zd;
	j2 = g->z * izd + h->z * zd;
	w1 = i1 * iyd + i2 * yd;
	w2 = j1 * iyd + j2 * yd;
	result->z = w1 * ixd + w2 * xd;
}

void fluid_stretch_tilt(fluid_sim *sim)
{
	float deltatime = 1.0f / 60.0f;
	vorton_list * this;
	vec3 p;
	int layer = fluid_tree_size(sim->depth - 1);
	int px,py,pz;
	float u,v,w;
	float x,y,z;
	int offset;

	vec3 dw, *vw;

	vec3 j1, j2, j3;

	this = sim->vortons;
	while(this)
	{
		p = sub(this->vort->p, sim->origin);

		px = (int)(p.x / sim->step.x);
		py = (int)(p.y / sim->step.y);
		pz = (int)(p.z / sim->step.z);

		offset = layer+cell_offset(px,py,pz);
		// FIXME: this line causes a segfault
		// compute jacobian matrix
		u = sim->tree[offset].v.x - // TODO: segfault
			sim->tree[layer+cell_offset(px+1,py,pz)].v.x;
		v = sim->tree[offset].v.y -
			sim->tree[layer+cell_offset(px,py+1,pz)].v.x;
		w = sim->tree[offset].v.z -
			sim->tree[layer+cell_offset(px,py,pz+1)].v.x;

		x = sim->step.x;
		y = sim->step.y;
		z = sim->step.z;
//	1x 1y 1z   a d g
//  2x 2y 2z = b e h
//  3x 3y 3z   c f i
		j1.x = u/x;
		j1.y = u/y;
		j1.z = u/z;

		j2.x = v/x;
		j2.y = v/y;
		j2.z = v/z;

		j3.x = w/x;
		j3.y = w/y;
		j3.z = w/z;

		// multiply jacobian by vorticity vector
		vw = &this->vort->w;

		dw.x = j1.x * vw->x + j1.y * vw->y + j1.z * vw->z;
		dw.y = j2.x * vw->x + j2.y * vw->y + j2.z * vw->z;
		dw.z = j3.x * vw->x + j3.y * vw->y + j3.z * vw->z;

#define TILT_FUDGE 0.2f
		// integrate with euler
		this->vort->w.x += TILT_FUDGE * dw.x * deltatime;
		this->vort->w.y += TILT_FUDGE * dw.y * deltatime;
		this->vort->w.z += TILT_FUDGE * dw.z * deltatime;
#undef TILT_FUDGE
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
		fluid_interpolate_velocity(sim, &velocity, &this->vort->p);
		velocity = mul(velocity, deltatime);
		this->vort->p = add(this->vort->p, velocity);
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

