
#include <stdint.h>
#include "3dmaths.h"

typedef struct vorton {
	vec3 p;	// position
	vec3 w;	// vorticity
	vec3 v;	// velocity
	float weight;	// for weighted average
	float magnitude;
	struct vorton* next;
} vorton;

typedef struct vorton_list {
	vorton * vort;
	struct vorton_list * next;
	struct vorton_list * prev;
} vorton_list;

typedef struct fluid_sim {
	int depth;
	vec3 origin;
	vec3 size;
	int cells;
	vec3 step;
	vec3 oneOverStep;
	vorton* tree;
	vorton_list* vortons;
} fluid_sim;

struct particle {
	vec3 p;
	uint8_t r, g, b, a;
};

fluid_sim * fluid_init(float x, float y, float z, int depth);
void fluid_end(fluid_sim *sim);
void fluid_tick(fluid_sim *sim);
void fluid_advect_tracers(fluid_sim *sim, struct particle *particles, int count);
void fluid_bound(fluid_sim *sim, struct particle *pos);
void fluid_update_box(fluid_sim *sim);

