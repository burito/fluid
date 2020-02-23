/*
Copyright (c) 2020 Daniel Burke

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

#include <stdint.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#include "log.h"
#include "octtree.h"
#include "3dmaths.h"

struct octtree* octree_init(uint32_t size)
{
	struct octtree *ret;
	ret = malloc(sizeof(struct octtree));
	if(ret == NULL)
	{
		log_error("malloc(octtree) %s", strerror(errno));
		return NULL;
	}
	ret->origin = (vec3){{0.0, 0.0, 0.0}};
	ret->volume = (vec3){{1.0, 1.0, 1.0}};
	ret->node_pool_size = size;
	ret->node_pool = malloc(size * sizeof(struct octtree_node));
	if(ret->node_pool == NULL)
	{
		log_error("malloc(node_pool) %s", strerror(errno));
		free(ret);
		return NULL;
	}
	memset(ret->node_pool, 0, size * sizeof(struct octtree_node));
	return ret;
}

void octtree_init(struct octtree* octtree)
{
	free(octtree->node_pool);
	free(octtree);
}
