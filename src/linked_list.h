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

struct linked_list_node {
	int32_t used;
	int32_t next;
	int32_t prev;
	uint32_t index;
};

struct linked_list {
	uint32_t node_pool_size;
	uint32_t node_count;
	struct linked_list_pool *node_pool;
};

struct linked_list* linked_list_init(uint32_t size);
void linked_list_free(struct linked_list* linked_list);
int linked_list_add(struct linked_list *linked_list, int parent, int index);
void linked_list_remove(struct linked_list *linked_list, int index);
