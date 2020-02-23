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
#include <stdlib.h
#include <errno.h>
#include <string.h>

#include "log.h"
#include "linked_list.h"


struct linked_list* linked_list_init(uint32_t size, uint32_t heads)
{
	struct linked_list *ret;
	ret = malloc(sizeof(struct linked_list));
	if(ret == NULL)
	{
		log_error("malloc(linked_list) %s", strerror(errno));
		return NULL;
	}
	ret->node_pool_size = size;
	ret->node_pool = malloc(size * sizeof(struct linked_list_node));
	if(ret->node_pool == NULL)
	{
		log_error("malloc(node_pool) %s", strerror(errno));
		free(ret);
		return NULL;
	}
	memset(ret->node_pool, 0, size * sizeof(struct linked_list_node));
	return ret;
}

void linked_list_free(struct linked_list* linked_list)
{
	free(linked_list->node_pool);
	free(linked_list);
}

int linked_list_add(struct linked_list *linked_list, int parent, int index)
{
	if(linked_list->node_count >= linked_list->node_pool_size)
	{
		log_warning("linked list is full");
		return -1;
	}
	int i = linked_list->node_count++;
	struct linked_list_node *pool = linked_list->node_pool;
	struct linked_list_node *node = &pool[i];
	node->used = 1;
	node->index = index;
	node->prev = -1;
	node->next = -1;
	if(parent<0)
	{
		return i;
	}
	struct linked_list_node *parent_node = &pool[parent];
	if(parent_node->used==0)
	{
		log_warning("requested empty parent");
		return -1;
	}
	node->prev = parent;
	node->next = parent_node->next;
	parent_node->next = i;
	if(node->next == -1)
	{
		return i;
	}
	struct linked_list_node *next = &pool[node->next];
	next->prev = i;
	return i;
}

void linked_list_remove(struct linked_list *linked_list, int index)
{
	struct linked_list_node *pool = linked_list->node_pool;
	struct linked_list_node *node = &pool[index];
	if(node->used != 1)
	{
		log_warning("Attempted to remove unused node");
		return;
	}
	if(node->prev != -1)
	{
		struct linked_list_node *prev = &pool[node->prev];
		prev->next = node->next;
	}
	if(node->next != -1)
	{
		struct linked_list_node *next = &pool[node->next];
		next->prev = node->prev;
	}
}