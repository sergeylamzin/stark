/*
    Copyright (C) 2014  Sergey Lamzin, https://github.com/sergeylamzin/stark

    This file is part of the StarK genome assembler.

    StarK is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    StarK is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

*/


#include "UbigraphThreaded.h"
#include <stdlib.h>
#include <string.h>

thread_pool_t thread_pool = 0;

struct one_id {
	int32_t (*f)(int32_t);
	int32_t a;
};

void* stark_ubigraph_1arg_wrapper(void * _arg) {
	struct one_id* arg = _arg;
	arg->f(arg->a);
	//free(arg);
	return NULL;
}

struct two_ids {
	int32_t (*f)(int32_t, int32_t);
	int32_t a;
	int32_t b;
};

void* stark_ubigraph_2arg_wrapper(void * _arg) {
	struct two_ids* arg = _arg;
	arg->f(arg->a, arg->b);
	//free(arg);
	return NULL;
}

struct three_ids {
	int32_t (*f)(int32_t, int32_t, int32_t);
	int32_t a;
	int32_t b;
	int32_t c;
};


void* stark_ubigraph_3arg_wrapper(void * _arg) {
	struct three_ids* arg = _arg;
	arg->f(arg->a, arg->b, arg->c);
	//free(arg);
	return NULL;
}

struct id_two_pointers {
	int32_t (*f)(int32_t, const char*, const char*);
	int32_t a;
	char b[256];
	char c[256];
};


void * stark_ubigraph_id_two_pointers_wrapper(void * _arg) {
	struct id_two_pointers* arg = _arg;
	arg->f(arg->a, arg->b, arg->c);
	return NULL;
	// free(arg);
}

/* Basic API methods */
void    stark_ubigraph_remove_vertex(vertex_id_t x) {
	struct one_id* arg = malloc(sizeof(struct one_id));
	arg->a = x;
	arg->f = ubigraph_remove_vertex;
	thread_pool_dispatch(thread_pool, stark_ubigraph_1arg_wrapper, arg, free, 0);
}
void    stark_ubigraph_remove_edge(edge_id_t e) {
	struct one_id* arg = malloc(sizeof(struct one_id));
	arg->a = e;
	arg->f = ubigraph_remove_edge;
	thread_pool_dispatch(thread_pool, stark_ubigraph_1arg_wrapper, arg, free, 0);
}

/* Vertex/edge creation when user wants to use their own id's */
void    stark_ubigraph_new_vertex_w_id(vertex_id_t x) {
	struct one_id* arg = malloc(sizeof(struct one_id));
	arg->a = x;
	arg->f = ubigraph_new_vertex_w_id;
	thread_pool_dispatch(thread_pool, stark_ubigraph_1arg_wrapper, arg, free, 0);
	// fprintf(stderr, "Creating vertex with id %d\n", x);
}


void    stark_ubigraph_new_edge_w_id(edge_id_t e, vertex_id_t x, vertex_id_t y) {
	struct three_ids* arg = malloc(sizeof(struct three_ids));
	arg->a = e;
	arg->b = x;
	arg->c = y;
	arg->f = ubigraph_new_edge_w_id;
	thread_pool_dispatch(thread_pool, stark_ubigraph_3arg_wrapper, arg, free, 0);
	// fprintf(stderr, "Creating edge with id %d\n", e);
}

/* Delete all vertices and edges */
void        stark_ubigraph_clear() {
	if (thread_pool) {
		thread_pool_purge(thread_pool);
	} else
		thread_pool = thread_pool_create(1);
		
	ubigraph_clear();
}

/* Set a vertex attribute */
void    stark_ubigraph_set_vertex_attribute(vertex_id_t x,
              const char* attribute, const char* value) {
	struct id_two_pointers* arg = malloc(sizeof(struct id_two_pointers));
	arg->a = x;
	memcpy(arg->b, attribute, strlen(attribute)+1);
	memcpy(arg->c, value, strlen(value)+1);
	arg->f = ubigraph_set_vertex_attribute;
	thread_pool_dispatch(thread_pool, stark_ubigraph_id_two_pointers_wrapper, arg, free, 0);
}

/* Vertex styles */
void    stark_ubigraph_change_vertex_style(vertex_id_t x, style_id_t s) {
	struct two_ids* arg = malloc(sizeof(struct two_ids));
	arg->a = x;
	arg->b = s;
	arg->f = ubigraph_change_vertex_style;
	thread_pool_dispatch(thread_pool, stark_ubigraph_2arg_wrapper, arg, free, 0);
}
void    stark_ubigraph_new_vertex_style_w_id(style_id_t s, 
              style_id_t parent_style) {
	struct two_ids* arg = malloc(sizeof(struct two_ids));
	arg->a = s;
	arg->b = parent_style;
	arg->f = ubigraph_new_vertex_style_w_id;
	thread_pool_dispatch(thread_pool, stark_ubigraph_2arg_wrapper, arg, free, 0);
}
void    stark_ubigraph_set_vertex_style_attribute(style_id_t s,
              const char* attribute, const char* value) {
	struct id_two_pointers* arg = malloc(sizeof(struct id_two_pointers));
	arg->a = s;
	memcpy(arg->b, attribute, strlen(attribute)+1);
	memcpy(arg->c, value, strlen(value)+1);
	arg->f = ubigraph_set_vertex_style_attribute;
	thread_pool_dispatch(thread_pool, stark_ubigraph_id_two_pointers_wrapper, arg, free, 0);
}

/* Set an edge attribute */
void    stark_ubigraph_set_edge_attribute(edge_id_t x,
              const char* attribute, const char* value) {
	struct id_two_pointers* arg = malloc(sizeof(struct id_two_pointers));
	arg->a = x;
	memcpy(arg->b, attribute, strlen(attribute)+1);
	memcpy(arg->c, value, strlen(value)+1);
	arg->f = ubigraph_set_edge_attribute;
	thread_pool_dispatch(thread_pool, stark_ubigraph_id_two_pointers_wrapper, arg, free, 0);
}

/* Edge styles */
void    stark_ubigraph_change_edge_style(edge_id_t x, style_id_t s) {
	struct two_ids* arg = malloc(sizeof(struct two_ids));
	arg->a = x;
	arg->b = s;
	arg->f = ubigraph_change_edge_style;
	thread_pool_dispatch(thread_pool, stark_ubigraph_2arg_wrapper, arg, free, 0);
}
void    stark_ubigraph_new_edge_style_w_id(style_id_t s,
              style_id_t parent_style) {
	struct two_ids* arg = malloc(sizeof(struct two_ids));
	arg->a = s;
	arg->b = parent_style;
	arg->f = ubigraph_new_edge_style_w_id;
	thread_pool_dispatch(thread_pool, stark_ubigraph_2arg_wrapper, arg, free, 0);
}
void    stark_ubigraph_set_edge_style_attribute(style_id_t s,
              const char* attribute, const char* value) {
	struct id_two_pointers* arg = malloc(sizeof(struct id_two_pointers));
	arg->a = s;
	memcpy(arg->b, attribute, strlen(attribute)+1);
	memcpy(arg->c, value, strlen(value)+1);
	arg->f = ubigraph_set_edge_style_attribute;
	thread_pool_dispatch(thread_pool, stark_ubigraph_id_two_pointers_wrapper, arg, free, 0);
}
