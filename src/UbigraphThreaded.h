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


#ifndef UBIGRAPITHREADED_H
#define UBIGRAPITHREADED_H

#include "UbigraphAPI.h"
#include "threadpool.h"

extern thread_pool_t thread_pool;

/* Basic API methods */
#define stark_ubigraph_new_vertex() ubigraph_new_vertex() 
#define stark_ubigraph_new_edge(x, y) ubigraph_new_edge(x, y) 

void    stark_ubigraph_remove_vertex(vertex_id_t x);
void    stark_ubigraph_remove_edge(edge_id_t e);

/* Vertex/edge creation when user wants to use their own id's */
void    stark_ubigraph_new_vertex_w_id(vertex_id_t x);
void    stark_ubigraph_new_edge_w_id(edge_id_t e, vertex_id_t x, vertex_id_t y);

/* Delete all vertices and edges */
void        stark_ubigraph_clear();

/* Set a vertex attribute */
void    stark_ubigraph_set_vertex_attribute(vertex_id_t x,
              const char* attribute, const char* value);

/* Vertex styles */
void    stark_ubigraph_change_vertex_style(vertex_id_t x, style_id_t s);

#define stark_ubigraph_new_vertex_style(parent_style) ubigraph_new_vertex_style(parent_style) 

void    stark_ubigraph_new_vertex_style_w_id(style_id_t s, 
              style_id_t parent_style);
void    stark_ubigraph_set_vertex_style_attribute(style_id_t s,
              const char* attribute, const char* value);

/* Set an edge attribute */
void    stark_ubigraph_set_edge_attribute(edge_id_t x,
              const char* attribute, const char* value);

/* Edge styles */
void    stark_ubigraph_change_edge_style(edge_id_t x, style_id_t s);

#define stark_ubigraph_new_edge_style(parent_style) ubigraph_new_edge_style(parent_style)

void    stark_ubigraph_new_edge_style_w_id(style_id_t s,
              style_id_t parent_style);
void    stark_ubigraph_set_edge_style_attribute(style_id_t s,
              const char* attribute, const char* value);


#endif

