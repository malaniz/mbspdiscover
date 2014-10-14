/*

The MIT License (MIT)

Copyright (c) 2014 Marcelo Alaniz

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include <hwloc.h>
#include "mbsp-discover.h"


multibsp_tree_node_t make_mapping(hwloc_topology_t topology, hwloc_obj_t obj, int depth ){
    char string[128];
    unsigned i;
    multibsp_tree_node_t n = NULL;
    multibsp_tree_node_t ntemp = NULL;

    if (obj->arity == 1) {
      // nodes with arity 1 are unnecesary for mapping 
      n = make_mapping(topology, obj->children[0], depth + 1);

    } else { // start: needed nodes

      n = (multibsp_tree_node_t) malloc(sizeof(multibsp_tree_node));

      n->length = obj->arity;
      hwloc_obj_snprintf(string, sizeof(string), topology, obj, "#", 0);
      sprintf(n->description, "%*s%s\n", 2*depth, "", string);


      // enter in to recursion
      for (i = 0; i < n->length; i++) {
        n->sons[i] = make_mapping(topology, obj->children[i], depth + 1);
        n->sons[i]->parent = n;
      }

      // back from recursion
      if (n->length == 0) {
        n->index = obj->os_index;
        n->level = 0;
      } else {
        n->index = n->sons[0]->index;
        n->level = n->sons[0]->level+ 1;
      }
    } // end: needed nodes
    return n;
}

multibsp_tree_node_t multibsp_discover_new()
{
  multibsp_tree_node_t mapping = NULL;


  // start: HWLOC detection code 
  hwloc_topology_t  topology;

  hwloc_topology_init(&topology);
  hwloc_topology_ignore_type(topology, HWLOC_OBJ_SOCKET);
  hwloc_topology_load(topology);
  // end: hwloc detection code 


  mapping = make_mapping(topology, hwloc_get_root_obj(topology), 0);
  return mapping;
}


char* multibsp_discover_print (multibsp_tree_node_t n) {
  char* result = (char*) malloc (sizeof(char)*4096);
  int i;

  if (n->length == 0 ){
    sprintf(result, "%d", n->index);
    return result;
  } else {
    sprintf(result, "%d [ ", n->index);
    for (i=0; i< n->length; i++) {
      sprintf(result, "%s %s ", result, multibsp_discover_print(n->sons[i]));
    }
    sprintf(result, "%s ]", result);

    return result;
  }
}


#ifdef DISCOVER_DEBUG
int main () {
  char* strp;
  multibsp_tree_node_t tree = NULL;
  tree = multibsp_discover_new();
  strp = multibsp_discover_print(tree);
  printf("%s\n", strp);
  return 0;
}
#endif
