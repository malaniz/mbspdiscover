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



/* start: TREE Structure module */
typedef struct multibsp_tree_node {
  char                          description[128];
  struct multibsp_tree_node*    sons[100];
  struct multibsp_tree_node*    parent;
  int                           index;
  int                           length;
  int                           level; // if length == 0 then is a leaf
  char                          type[128];
} multibsp_tree_node;

typedef multibsp_tree_node* multibsp_tree_node_t;


// start: build mapping structure <<constructor>>
multibsp_tree_node_t multibsp_discover_new();
// end: build mapping structure 

// start: print mapping struct  
char* multibsp_discover_print(multibsp_tree_node_t n);
// end: print mapping struct  



/* end: TREE Structure module */
