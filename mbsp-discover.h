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
