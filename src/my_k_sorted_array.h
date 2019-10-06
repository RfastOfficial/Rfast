//Author: Stefanos Fafalios

#ifndef _my_sorted_array_
#define _my_sorted_array_

using namespace std;

typedef struct ARRAY_NODE {
  int index;
  double value;

} a_node;

a_node* init_array(int);
void clear_array(a_node* my_ar);
a_node* refresh_array(a_node*, int);
void k_sorted_put(a_node*, int, int, double);

#endif
