//Author: Stefanos Fafalios

#include <cstdlib>
#include "my_k_sorted_array.h"

using namespace std;

void (*_place_vals_)(a_node*, int,int, double);

void clear_array(a_node* my_ar){
  delete [] my_ar;
}

a_node* refresh_array(a_node* my_arr, int size){
  // make all indexes of the array "invalid"
  a_node* first = my_arr;
  for(int i=0; i<size;i++,my_arr++){
    my_arr->index = -1;
  }

  return first;
}

void binary_place_new_values(a_node* my_arr, int imin, int imax, double value, double index, int size){
  // calculate midpoint to cut set in half
  int imid = (imax+imin)/2;
  if(imid>size-1){
    imid=size-1;
  }

  // three-way comparison
  if (my_arr[imid].value > value){
    if(imid == 0){
      for(int i = size-1; i > 0;i--){
        my_arr[i] = my_arr[i-1];
      }
      my_arr[0].index = index;
      my_arr[0].value = value;
      return;
    }

    if(my_arr[imid-1].value <= value){
      for(int i = size-1; i > imid;i--){
        my_arr[i] = my_arr[i-1];
      }
      my_arr[imid].index = index;
      my_arr[imid].value = value;

      return;
    }
    // key is in lower subset
    binary_place_new_values(my_arr, imin, imid-1, value, index, size);
  }
  else if (my_arr[imid].value < value){

    if(my_arr[imid+1].value >= value){
      for(int i = size-1; i > imid+1;i--){
        my_arr[i] = my_arr[i-1];
      }
      my_arr[imid+1].index = index;
      my_arr[imid+1].value = value;

      return;
    }
    // key is in upper subset
    binary_place_new_values(my_arr, imid+1, imax, value, index, size);
  }
  else{
    // key has been found so insert it after imax
    for(int i = size-1; i > imid+1;i--){
      my_arr[i] = my_arr[i-1];
    }
    my_arr[imid+1].index = index;
    my_arr[imid+1].value = value;
  }
  return;
}

// suitable for big size eg more than 100
void place_new_values(a_node* my_arr, int size,int index, double value){
  // case value should not be inserted in the list at all
  if(my_arr[size-1].index!=-1 && my_arr[size-1].value <= value){
    return;
  }

  if(index == 0){
    my_arr[0].index = 0;
    my_arr[0].value = value;
    return;
  }
  //case one element in array
  if(index-1 == 0){
    if(my_arr[0].value <= value){
      my_arr[1].index = 1;
      my_arr[1].value = value;
    }
    else{
      my_arr[1] = my_arr[0];
      my_arr[0].value = value;
      my_arr[0].index = 1;
    }
    return;
  }
  // case array not full
  if(index < size){
    //case empty array
    //case element should be placed after the last valid value
    if(my_arr[index-1].value <= value){
      my_arr[index].index = index;
      my_arr[index].value = value;
      return;
    }
    //case the element should be placed between 0 and index-1
    binary_place_new_values(my_arr, 0, index-1, value, index, size);
  }
  else{
    binary_place_new_values(my_arr, 0, size-1, value, index, size);
  }

  return;
}

// suitable for smaller size eg less than 100
void place_new_values2(a_node* my_arr, int size,int index, double value){
  // case value should not be inserted in the list at all
  if(my_arr[size-1].index!=-1 && my_arr[size-1].value <= value){
    return;
  }

  if(index == 0){
    my_arr[0].index = 0;
    my_arr[0].value = value;
    return;
  }
  //case one element in array
  if(index-1 == 0){
    if(my_arr[0].value <= value){
      my_arr[1].index = 1;
      my_arr[1].value = value;
    }
    else{
      my_arr[1] = my_arr[0];
      my_arr[0].value = value;
      my_arr[0].index = 1;
    }
    return;
  }
  // case array not full
  if(index < size){
    //case element should be placed after the last valid value
    if(my_arr[index-1].value <= value){
      my_arr[index].index = index;
      my_arr[index].value = value;
      return;
    }
    //case the element should be placed between 0 and index-1
    //Rcout<<"Calling binary 0 "<<index-1<<" val "<<value<<endl;
    int i;
    for(i = index-1; i>=0;i--){
      if(my_arr[i].value <= value){
        break;
	  }
	}
    //value should be written in i
    for(int j = index; j > i+1;j--){
      my_arr[j] = my_arr[j-1];
	}

    my_arr[i+1].value = value;
    my_arr[i+1].index = index;
    return;

  }
  else{
    int i;

    for(i = size-1; i>=0;i--){
      if(my_arr[i].value <= value){
        break;
	  }
	}
        //value should be written in i
      for(int j = size-1; j > i+1;j--){
        my_arr[j] = my_arr[j-1];
	  }

      my_arr[i+1].value = value;
      my_arr[i+1].index = index;

      return;
  }

  return;
}

a_node* init_array(int K){
  a_node* my_ar = new a_node[K];
  a_node* first = my_ar;
  for(int i=0;i<K;i++,my_ar++){
    my_ar->index = -1;
  }

  if(K<100){
    _place_vals_ = &place_new_values2;
  }
  else{
    _place_vals_ = &place_new_values;
  }

  return first;
}

void k_sorted_put(a_node* my_arr, int size,int index, double value){
  _place_vals_(my_arr,size,index,value);
}
