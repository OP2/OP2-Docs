// Kernel computing the center of mass of an element
void kernel( const double x[3][2], double c[2]) {
  c[0] = (x[0][0] + x[1][0] + x[2][0])/3.0f;
  c[1] = (x[0][1] + x[1][1] + x[2][1])/3.0f;
}

// OP2 main program
#include <op_lib_cpp.h>
#include <op_seq.h>

int main(int argc, char **argv) {

  // Vertex coordinates
  double x[] = { 0.0f, 0.0f,
                 0.9f, 0.1f,
                 0.1f, 0.9f,
                 1.0f, 1.0f };

  // Center of mass (initialised to zero)
  double c[] = { 0.0f, 0.0f,
                 0.0f, 0.0f };

  // Element to vertex mapping
  int p1_ele_i[] = { 1,2,3,
                     3,2,4 };

  // OP2 initialisation
  op_init(argc, argv, 2);

  // Declare sets, maps, and datasets
  op_set p1       = op_decl_set(4, "nodes");
  op_set elements = op_decl_set(2, "elements");

  op_map p1_ele = op_decl_map(elements, p1, 3, p1_ele_i, "element-node");

  op_dat coordinates = op_decl_dat(p1,       2, "double", x, "coordinates");
  op_dat com         = op_decl_dat(elements, 2, "double", c, "com");

  // Parallel loop over the elements
  op_par_loop ( kernel, "kernel", elements,
                op_arg_dat(coordinates, OP_ALL, p1_ele, 2, "double", OP_READ),
                op_arg_dat(com,         -1,     OP_ID,  2, "double", OP_INC)
              );
}

