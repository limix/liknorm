#include "liknorm.h"

int main()
{
  LikNormMachine *machine = liknorm_create_machine(1000);
  liknorm_destroy_machine(machine);

  return 0;
}
