#include "liknorm/liknorm.h"

int main()
{
    struct LikNormMachine *machine = liknorm_create_machine(1000);
    liknorm_destroy_machine(machine);

    return 0;
}
