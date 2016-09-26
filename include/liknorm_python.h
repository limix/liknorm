void py_create_liknorm_machine(int n, double precision);
void py_destroy_liknorm_machine();
void py_integrate(int likname_id, double *y, double *aphi,
                double *normal_tau, double *normal_eta, int n,
                double *mean, double *variance);
