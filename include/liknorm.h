#ifndef LIKNORM_H
#define LIKNORM_H

typedef void log_partition (double  theta,
                            double *b0,
                            double *logb1,
                            double *logb2,
                            double *sign);

typedef struct
{
        double y;
        double aphi;
        log_partition* lp;
} ExpFam;

typedef struct
{
        double eta;
        double tau;
} Normal;

typedef struct
{
        double *log_zeroth;
        double *u;
        double *v;
        int n;
        double precision;
} LikNormMachine;

LikNormMachine* create_liknorm_machine(int n, double precision);
void destroy_liknorm_machine(LikNormMachine* machine);
void integrate(LikNormMachine *machine, ExpFam ef, Normal normal,
               double* mean, double *variance);
log_partition* get_log_partition(char *name);

#endif /* ifndef LIKNORM_H */
