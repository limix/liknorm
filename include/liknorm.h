#ifndef LIKNORM_H
#define LIKNORM_H

typedef void log_partition (double x,
                            void* lp_data,
                            double *A0,
                            double *logA1,
                            double *logA2,
                            double* sign);

typedef struct
{
        double Ty;
        log_partition* lp;
        void* lp_data;
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

#endif /* ifndef LIKNORM_H */
