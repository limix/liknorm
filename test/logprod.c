#include "cass.h"
#include "liknorm.h"
#include "stdio.h"

int main() {

    struct LikNormMachine *machine = liknorm_create_machine(10);
    double y = 5.0;
    double tau = 0.15;
    double eta = 2.1;
    liknorm_set_poisson(machine, y);
    liknorm_set_prior(machine, tau, eta);

    assert_almost_equal(liknorm_logprod(machine, 700),
                        -1.0142320547350045e+304);
    assert_almost_equal(liknorm_logprod(machine, 50), -5.184705528587072e+21);
    assert_almost_equal(liknorm_logprod(machine, 12.0), -162701.74640927234);
    assert_almost_equal(liknorm_logprod(machine, 0.59), -18.996086183827522);
    assert_almost_equal(liknorm_logprod(machine, -3.0), -43.37977733679752);
    assert_almost_equal(liknorm_logprod(machine, -10.9), -107.65575872666366);
    assert_almost_equal(liknorm_logprod(machine, -20), -193.35499027049082);
    assert_almost_equal(liknorm_logprod(machine, -50), -563.8549902684297);
    assert_almost_equal(liknorm_logprod(machine, -500), -22321.35499026843);

    double k = 1;
    double n = 10;
    liknorm_set_binomial(machine, k, n);
    assert_almost_equal(liknorm_logprod(machine, 0.0), -21.196385238253026);
    assert_almost_equal(liknorm_logprod(machine, 1.0), -24.372530307835795);
    assert_almost_equal(liknorm_logprod(machine, 10.0), -90.76536742165445);
    assert_almost_equal(liknorm_logprod(machine, 50.0), -546.7649134326535);
    assert_almost_equal(liknorm_logprod(machine, 100.0), -1454.2649134326534);
    assert_almost_equal(liknorm_logprod(machine, 500.0), -22214.264913432653);
    assert_almost_equal(liknorm_logprod(machine, 5000.0), -1909514.2649134325);

    assert_almost_equal(liknorm_logprod(machine, -10.0), -52.76536742164575);
    assert_almost_equal(liknorm_logprod(machine, -50.0), -356.7649134326536);
    assert_almost_equal(liknorm_logprod(machine, -100.0), -1074.2649134326537);
    assert_almost_equal(liknorm_logprod(machine, -500.0), -20314.264913432653);
    assert_almost_equal(liknorm_logprod(machine, -5000.0), -1890514.2649134325);

    liknorm_destroy_machine(machine);

    return 0;
}
