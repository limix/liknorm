function(define_sources BASE_DIR _SOURCES)

    list(APPEND REL_SOURCES normal.c
        optimizer/zero.c
        optimizer/brent.c
        optimizer/bracket.c
        interval.c
        integrate.c
        gfunc.c
        partition/gamma.c
        partition/geometric.c
        partition/exponential.c
        partition/binomial.c
        partition/poisson.c
        liknorm.c)

    limix_convert_rel_to_full(FULL_SOURCES ${BASE_DIR} "${REL_SOURCES}")
    set(${_SOURCES} ${FULL_SOURCES} PARENT_SCOPE)

endfunction(define_sources)
