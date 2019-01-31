function(define_public_headers BASE_DIR _PUBLIC_HEADERS _PRIVATE_HEADERS)
    set(${_PUBLIC_HEADERS} "${BASE_DIR}/liknorm.h"
        PARENT_SCOPE)
endfunction(define_public_headers)
