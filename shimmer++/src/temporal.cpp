#include "../src/temporal.h"

namespace shimmer{

variable::variable(){};
variable::variable(const vector_t&p,const vector_t&f,const vector_t&l)
    {
        pressure = p;
        flux = f;
        L_rate = l;
    };
}