#include "../src/temporal.h"

namespace shimmer{

variable::variable(){};
variable::variable(const vector_t&p,const vector_t&f,const vector_t&l)
    {
        pressure = p;
        flux = f;
        L_rate = l;
    };


vector_t
variable::make_vector() const
{
    size_t num_pipes = flux.size();
    size_t num_nodes = pressure.size();

    vector_t vec(2 * num_nodes + num_pipes);
    vec.head(num_nodes) = pressure;
    vec.segment(num_nodes, num_pipes) = flux;
    vec.tail(num_nodes) = L_rate;

    return vec;
}

}