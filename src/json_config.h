#ifndef __JSON_CONFIG_H_
#define __JSON_CONFIG_H_

#ifndef NDEBUG
#define VERBOSE true
#else
#define VERBOSE false
#endif

#include <omp.h>
#include "constants.h"
#include "json.h"

//
// - Definitions
//
typedef std::tuple<unsigned int, // Nodes expanded
                   unsigned int, // Nodes inserted
                   unsigned int, // Nodes touched
                   unsigned int, // Nodes updated
                   unsigned int, // Surplus nodes
                   unsigned int, // Total length of paths
                   double,       // Difference between perturbed and opti
                   double,       // Difference between free flow and opti
                   double,       // Difference between free flow and perturbed
                   double,       // Time spent in A*
                   double>       // Time to do the search
    t_results;

typedef struct config
{
    double time = DBL_MAX;
    uint32_t itrs = warthog::INF32;
    unsigned char threads = 0;
    bool verbose = VERBOSE;
    bool debug = false;
    bool thread_alloc = false;
} config;

void
to_json(nlohmann::json& j, const config& c)
{
    j = {
        {"time", c.time},
        {"itrs", c.itrs}, {"threads", c.threads},
        {"verbose", c.verbose}, {"debug", c.debug},
        {"thread_alloc", c.thread_alloc}
    };
}

void
from_json(const nlohmann::json& j, config &c)
{
    j.at("time").get_to(c.time);
    j.at("itrs").get_to(c.itrs);
    j.at("threads").get_to(c.threads);
    j.at("verbose").get_to(c.verbose);
    j.at("debug").get_to(c.debug);
    j.at("thread_alloc").get_to(c.thread_alloc);
}

std::ostream&
operator<<(std::ostream& os, config &c)
{
    nlohmann::json j = c;

    os << j;

    return os;
}

config&
operator>>(std::istream& is, config &c)
{
    nlohmann::json j;

    is >> j;
    c = j.get<config>(); // this will only read a complete object

    return c;
}

/**
 * Takes care of "default parameters" as we use a bunch of wildcards to
 * represent different unbounded values.
 *
 * TODO Should this be part of cpd search directly?
 */
inline void
sanitise_conf(config& conf)
{

    if (conf.itrs == 0)
    { conf.itrs = warthog::INF32; }

    if (conf.time == 0)
    { conf.time = DBL_MAX; }

    // Enforce single threaded or use max threads
#ifdef SINGLE_THREADED
    conf.threads = 1;
#else
    if (conf.threads == 0)
    {
        conf.threads = omp_get_max_threads();
    }
#endif

}


#endif // __JSON_CONFIG_H_
