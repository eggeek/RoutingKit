#include <sys/stat.h>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <csignal>
#include <omp.h>
#include <fstream>
#include <iostream>

#include <routingkit/vector_io.h>
#include <routingkit/permutation.h>
#include <routingkit/inverse_vector.h>
#include <routingkit/customizable_contraction_hierarchy.h>
#include <routingkit/nested_dissection.h>
#include <routingkit/graph_util.h>
#include <routingkit/min_max.h>
#include <routingkit/timer.h>
#include <routingkit/permutation.h>
#include "cfg.h"
#include "log.h"
#include "warthog_timer.h"
#include "json_config.h"
#include "json.h"
using namespace RoutingKit;
using namespace std;
// get the index of an edge <u, v>
map<pair<unsigned, unsigned>, unsigned> idx;
vector<unsigned> first_out, tail, head, weight, order;
vector<CustomizableContractionHierarchyQuery> algos;
string first_out_file, head_file, weight_file, order_file;
string fifo = "/tmp/cchfifo";
CustomizableContractionHierarchy cch;
CustomizableContractionHierarchyMetric metric;
CustomizableContractionHierarchyPartialCustomization partial_update;
warthog::util::cfg cfg;

void
signalHandler(int signum)
{
    warning(true, "Interrupt signal", signum, "received.");

    remove(fifo.c_str());

    exit(signum);
}

struct Query {
  int s, t;
};

struct Perturb {
  unsigned u, v, w;
};

double run_customization(config& conf, vector<Perturb>& diffs) {
  trace(conf.verbose, "Preparing to perturb", diffs.size(), "edges.");
  warthog::timer t;
  double custnano = 0;
  t.start();
  partial_update.reset();
  for (auto& diff: diffs) {
    auto eid = idx[{diff.u, diff.v}];
    weight[eid] = diff.w;
    partial_update.update_arc(eid);
  }
  partial_update.customize(metric);
  t.stop();
  custnano = t.elapsed_time_nano();
  trace(conf.verbose, "Processed", diffs.size(), "edges in ", t.elapsed_time_nano(), "us.");
  return custnano;
}

void run_experiment(config& conf,
    vector<Query>& queries, string& fifo_out,
    const double& custnano, const double& pertubnano, const double& t_read) {

  size_t n_results = queries.size();
  vector<unsigned> path;
  warthog::timer t, tquery;
  long long t_query = 0, t_ext = 0;
  streambuf* buf;
  ofstream of;
  if (fifo_out == "-")
  {
      buf = std::cout.rdbuf();
  }
  else
  {
      of.open(fifo_out);
      buf = of.rdbuf();
  }
  ostream out(buf);

#ifdef SINGLE_THREADED
      unsigned int threads = 1;
#else
      unsigned int threads = conf.threads;
#endif


  user(conf.verbose, "Preparing to process", queries.size(), "queries using",
       (int)threads, "threads.");

  t.start();

#pragma omp parallel num_threads(threads)                               \
    reduction(+ : t_query, t_ext)
    {

      // Parallel data
      unsigned int thread_count = omp_get_num_threads();
      unsigned int thread_id = omp_get_thread_num();
      warthog::timer t_thread;

      size_t from = 0;
      size_t to = n_results;
      if (!conf.thread_alloc)
      {
          // Instead of bothering with manual conversion (think 'ceil()'), we
          // use the magic of "usual arithmetic" to achieve the right from/to
          // values.
          size_t step = n_results * thread_id;
          from = step / thread_count;
          to = (step + n_results) / thread_count;
      }

      CustomizableContractionHierarchyQuery& algo = algos[thread_id];
      t_thread.start();
      for (size_t i=from; i<to; i++) {
        const Query& q = queries[i];
        tquery.start();
        algo.reset().add_source(q.s).add_target(q.t).run();
        auto d = algo.get_distance();
        tquery.stop();

        t_query += tquery.elapsed_time_nano();

        tquery.start();
        path = algo.get_node_path();
        tquery.stop();
        t_query += tquery.elapsed_time_nano();
        t_ext += tquery.elapsed_time_nano();
        debug(conf.verbose, "[", thread_id, "]: query:", i, "from", q.s, "to", q.t, "dist:", d);
      }
      t_thread.stop();

#pragma omp critical
          trace(conf.verbose, "[", thread_id, "] Processed", to - from,
                "trips in", t_thread.elapsed_time_nano(), "us.");
    }
    t.stop();
    user(conf.verbose, "Processed", n_results, "in", t.elapsed_time_nano(), "us.");
  out << (long long)t_query << ","  // time cost on query
      << (long long)t_ext << "," // time cost on path extraction
      << (long long)t.elapsed_time_nano() << "," // time cost on entire function
      << (long long)custnano << "," // time cost on customize
      << (long long)pertubnano << "," // time cost on perturbation
      << (long long)t_read // time cost on reading
      << endl;
  if (fifo_out != "-") { of.close(); }
}

void load_queries(string& qfile, config& conf, vector<Query>& queries) {
  if (qfile == "-") {
    debug(conf.verbose, "No query file, skip.");
    return;
  }
  debug(conf.verbose, "Reading queries from", qfile);
  size_t s = 0;
  ifstream fd(qfile);
  fd >> s;
  debug(conf.verbose, "Preparing to read", s, "items.");
  queries.resize(s);
  for (size_t i=0; i<s; i++) 
    fd >> queries[i].s >> queries[i].t;
  trace(conf.verbose, "Read", queries.size(), "queries.");
  fd.close();
}

void load_diff(string& dfile, config& conf, vector<Perturb>& diffs) {
  if (dfile == "-") {
    debug(conf.verbose, "No diff file, skip.");
    return;
  }
  debug(conf.verbose, "Reading diff from", dfile);
  size_t s = 0;
  ifstream fd(dfile);
  fd >> s;
  debug(conf.verbose, "Preparing to read", s, "perturbations.");
  diffs.resize(s);
  for (size_t i=0; i<s; i++) {
    fd >> diffs[i].u >> diffs[i].v >> diffs[i].w;
  }
  trace(conf.verbose, "Read", s, "perturbations.");
  fd.close();
}


void run_cch() {
  warthog::timer t;
  fstream fd;
  config conf;
  string fifo_out, qfile, dfile;
  vector<Query> queries;
  vector<Perturb> diffs;
  double custnano = 0;
  double pertubnano = 0;
  cch = CustomizableContractionHierarchy(order, tail, head, [](const std::string&){}, true);
	metric = CustomizableContractionHierarchyMetric(cch, weight);
  partial_update = CustomizableContractionHierarchyPartialCustomization(cch);

  t.start();
  metric.customize();
  t.stop();
  custnano = t.elapsed_time_nano();

  for (auto& algo: algos) {
    algo = CustomizableContractionHierarchyQuery(metric);
  }

  debug(true, "Customization:", t.elapsed_time_nano(), " us");
  while (true) {
    fd.open(fifo);
    debug(VERBOSE, "waiting for writers...");

    if (fd.good())
    {
        debug(VERBOSE, "Got a writer");
    }

    t.start();
    // Start by reading config
    try
    {
        fd >> conf;
    } // Ignore bad parsing and fall back on default conf
    catch (std::exception& e)
    {
        debug(conf.verbose, e.what());
    }
    trace(conf.verbose, conf);
    fd >> qfile >> dfile >> fifo_out;
    debug(conf.verbose, "Output to", fifo_out);
    load_queries(qfile, conf, queries);
    load_diff(dfile, conf, diffs);
    t.stop();
    if (diffs.size()) {
      pertubnano = run_customization(conf, diffs);
    }
    if (queries.size()) {
      run_experiment(conf, queries, fifo_out, custnano, pertubnano, t.elapsed_time_nano());
    }
  }
}


int main(int argc, char* argv[]) {

  // run command
  // ./fifo --input <first_out_file> [head_file] [weight_file] [order_file] [--fifo fifo name]
  // parse arguments
  warthog::util::param valid_args[] =
      {
          {"input", required_argument, 0, 1},
          {"fifo",  required_argument, 0, 1},
          {"alg",   required_argument, 0, 1},
          {"div",   required_argument, 0, 1},
          {"mod",   required_argument, 0, 1},
          {"num",   required_argument, 0, 1},
          // {"problem",  required_argument, 0, 1},
          {0,  0, 0, 0}
      };

  cfg.parse_args(argc, argv, "-f", valid_args);
  first_out_file = cfg.get_param_value("input");
  if (first_out_file == "") {
    std::cerr << "parameter is missing: --input\n";
    return EXIT_FAILURE;
  }
  head_file = cfg.get_param_value("input");
  if (head_file == "") {
    head_file = first_out_file.substr(0, first_out_file.find_last_of(".")) + ".head";
  }
  weight_file = cfg.get_param_value("input");
  if (weight_file== "") {
    weight_file = first_out_file.substr(0, first_out_file.find_last_of(".")) + ".weight";
  }
  order_file = cfg.get_param_value("input");
  if (order_file == "") {
    order_file = first_out_file.substr(0, first_out_file.find_last_of(".")) + ".cch_order";
  }

#ifdef SINGLE_THREADED
    algos.resize(1);
#else
    algos.resize(omp_get_max_threads());
#endif

  string other = cfg.get_param_value("fifo");
  if (other != "") {
    fifo = other;
  }

  // create fifo
  int status = mkfifo(fifo.c_str(), S_IFIFO | 0666);

  if (status < 0)
  {
      perror("mkfifo");
      return EXIT_FAILURE;
  }

  debug(true, "Reading from", fifo);

  // Register signal handlers
  signal(SIGINT, signalHandler);
  signal(SIGTERM, signalHandler);
  signal(SIGABRT, signalHandler);

  first_out = load_vector<unsigned>(first_out_file);
	tail = invert_inverse_vector(first_out);
  head = load_vector<unsigned>(head_file);
  weight = load_vector<unsigned>(weight_file);
  order = load_vector<unsigned>(order_file);
  // init index
  for (size_t i=0; i<tail.size(); i++) {
    idx[{tail[i], head[i]}] = i;
  }
  run_cch();
}
