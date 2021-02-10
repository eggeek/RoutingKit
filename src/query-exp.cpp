#include <routingkit/vector_io.h>
#include <routingkit/permutation.h>
#include <routingkit/inverse_vector.h>
#include <routingkit/customizable_contraction_hierarchy.h>
#include <routingkit/nested_dissection.h>
#include <routingkit/graph_util.h>
#include <routingkit/min_max.h>
#include <routingkit/timer.h>
#include <routingkit/permutation.h>

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <vector>
#include <random>
#include <map>

using namespace RoutingKit;
using namespace std;

double 
fRand(double fMin, double fMax)
{
  double f = (double)rand() / RAND_MAX;
  return fMin + f * (fMax - fMin);
}

struct QueryRow
{
  string alg;
  unsigned id, pcost, plen, custnum;
  double ref_dist, tcost, custcost, extcost;
  string header = "id\talg\tnanos\textcost\tpcost\tplen\tref_dist\tcustcost\tcustnum";

  string to_string()
  {
    string res = "";
    res += std::to_string(id) + "\t";
    res += alg + "\t";
    res += std::to_string(tcost) + "\t";
    res += std::to_string(extcost) + "\t";
    res += std::to_string(pcost) + "\t";
    res += std::to_string(plen) + "\t";
    res += std::to_string(ref_dist) + "\t";
    res += std::to_string(custcost) + "\t";
    res += std::to_string(custnum);
    return res;
  }
};

vector<pair<unsigned, unsigned>>
load_w_change(string custfile, 
    vector<unsigned>& tail, 
    vector<unsigned>& head,
    vector<unsigned>& weight)
{
  map<pair<unsigned, unsigned>, unsigned> idx;

  for (int i=0; i<(int)tail.size(); i++)
    idx[{tail[i], head[i]}] = i;

  vector<double> r;

  vector<pair<unsigned, unsigned>> res;
  ifstream ifs(custfile);
  int num;
  ifs >> num;
  res.resize(num);
  if (num == 0) return res;

  for (int i=0; i<num; i++)
  {
    unsigned u_, v_;
    double w_;
    ifs >> u_ >> v_ >> w_;
    unsigned eid = idx[{u_, v_}];
    res[i] = {eid, (unsigned)w_};
    r.push_back((double) w_ / (double)weight[eid]);
  }
  sort(r.begin(), r.end());
  cerr
    << "changed ratio: "
    << r[0] << "\t"
    << r[r.size() * 0.01] << "\t"
    << r[r.size() / 2] << "\t"
    << r[r.size() * 0.75] << "\t"
    << r.back() << endl;
  return res;
}



vector<unsigned> first_out, tail, head, weight, order, source, target;
vector<double> ref_dist;
vector<float> lat, lng;
map<pair<unsigned, unsigned>, unsigned> mat;

void load_queries(string qfile)
{
  ifstream fin(qfile);
  string line, desc;
  int num, id = 0;
  while (getline(fin, line))
  {
    if (line.empty() || line[0] == '#' || line[0] == 'c')
      continue;
    stringstream sin(line);
    sin >> desc;
    if (desc == "p")
    {
      sin >> desc >> desc >> desc >> num >> desc;
      source.resize(num);
      target.resize(num);
      ref_dist.resize(num);
      id = 0;
    }
    else if (desc == "q")
    {
      sin >> source[id] >> target[id] >> ref_dist[id];
      id++;
    }
  }
}

bool sanity_path_checking(vector<unsigned>& path, unsigned dist)
{
  unsigned tot = 0;
  for (int i=1; i<(int)path.size(); i++)
  {
    if (mat.find({path[i-1], path[i]}) == mat.end())
    {
      cerr << "invalid edge on path" << endl;
      return false;
    }
    tot += mat[{path[i-1], path[i]}];
  }
  if (tot != dist)
  {
    cerr << "path length disagree with dist" << endl;
    return false;
  }
  return true;
}

void run_experiment(CustomizableContractionHierarchyMetric& metric, string alg, int custnum, double custcost, int rep=10, int start=0, int end=-1)
{
  CustomizableContractionHierarchyQuery query(metric);

  vector<unsigned> path;
  long long tot = 0, ext = 0;
  vector<QueryRow> rows;
  int qnum = source.size();
  rows.resize(qnum);
  if (end == -1) end = qnum;
  cerr << "Running queries from " << start << " to " << end << " in average " << rep << " runs" << endl;
  for (int cnt=0; cnt<rep; cnt++)
  {
    for (int i=start; i<end; i++)
    {
      unsigned u = source[i], v = target[i];
      tot = -get_micro_time();
      query.reset().add_source(u).add_target(v).run();

      unsigned dist = query.get_distance();

      ext = -get_micro_time();
      path = query.get_node_path();
      ext += get_micro_time();

      tot += get_micro_time();

      QueryRow& r = rows[i];
      r.id = i, r.alg = alg;
      r.pcost = dist, r.plen = path.size();
      r.ref_dist = ref_dist[i];
      if (cnt == 0) {
        r.tcost = tot;
        r.extcost = ext;
      }
      else {
        r.tcost += tot;
        r.extcost += ext;
      }
    }
  }
  for (int i=start; i<end; i++)
  {
    rows[i].tcost /= (double)rep;
    rows[i].extcost /= (double)rep;
    // micro to nanos
    rows[i].tcost *= 1000;
    rows[i].extcost *= 1000;
    rows[i].custcost = custcost * 1000;
    rows[i].custnum = custnum;
    cout << rows[i].to_string() << endl;
  }
}

void CCH()
{
  CustomizableContractionHierarchy cch(order, tail, head, [](const std::string&){}, true);
  CustomizableContractionHierarchyMetric reference_metric(cch, weight);
  long long tot = 0, t;
  for (int i = 0; i < 10; i++) {
    reference_metric.reset(weight);
    t = -get_micro_time();
    reference_metric.customize();
    t += get_micro_time();
    cerr << "Customization in " << t << " micros" << endl;
    tot += t;
  }


  cerr << "Running queries ..." << endl;
  QueryRow row;
  cout << row.header << endl;
  run_experiment(reference_metric, "cch", 0, tot / 10);
}

void CCH_perturb(string diff_file)
{
  auto c = load_w_change(diff_file, tail, head, weight);
  CustomizableContractionHierarchy cch(order, tail, head, [](const std::string&){}, true);
	CustomizableContractionHierarchyMetric partial_metric(cch, weight);
  CustomizableContractionHierarchyPartialCustomization partial_update(cch);

  long long tot = 0, t;
  for (int r=0; r<10; r++) {
    partial_metric.reset(weight);
    partial_update.reset();
    partial_metric.customize();
    cerr << "update " << c.size() << " arcs, ";
    for (unsigned i=0; i<c.size(); i++)
    {
      weight[c[i].first] = c[i].second;
      partial_update.update_arc(c[i].first);
    }

    t = -get_micro_time();
    partial_update.customize(partial_metric);
    t += get_micro_time();
    cerr << "done in " << t << " micros" << endl;
    tot += t;
  }

  cerr << "Running queries ..." << endl;
  QueryRow row;
  cout << row.header << endl;
  run_experiment(partial_metric, "cch_perturb", c.size(), tot / 10);
}

void CCH_random_perturb(double ratio=0.1)
{
  CustomizableContractionHierarchy cch(order, tail, head, [](const std::string&){}, true);
	CustomizableContractionHierarchyMetric partial_metric(cch, weight);
  CustomizableContractionHierarchyPartialCustomization partial_update(cch);

  long long t;
  int num = weight.size() * ratio;
  cerr << "update " << num << " arcs, ";
  partial_metric.customize();
  for (int i=0; i<num; i++)
  {
    unsigned id = rand() % weight.size();
    weight[id] = (int)(weight[id] * fRand(0.3, 1.5));
    partial_update.update_arc(id);
  }
  t = -get_micro_time();
  partial_update.customize(partial_metric);
  t += get_micro_time();
  cerr << "done in " << t << " micros" << endl;

  cerr << "Running queries ..." << endl;
  QueryRow row;
  cout << row.header << endl;
  run_experiment(partial_metric, "cch_rd_cust", num, t);
}

void CCH_cumulative_perturb(string diff_file)
{
  auto c = load_w_change(diff_file, tail, head, weight);
  int bnum = 10;
  int bsize = c.size() / bnum;

  CustomizableContractionHierarchy cch(order, tail, head, [](const std::string&){}, true);
	CustomizableContractionHierarchyMetric partial_metric(cch, weight);
  CustomizableContractionHierarchyPartialCustomization partial_update(cch);
  partial_metric.customize();

  QueryRow row;
  cout << row.header << endl;

  int totq = 0;
  for (int i=0; i<bnum; i++)
  {
    int lb = i * bsize, ub = i<bnum-1?(i+1)*bsize: c.size();
    partial_update.reset();
    for (int j=lb; j<ub; j++)
    {
      partial_update.update_arc(c[j].first);
      weight[c[j].first] = c[j].second;
    }
    long long t = -get_micro_time();
    partial_update.customize(partial_metric);
    t += get_micro_time();
    int qnums = min((1 << i) * 100, (int)source.size() - totq);
    run_experiment(partial_metric, "cch_cumulative-" + to_string(i), ub-lb, t, 10, totq, totq + qnums);
    totq += qnums;
  }
}

void init_mat()
{
  mat.clear();
  for (int i=0; i<(int)tail.size(); i++)
  {
    mat[{tail[i], head[i]}] = weight[i];
  }
}

int main(int argc, char* argv[])
{
  string first_out_file, head_file, weight_file, diff_file, query_file, opt = "normal";
  string order_file;
  if (argc < 7) {
    cerr << argv[0] 
         << "first_out head weight order query diff" << endl;
    exit(1);
  }
  else {
    first_out_file = argv[1];
    head_file = argv[2];
    weight_file = argv[3];
    order_file = argv[4];
    query_file = argv[5];
    diff_file = argv[6];
    if (argc >= 8)
      opt = argv[7];
  }

  first_out = load_vector<unsigned>(first_out_file);
  // unsigned node_count = first_out.size() - 1;
	tail = invert_inverse_vector(first_out);
  head = load_vector<unsigned>(head_file);
  weight = load_vector<unsigned>(weight_file);
  order = load_vector<unsigned>(order_file);
  init_mat();
  load_queries(query_file);
  if (opt == "normal") 
    CCH();
  else if (opt == "diff")
    CCH_perturb(diff_file);
  else if (opt == "random")
    CCH_random_perturb();
  else if (opt == "cumulative")
    CCH_cumulative_perturb(diff_file);
}
