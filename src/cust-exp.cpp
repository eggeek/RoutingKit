#include <routingkit/vector_io.h>
#include <routingkit/inverse_vector.h>
#include <routingkit/customizable_contraction_hierarchy.h>
#include <routingkit/timer.h>

#include <vector>
#include <map>
#include <iostream>

using namespace RoutingKit;
using namespace std;

vector<pair<unsigned, unsigned>>
load_w_change(string custfile, 
    vector<unsigned>& tail, 
    vector<unsigned>& head,
    vector<unsigned>& weight,
    double &ratio)
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
  ratio = r[r.size() / 2];
  cerr
    << "changed ratio: "
    << r[0] << "\t"
    << r[r.size() * 0.01] << "\t"
    << r[r.size() / 2] << "\t"
    << r[r.size() * 0.75] << "\t"
    << r.back() << endl;
  return res;
}

vector<unsigned> first_out, tail, head, weight, order;
string first_out_file, head_file, weight_file, diff_file, order_file;


long long normal_cust()
{
  double ratio;
  auto weight = load_vector<unsigned>(weight_file);
  auto c = load_w_change(diff_file, tail, head, weight, ratio);
  for (auto it: c) weight[it.first] = it.second;

  CustomizableContractionHierarchy cch(order, tail, head, [](const std::string&){}, true);
	CustomizableContractionHierarchyMetric metric(cch, weight);

  long long t = -get_micro_time();
  metric.customize();
  t += get_micro_time();
  return t;
}

long long partial_cust()
{
  double ratio;
  auto weight = load_vector<unsigned>(weight_file);
  auto c = load_w_change(diff_file, tail, head, weight, ratio);
  CustomizableContractionHierarchy cch(order, tail, head, [](const std::string&){}, true);
	CustomizableContractionHierarchyMetric partial_metric(cch, weight);
  CustomizableContractionHierarchyPartialCustomization partial_update(cch);
  partial_metric.customize();
  for (unsigned i=0; i<c.size(); i++)
  {
    weight[c[i].first] = c[i].second;
    partial_update.update_arc(c[i].first);
  }

  long long t = -get_micro_time();
  partial_update.customize(partial_metric);
  t += get_micro_time();
  return t;
}

long long seq_cust(int nruns)
{
  double ratio;
  auto weight = load_vector<unsigned>(weight_file);
  auto c = load_w_change(diff_file, tail, head, weight, ratio);
  CustomizableContractionHierarchy cch(order, tail, head, [](const std::string&msg){cout << msg << endl;}, true);
	CustomizableContractionHierarchyMetric partial_metric(cch, weight);
  CustomizableContractionHierarchyPartialCustomization partial_update(cch);
  partial_metric.customize();

  int bsize = c.size() / nruns;
  long long tot = 0;

  for (int i=0; i<nruns-1; i++) {

    // int i = 0;
    partial_update.reset();
    for (int j=i*bsize; j<(i+1)*bsize; j++) {
      weight[c[j].first] = c[j].second;
      partial_update.update_arc(c[j].first);
    }
    long long t = -get_micro_time();
    partial_update.customize(partial_metric);
    t += get_micro_time();

    tot += t;
  }

  partial_update.reset();
  for (int j=(nruns-1)*bsize; j<(int)c.size(); j++) {
    weight[c[j].first] = c[j].second;
    partial_update.update_arc(c[j].first);
  }
  long long t = -get_micro_time();
  partial_update.customize(partial_metric);
  t += get_micro_time();
  tot += t;
  return tot;
}

void run_customization()
{
  long long t_normal = normal_cust();
  long long t_partial = partial_cust();
  long long t_seq = seq_cust(10);

  auto weight = load_vector<unsigned>(weight_file);
  double ratio;
  auto c = load_w_change(diff_file, tail, head, weight, ratio);

  string header = "perturb_num\tnormal\tpartial\tsequential\tratio";
  cout << header << endl;
  cout << c.size() << "\t"
       << t_normal << "\t"
       << t_partial << "\t"
       << t_seq << "\t" 
       << diff_file << endl;
}

int main(int argc, char** argv)
{
  if (argc != 6) {
    cerr << argv[0] 
         << "first_out head weight order diff" << endl;
    exit(1);
  }
  else {
    first_out_file = argv[1];
    head_file = argv[2];
    weight_file = argv[3];
    order_file = argv[4];
    diff_file = argv[5];
  }
  first_out = load_vector<unsigned>(first_out_file);
	tail = invert_inverse_vector(first_out);
  head = load_vector<unsigned>(head_file);
  order = load_vector<unsigned>(order_file);
  run_customization();
  return 0;
}
