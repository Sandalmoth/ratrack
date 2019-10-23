#include <cassert>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <regex>
#include <string>
#include <tclap/CmdLine.h>
#include <valarray>
#include <vector>


using namespace std;


const string VERSION = "0.0.1";


struct Arguments {
  int n0;
  vector<double> birth_rate;
  double death_rate;
  double interaction_birth_rate;
  double interaction_death_rate;
  double t_end;
};


class Cell {
public:
  Cell(vector<double> b2, double d, double p, double q, double t_max)
    : b_base(b2)
    , d(d)
    , p(p)
    , q(q)
    , t_max(t_max) {


    // TODO Switch to monotonic cubic interpolation
    // cout << "CELL CONSTRUCTOR" << endl;


    // cout << "base vector " << b.size() << " \n";
    // for (auto x: b_base) cout << x << endl;

    // std::cout.precision(3);

    for (int z = 0; z < t_max * INTERPOLATIONS; ++z) {
      size_t i = static_cast<double>(z) / static_cast<double>(INTERPOLATIONS) / t_max * (b_base.size() - 1);
      double x = static_cast<double>(z) / static_cast<double>(INTERPOLATIONS) / t_max * (b_base.size() - 1) - i;
      b.push_back(b_base[i] * (1 - x) + b_base[i+1] * x);
      // cout << z << '\t' << i << '\t' << x << '\t' << b.back() << endl;
    }
    b.push_back(b_base.back());

    // cout << "interpolation vector " << b.size() << " \n";
    // for (auto x: b) cout << x << endl;
    // cout << "CELL CONSTRUCTOR END" << endl;
  }


  // const Cell &operator=(const Cell &other) {
  // }

  double get_birth_rate(double t) {
    // cout << endl << static_cast<size_t>(t*INTERPOLATIONS) << '\t' << b.size() << endl;
    return b[static_cast<size_t>(t*INTERPOLATIONS)];
  }

  const int INTERPOLATIONS = 100.0;
  vector<double> b_base;
  vector<double> b;
  double d;
  double p;
  double q;
  double t_max;
};


template <typename TCell, typename TRng=std::mt19937>
class LB {
public:
  LB(TCell wt) {
    urd = std::uniform_real_distribution<double>(std::nextafter(0.0, 1.0), 1.0);
    std::random_device rd;
    rng.seed(rd());
    type_count = 1;
    X.resize(type_count);
    a.resize(type_count * 2);
    T.resize(type_count * 2);
    P.resize(type_count * 2);
    r.resize(type_count * 2);
    dt.resize(type_count * 2);
    X = 0;
    // a has to be calculated, so no reason to waste time initializing
    T = 0;
    P = 0;
    for (size_t i = 0; i < type_count; ++i) {
      cells.push_back(wt);
    }
  }

  void set_cell_count(size_t count) {
    X[0] = count;
  }

  void get_birth_rates() {
    for (size_t i = 0; i < type_count; ++i) {
      double rate = cells[i].get_birth_rate(t);
      double interaction = cells[i].p * (cells[i].get_birth_rate(t) - cells[i].d); // Constant carrying capacity
      // double interaction = cells[i].p; // Carrying capacity varies with effective reproduction rate
      int sizemult = X.sum() - 1;
      double ai = X[i] * rate - sizemult * X[i] * interaction;
      if (ai >= 0.0)
        a[i*2] = ai;
      else {
        a[i*2] = 0;
        a[i*2 + 1] -= ai;
      }
    }
  }
  void get_death_rates() {
    for (size_t i = 0; i < type_count; ++i) {
      double rate = cells[i].d;
      double interaction = cells[i].q * (cells[i].get_birth_rate(t) - cells[i].d); // Constant carrying capacity
      // double interaction = cells[i].q; // Carrying capacity varies with effective reproduction rate
      int sizemult = X.sum() - 1;
      a[i*2 + 1] = X[i] * rate + sizemult * X[i] * interaction;
    }
  }

  void update_rates() {
    get_death_rates(); // order is important!
    get_birth_rates();
  }

  void init() {
    for (auto &rr: r) {
      rr = urd(rng);
    }
    P = log(1.0/r);
    // calculate all propensity values (a)
    // These functions should probably be optimized with sfinae or something
    update_rates();
  }


  void simulate(double interval) {
    double t_end = t + interval;
    init();

    // std::cout.precision(3);
    std::cout << "time\tsize\trate\n";
    std::cout << t;
    for (auto x: X) std::cout << '\t' << x;
    std::cout << '\t' << cells[0].get_birth_rate(t) << '\n';

    double print_interval = 0.1;
    double next_print = 0.0 + print_interval;


    while (t < t_end) {
      dt = (P - T) / a;
      // double d = dt.min();

      if (t > next_print) {
        std::cout << t;
        for (auto x: X) std::cout << '\t' << x;
        std::cout << '\t' << cells[0].get_birth_rate(t) << '\n';
        next_print += print_interval;
      }

      size_t u = std::min_element(std::begin(dt), std::end(dt)) - std::begin(dt);
      double d = dt[u];
      t += d;
      int event_celltype = u/2;
      // std::cout << u << '\n';
      if (u%2 == 0) {
        ++X[event_celltype];
      } else {
        --X[event_celltype];
      }
      T = T + a*d;
      r[u] = urd(rng);
      P[u] += log(1.0/r[u]);
      update_rates();

    }
    // also print end state
    std::cout << t;
    for (auto x: X) std::cout << '\t' << x;
    std::cout << '\t' << cells[0].get_birth_rate(t) << '\n';
    next_print += print_interval;

  }

private:
  double t = 0;
  size_t type_count;
  std::valarray<double> X;
  std::valarray<double> a;
  std::valarray<double> T;
  std::valarray<double> P;
  std::valarray<double> r;
  std::valarray<double> dt;

  std::vector<TCell> cells;

  TRng rng;
  std::uniform_real_distribution<double> urd;

};


int main(int argc, char **argv) {

  // ### Argument parsing ### //

  Arguments a;
  try {
    TCLAP::CmdLine cmd("General treatment simulator", ' ', VERSION);

    TCLAP::ValueArg<int> a_n0("n", "n0", "Starting cell count", true, 100, "integer", cmd);
    TCLAP::ValueArg<string> a_birth_rate("b", "birth-rate", "Birth rate", true, "", "[0, 1, 2]", cmd);
    TCLAP::ValueArg<double> a_death_rate("d", "death_rate", "Death rate", true, 100, "double", cmd);
    TCLAP::ValueArg<double> a_interaction_birth_rate("p", "interaction-birth-rate", "Interaction Birth rate", true, 100, "double", cmd);
    TCLAP::ValueArg<double> a_interaction_death_rate("q", "interaction-death_rate", "Interaction Death rate", true, 100, "double", cmd);
    TCLAP::ValueArg<double> a_t("t", "t-max", "Simulation time", true, 100, "double", cmd);

    cmd.parse(argc, argv);

    a.n0 = a_n0.getValue();
    string bstring = a_birth_rate.getValue();
    smatch m_b;
    regex re_b("\\d+\\.\\d+");
    while (regex_search(bstring, m_b, re_b)) {
      for (auto x: m_b) {
        a.birth_rate.push_back(stod(x));
      }
      bstring = m_b.suffix().str();
    }
    a.death_rate = a_death_rate.getValue();
    a.interaction_birth_rate = a_interaction_birth_rate.getValue();
    a.interaction_death_rate = a_interaction_death_rate.getValue();
    a.t_end = a_t.getValue();

    // Sanity checks
    assert(a.n0 > 0);
    assert(a.death_rate >= 0);
    assert(a.interaction_birth_rate >= 0);
    assert(a.interaction_death_rate >= 0);
    assert(a.t_end > 0);

  } catch (TCLAP::ArgException &e) {
    cerr << "TCLAP Error: " << e.error() << endl << "\targ: " << e.argId() << endl;
    return 1;
  }

  // ### Simulation ### //

  Cell wt(a.birth_rate, a.death_rate, a.interaction_birth_rate, a.interaction_death_rate, a.t_end);
  // LB<Cell> lb(Cell(a.birth_rate, a.death_rate, a.interaction_birth_rate, a.interaction_death_rate, a.t_end));
  LB<Cell> lb(wt);
  lb.set_cell_count(a.n0);
  lb.simulate(a.t_end);
}

