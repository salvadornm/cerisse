#ifndef CNS_RANDOM_H_
#define CNS_RANDOM_H_

#include <AMReX_IntVect.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Random.H>
#include <AMReX_SPACE.H>
#include <AMReX_Vector.H>

#include "index_macros.H"

#if NUM_FIELD > 0

/**
 * @brief A class for generating, broadcasting, and storing random numbers across
 * levels
 */
class UniqueRand
{
public:
  UniqueRand() = default;

  void init(int num_rand, int lev, amrex::Vector<amrex::IntVect> ref_ratio)
  {
    _max_lev = ref_ratio.size();
    _num_rand = num_rand;
    _lev = lev;

    // how many finest steps per coarsest step
    _num_step = 1;
    _ref_ratio.resize(_max_lev);
    for (int i = 0; i < _max_lev; ++i) {
      _ref_ratio[i] = ref_ratio[i].max();
      _num_step *= ref_ratio[i].max();
    }

    _data.resize(_num_step * _num_rand * NUM_FIELD);

    _num_per_coarse_step = 1;
    for (int i = 0; i < lev; ++i) { _num_per_coarse_step *= _ref_ratio[i]; }

    _num_sum = 1;
    for (int i = _max_lev - 1; i >= _lev; --i) { _num_sum *= _ref_ratio[i]; }

    amrex::Print() << "** Called UniqueRand.init()" << std::endl;
  }

  /**
   * @brief Generate random numbers for a coarse step
   * @param sqrtdt_coarse sqrt(dt) of the coarest level
   * @param rand_mode 0: unmodified random normal; 1: normalised randn; 2: [1,-1]; 
   *                  3: generate N/2 randn and flip sign to produce the other half (antithetic sampling)
   * @param sync_mode 0: randomise at finest level, sum up to coarser levels; 1:
   * randomise at coarsest level, divide to finer levels
   */
  void generate_new(amrex::Real sqrtdt_coarse, int rand_mode, int sync_mode)
  {
    if (_lev > 0) return; // only do this at the beginning of a coarse timestep

    if (amrex::ParallelDescriptor::IOProcessor()) {
      const amrex::Real sqrtdt =
        sqrtdt_coarse / sqrt(_num_step); // sqrtdt at this level

      if (rand_mode == 0 || rand_mode == 1) {
        // generate random normal //
        amrex::Real mean, stdd; // store for normalisation
        for (int s = 0; s < _num_step; ++s) {
          for (int d = 0; d < _num_rand; ++d) {
            mean = 0.0;
            stdd = 0.0;
            for (int nf = 0; nf < NUM_FIELD; ++nf) {
              amrex::Real rand = amrex::RandomNormal(0.0, sqrtdt);
              mean += rand;
              stdd += rand * rand;
              raw_data(s, d, nf) = rand;
            }

            // normalisation //
            if (rand_mode == 1) {
              mean /= NUM_FIELD;
              stdd = sqrt(stdd / NUM_FIELD - mean * mean);
              for (int nf = 0; nf < NUM_FIELD; ++nf) {
                raw_data(s, d, nf) = (raw_data(s, d, nf) - mean) * sqrtdt / stdd;
              }
            }
          }
        }
      } else if (rand_mode == 2) {
        // generate random [1,-1] //
        amrex::Real corr_fac =
          (NUM_FIELD % 2 == 0)
            ? 1.0
            : sqrt(amrex::Real(NUM_FIELD) /
                    amrex::Real(NUM_FIELD - 1)); // to give a variance of exactly 1
        for (int s = 0; s < _num_step; ++s) {
          for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            for (int nf = 0; nf < NUM_FIELD - 1; nf += 2) {
              raw_data(s, d, nf) = sqrtdt * corr_fac *
                                    (2.0 * amrex::Real(amrex::Random_int(2)) - 1.0);
              raw_data(s, d, nf + 1) = -raw_data(s, d, nf);
            }
            if (NUM_FIELD % 2 == 1) {
              int swap_id = amrex::Random_int(NUM_FIELD - 2);
              raw_data(s, d, NUM_FIELD - 1) = raw_data(s, d, swap_id);
              raw_data(s, d, swap_id) = 0.0;
            }
          }
        }
      } else if (rand_mode == 3) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(NUM_FIELD % 2 == 0, 
                                          "ERROR: Antithetic sampling must have even NUM_FIELD");
        // generate random normal //
        for (int s = 0; s < _num_step; ++s) {
          for (int d = 0; d < _num_rand; ++d) {
            for (int nf = 0; nf < NUM_FIELD / 2; ++nf) {
              amrex::Real rand = amrex::RandomNormal(0.0, sqrtdt);
              raw_data(s, d, nf) = rand;
              raw_data(s, d, NUM_FIELD / 2 + nf) = -rand;
            }
            // shuffle
            for (int nf = 0; nf < NUM_FIELD / 2; ++nf) {
              int swap_id = amrex::Random_int(NUM_FIELD - 1);
              amrex::Real tmp = raw_data(s, d, nf);
              raw_data(s, d, nf) = raw_data(s, d, swap_id);
              raw_data(s, d, swap_id) = tmp;
            }
          }
        }
      }

      // replace all random numbers with the first ones if sync_mode == 1 //
      if (sync_mode == 1) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          for (int nf = 0; nf < NUM_FIELD; ++nf) {
            // raw_data(0, d, nf) /=
            //   sqrt(_num_step); // this is to account for the fact that each step
            //                    // is not independent now
            // for (int s = 1; s < _num_step; ++s) {
            //   raw_data(s, d, nf) = raw_data(0, d, nf);
            // }

            // or

            amrex::Real sum = 0;
            for (int s = 0; s < _num_step; ++s) {
              sum += raw_data(s, d, nf);
            }
            sum /= _num_step;
            for (int s = 0; s < _num_step; ++s) {
              raw_data(s, d, nf) = sum;
            }
          }
        }
      }
    }

    // broadcast to all ranks //
    amrex::ParallelDescriptor::Bcast(
      _data.data(), _data.size(), amrex::ParallelDescriptor::IOProcessorNumber());

    // amrex::Print() << "** Called WienerProcess.generate_new()" << std::endl;
  }

  /**
   * @brief Get random numbers for a step
   * @param step AMRLevel::nStep()
   * @param out  output array
   */
  void get_rand(int step, amrex::Real out[][NUM_FIELD])
  {
    int cycle = step % _num_per_coarse_step;

    int start_sum = cycle * _num_sum;

    for (int d = 0; d < _num_rand; ++d) {
      for (int nf = 0; nf < NUM_FIELD; ++nf) {
        out[d][nf] = 0.0;
        for (int s = start_sum; s < start_sum + _num_sum; ++s) {
          out[d][nf] += raw_data(s, d, nf);
        }
      }
    }
  }

  /**
   * @brief Access data vector in 3D indexing
   */
  amrex::Real& raw_data(size_t i, size_t j, size_t k)
  {
    assert(i >= 0 && i < _num_step && j >= 0 && j < _num_rand && k >= 0 &&
           k < NUM_FIELD);
    return _data[i * _num_rand * NUM_FIELD + j * NUM_FIELD + k];
  }

private:
  int _lev;                      // level at which this instance is initialised
  int _max_lev;                  // max refinement level in the simulation
  amrex::Vector<int> _ref_ratio; // refinement ratio of all levels
  int _num_per_coarse_step; // number of steps in this level per coarsest level step
  int _num_sum;             // number of finest steps per step of this level
  int _num_step;            // number of finest steps per coarsest level step
  int _num_rand;            // number of random numbers requested (usually NUM_FIELD)

  static amrex::Vector<amrex::Real>
    _data; // data object, static so only one copy exists
};

#else

/**
 * @brief A dummy UniqueRand class if NUM_FIELD == 0.
 */
class UniqueRand
{
public:
  UniqueRand() = default;

  void init(int num_rand, int lev, amrex::Vector<amrex::IntVect> ref_ratio) {}

  void generate_new(amrex::Real sqrtdt_coarse, int rand_mode, int sync_mode) {}

  void get_rand(int step, amrex::Real out[][NUM_FIELD]) {}

private:
  static amrex::Vector<amrex::Real>
    _data; // data object, static so only one copy exists
};

#endif // end if NUM_FIELD > 0

#endif // end ifndef