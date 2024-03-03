/**
 * @file eos_test.cpp
 * @brief Unit tests for the EoS implementation. Compile with -DUSE_GAMMALAW_EOS
 * to compare with perfect gas EoS.
 */

#ifndef USE_GAMMALAW_EOS
static_assert(false, "USE_GAMMALAW_EOS not defined.");
#endif

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

#include <AMReX.H>
#include <PelePhysics.H>

int main(int argc, char **argv) {
  amrex::Initialize(argc, argv);
  int result = Catch::Session().run(argc, argv);
  amrex::Finalize();
  return result;
}

///////////////////////////////////////////////////////////////////////////////
// Helper functions to compare floating point numbers
// https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon

template <class T>
std::enable_if_t<not std::numeric_limits<T>::is_integer, bool>
equal_within_ulps(T x, T y, std::size_t n = 1) {
  // Since `epsilon()` is the gap size (ULP, unit in the last place)
  // of floating-point numbers in interval [1, 2), we can scale it to
  // the gap size in interval [2^e, 2^{e+1}), where `e` is the exponent
  // of `x` and `y`.

  // If `x` and `y` have different gap sizes (which means they have
  // different exponents), we take the smaller one. Taking the bigger
  // one is also reasonable, I guess.
  const T m = std::min(std::fabs(x), std::fabs(y));
  // Subnormal numbers have fixed exponent, which is `min_exponent - 1`.
  const int exp = m < std::numeric_limits<T>::min()
                      ? std::numeric_limits<T>::min_exponent - 1
                      : std::ilogb(m);
  // We consider `x` and `y` equal if the difference between them is within `n`
  // ULPs.
  return std::fabs(x - y) <=
         n * std::ldexp(std::numeric_limits<T>::epsilon(), exp);
}

template <class T>
std::enable_if_t<not std::numeric_limits<T>::is_integer, bool>
equal_within_abs(T x, T y, T abs) {  
  return std::fabs(x - y) <= abs;
}

///////////////////////////////////////////////////////////////////////////////
// Make an instance of the EoS classes

#include "Index.h"
#include "Thermodynamics.h"

calorifically_perfect_gas_t<indicies_t> perfect_gas_eos;
multispecies_gas_t<indicies_t> gamma_law_eos;

const amrex::Real R = 1.225;
const amrex::Real P = 101325;
const amrex::Real E = 206785.7142857143;
const amrex::Real Y[1] = {1.0};

///////////////////////////////////////////////////////////////////////////////
// Test cases follow

TEST_CASE("Test RYP2E()") {
  amrex::Real e1, e2;
  perfect_gas_eos.RYP2E(R, nullptr, P, e1);
  gamma_law_eos.RYP2E(R, Y, P, e2);

  REQUIRE(equal_within_ulps(e1, e2));
}

TEST_CASE("Test RYE2TP()") {
  amrex::Real t1, t2, p1, p2;
  perfect_gas_eos.RYE2TP(R, nullptr, E, t1, p1);
  gamma_law_eos.RYE2TP(R, Y, E, t2, p2);

  REQUIRE(equal_within_abs(t1, t2, 0.2)); // mw is different in PelePhysics and perfect_gas
  REQUIRE(equal_within_ulps(p1, p2));
}

TEST_CASE("Test RYE2Cs()") {
  amrex::Real cs1, cs2;
  perfect_gas_eos.RYE2Cs(R, nullptr, E, cs1);
  gamma_law_eos.RYE2Cs(R, Y, E, cs2);

  REQUIRE(equal_within_ulps(cs1, cs2));
}

TEST_CASE("Test RYE2TPCs()") {
  amrex::Real t1, t2, t3, p1, p2, p3, cs1, cs2, cs3;
  perfect_gas_eos.RYE2TP(R, nullptr, E, t1, p1);
  perfect_gas_eos.RYE2Cs(R, nullptr, E, cs1);
  gamma_law_eos.RYE2TPCs(R, Y, E, t2, p2, cs2);
  gamma_law_eos.RYE2TP(R, Y, E, t3, p3);
  gamma_law_eos.RYE2Cs(R, Y, E, cs3);

  REQUIRE(equal_within_abs(t1, t2, 0.2)); // mw is different in PelePhysics and perfect_gas
  REQUIRE(equal_within_ulps(p1, p2));
  REQUIRE(equal_within_ulps(cs1, cs2));
  REQUIRE(equal_within_ulps(t2, t3));
  REQUIRE(equal_within_ulps(p2, p3));
  REQUIRE(equal_within_ulps(cs2, cs3));
}