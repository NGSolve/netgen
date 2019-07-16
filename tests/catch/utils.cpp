
#include "catch.hpp"
#include <core/ngcore.hpp>
using namespace ngcore;
using namespace std;


long shuffle(long N, long i) {
    // Shuffle the numbers using multiplication with a prime number to force many updates of min, max
    constexpr long P = 101;
    return (N/2 + i*P) % N;
}

void testThreading(int n_threads)
{
  TaskManager::SetNumThreads(n_threads);
  n_threads = EnterTaskManager();

  constexpr long N = 100000;


  SECTION( "atomic operations" ) {
      long i_min = 2*N;
      long i_max = 0;
      long i_sum = 0;

      double d_min = 1e100;
      double d_max = 0.0;
      double d_sum = 0.0;

      ParallelFor( Range(N), [&] (long i) {
          AtomicMin(i_min, shuffle(N,i));
      });
      REQUIRE( i_min==0 );

      ParallelFor( Range(N), [&] (long i) {
          AtomicMax(i_max, shuffle(N,i));
      });
      REQUIRE( i_max==N-1 );

      ParallelFor( Range(N), [&] (long i) {
          AsAtomic(i_sum) += i;
      });
      REQUIRE( i_sum==N*(N-1)/2 );

      ParallelFor( Range(N), [&] (double i) {
          AtomicMin(d_min, static_cast<double>(shuffle(N,i)));
      });
      REQUIRE( d_min==0 );

      ParallelFor( Range(N), [&] (double i) {
          AtomicMax(d_max, static_cast<double>(shuffle(N,i)));
      });
      REQUIRE( d_max==N-1 );

      ParallelFor( Range(N), [&] (double i) {
          AtomicAdd(d_sum, i);
      });
      REQUIRE( d_sum==N*(N-1)/2 );

  }
  ExitTaskManager(n_threads);
}

TEST_CASE("Threading - 1 Thread") { testThreading(1); }
TEST_CASE("Threading - 2 Thread") { testThreading(2); }
TEST_CASE("Threading - 8 Thread") { testThreading(8); }
