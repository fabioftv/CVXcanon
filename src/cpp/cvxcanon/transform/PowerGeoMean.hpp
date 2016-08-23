#ifndef CVXCANON_TRANSFORM_POWER_GEO_MEAN_H
#define CVXCANON_TRANSFORM_POWER_GEO_MEAN_H

#include "cvxcanon/transform/ProblemTransform.hpp"

class PowerGeoMean : public ProblemTransform {
   public:
      PowerGeoMean();
      ~PowerGeoMean();

      bool accepts(const Problem& problem) override;
      Problem apply(const Problem& problem) override;
};

class GeoMeanIneq {
   public:
      GeoMeanIneq();
      ~GeoMeanIneq();

      long gcd(long a, long b);
      std::pair<int, int> fraction(double number);
      std::pair<std::pair<int, int>, std::pair<std::pair<int, int>, 
         std::pair<int, int>>> build_power_pos(double number);
      std::pair<std::pair<int, int>, std::pair<std::pair<int, int>, 
         std::pair<int, int>>> build_power_neg(double number);
      std::pair<std::pair<int, int>, std::pair<std::pair<int, int>, 
         std::pair<int, int>>> build_power_middle(double number);
      bool is_integer(double number);
      bool power_two_test(double number);
      bool dyadic_nonnegative_fraction_test(std::pair<double, double> fraction);
      bool weight_vector_test(std::vector<std::pair<double, double>> w);
      bool dyadic_weight_vector_test(std::vector<std::pair<double, double>> w);

      std::vector<int> sort(std::vector<double> num);
//TODO(fatbioftv): Implement Merge Sort
//    void merge(std::vector<double> num, int low, int middle, int high);
//    void merge_sort(std::vector<double> num, int low, int high);
      std::vector<std::pair<int, int>> make_frac(std::vector<double> a, 
                                                 int denominator);
      std::vector<std::pair<int, int>> 
         dyad_completion(std::vector<std::pair<int, int>> w);
      double approx_error(std::vector<std::pair<double, double>> a_orig, 
                          std::vector<std::pair<double, double>> w_approx);
      int next_power_two(int number);
      bool check_dyad(std::vector<std::pair<double, double>> w, 
                      std::vector<std::pair<double, double>> w_dyad);
      std::pair<std::vector<std::pair<int, int>>, std::vector<std::pair<int,
         int>>> split(std::vector<std::pair<double, double>> wdyad);
};

#endif  // CVXCANON_TRANSFORM_POWER_GEO_MEAN_H
