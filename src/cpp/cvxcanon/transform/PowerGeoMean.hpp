#ifndef CVXCANON_TRANSFORM_POWER_GEO_MEAN_H
#define CVXCANON_TRANSFORM_POWER_GEO_MEAN_H

#include <vector>

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
      std::vector<int> sort(std::vector<double> test);
      std::pair<std::vector<std::pair<double, double>>, 
         std::vector<std::pair<double, double>>> fracify
            (std::vector<std::pair<double, double>> a, int max_denom, 
               bool force_dyad);
      std::vector<std::pair<double, double>> make_frac(std::vector<double> a, 
                                                 int denominator);
      std::vector<std::pair<double, double>> 
         dyad_completion(std::vector<std::pair<double, double>> w);
      double approx_error(std::vector<std::pair<double, double>> a_orig, 
                          std::vector<std::pair<double, double>> w_approx);
      int next_power_two(int number);
      bool check_dyad(std::vector<std::pair<double, double>> w, 
                      std::vector<std::pair<double, double>> w_dyad);
      std::pair<std::vector<std::pair<double, double>>, 
         std::vector<std::pair<double, double>>> split 
            (std::vector<std::pair<double, double>> wdyad);
      std::vector<std::pair<std::vector<std::pair<double, double>>, 
         std::pair<std::vector<std::pair<double, double>>, 
            std::vector<std::pair<double, double>>>>> decompose
               (std::vector<std::pair<double, double>> w_dyad);
      double get_max_denom(std::vector<std::pair<double, double>> tup);
      unsigned get_number_of_digits (unsigned number);
      int int_to_binary(int number);
      double lower_bound(std::vector<std::pair<double, double>> w_dyad);
};

#endif  // CVXCANON_TRANSFORM_POWER_GEO_MEAN_H
