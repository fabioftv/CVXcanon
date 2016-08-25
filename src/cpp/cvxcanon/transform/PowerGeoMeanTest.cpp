#include "cvxcanon/transform/PowerGeoMean.hpp"

#include <cmath>
#include <unordered_map>
#include <vector>

#include "gtest/gtest.h"

TEST(PowerGeoMeanTest, fraction) {
   GeoMeanIneq geo_mean;

   std::pair<int, int> frac;
   frac = geo_mean.fraction(0.5);
   EXPECT_EQ(1, frac.first);
   EXPECT_EQ(2, frac.second);

   frac = geo_mean.fraction(0.75);
   EXPECT_EQ(3, frac.first);
   EXPECT_EQ(4, frac.second);

   frac = geo_mean.fraction(0.54);
   EXPECT_EQ(27, frac.first);
   EXPECT_EQ(50, frac.second);

   frac = geo_mean.fraction(2);
   EXPECT_EQ(2, frac.first);
   EXPECT_EQ(1, frac.second);
}

TEST(PowerGeoMeanTest, build_power_pos) {
   GeoMeanIneq geo_mean;

   std::pair<std::pair<int, int>, std::pair<std::pair<int, int>, 
      std::pair<int, int>>> pow_high;

   pow_high = geo_mean.build_power_pos(2.75);

   EXPECT_EQ(11, pow_high.first.first);
   EXPECT_EQ(4, pow_high.first.second);
   EXPECT_EQ(4, pow_high.second.first.first);
   EXPECT_EQ(11, pow_high.second.first.second);
   EXPECT_EQ(7, pow_high.second.second.first);
   EXPECT_EQ(11, pow_high.second.second.second);

   pow_high = geo_mean.build_power_pos(3.3);

   EXPECT_EQ(33, pow_high.first.first);
   EXPECT_EQ(10, pow_high.first.second);
   EXPECT_EQ(10, pow_high.second.first.first);
   EXPECT_EQ(33, pow_high.second.first.second);
   EXPECT_EQ(23, pow_high.second.second.first);
   EXPECT_EQ(33, pow_high.second.second.second);
}

TEST(PowerGeoMeanTest, build_power_neg) {
   GeoMeanIneq geo_mean;

   std::pair<std::pair<int, int>, std::pair<std::pair<int, int>, 
      std::pair<int, int>>> pow_low;

   pow_low = geo_mean.build_power_neg(-3.54);

   EXPECT_EQ(-177, pow_low.first.first);
   EXPECT_EQ(50, pow_low.first.second);
   EXPECT_EQ(177, pow_low.second.first.first);
   EXPECT_EQ(227, pow_low.second.first.second);
   EXPECT_EQ(50, pow_low.second.second.first);
   EXPECT_EQ(227, pow_low.second.second.second);
}

TEST(PowerGeoMeanTest, build_power_middle) {
   GeoMeanIneq geo_mean;

   std::pair<std::pair<int, int>, std::pair<std::pair<int, int>, 
      std::pair<int, int>>> pow_middle;

   pow_middle = geo_mean.build_power_middle(0.75);

   EXPECT_EQ(3, pow_middle.first.first);
   EXPECT_EQ(4, pow_middle.first.second);
   EXPECT_EQ(3, pow_middle.second.first.first);
   EXPECT_EQ(4, pow_middle.second.first.second);
   EXPECT_EQ(1, pow_middle.second.second.first);
   EXPECT_EQ(4, pow_middle.second.second.second);
}

TEST(PowerGeoMeanTest, is_integer) {
   GeoMeanIneq geo_mean;

   EXPECT_EQ(true, geo_mean.is_integer(4));
   EXPECT_EQ(false, geo_mean.is_integer(4.25));
   EXPECT_EQ(false, geo_mean.is_integer(4.00000001));
   EXPECT_EQ(true, geo_mean.is_integer(4.000000001));
}

TEST(PowerGeoMeanTest, power_two_test) {
   GeoMeanIneq geo_mean;

   EXPECT_EQ(true, geo_mean.power_two_test(4));
   EXPECT_EQ(true, geo_mean.power_two_test(2^10));
   EXPECT_EQ(true, geo_mean.power_two_test(1));
   EXPECT_EQ(true, geo_mean.power_two_test(1.0));
   EXPECT_EQ(false, geo_mean.power_two_test(0));
   EXPECT_EQ(false, geo_mean.power_two_test(-4));
}

TEST(PowerGeoMeanTest, Dyadic_Nonnegative_Fraction_Test) {
   GeoMeanIneq geo_mean;
   
   std::pair<double, double> frac;
   frac.first = 1;
   frac.second = 4;
   EXPECT_EQ(true, geo_mean.dyadic_nonnegative_fraction_test(frac));

   frac.second = 3;
   EXPECT_EQ(false, geo_mean.dyadic_nonnegative_fraction_test(frac));

   frac.first = 0;
   frac.second = 1;
   EXPECT_EQ(true, geo_mean.dyadic_nonnegative_fraction_test(frac));

   frac.first = 1;
   frac.second = 1;
   EXPECT_EQ(true, geo_mean.dyadic_nonnegative_fraction_test(frac));

   frac.first = -1;
   frac.second = 4;
   EXPECT_EQ(false, geo_mean.dyadic_nonnegative_fraction_test(frac));

   frac.first = 1;
   frac.second = 6;
   EXPECT_EQ(false, geo_mean.dyadic_nonnegative_fraction_test(frac));

   frac.first = 1.2;
   frac.second = 4;
   EXPECT_EQ(false, geo_mean.dyadic_nonnegative_fraction_test(frac));
}


TEST(PowerGeoMeanTest, weight_vector_test) {
   GeoMeanIneq geo_mean;

   std::vector<std::pair<double, double>> w(2);
   w[0].first = 1; w[0].second = 3; w[1].first = 2; w[1].second = 3;
   EXPECT_EQ(true, geo_mean.weight_vector_test(w));

   w[0].first = 2;
   EXPECT_EQ(false, geo_mean.weight_vector_test(w));

   w[0].first = 0.1; w[0].second = 1; w[1].first = 0.9; w[1].second = 1;
   EXPECT_EQ(false, geo_mean.weight_vector_test(w));

   std::vector<std::pair<double, double>> v(3);
   v[0].first = 0; v[0].second = 1; v[1].first = 0; v[1].second = 1;
   v[2].first = 1; v[2].second = 1;
   EXPECT_EQ(true, geo_mean.weight_vector_test(v));
}


TEST(PowerGeoMeanTest, dyadic_weight_vector_test) {
   GeoMeanIneq geo_mean;

   std::vector<std::pair<double, double>> w(2);
   w[0].first = 1; w[0].second = 2; w[1].first = 1; w[1].second = 2;
   EXPECT_EQ(true, geo_mean.dyadic_weight_vector_test(w));

   w[0].second = 3; w[1].first = 2; w[1].second = 3;
   EXPECT_EQ(false, geo_mean.dyadic_weight_vector_test(w));

   std::vector<std::pair<double, double>> v(3);
   v[0].first = 0; v[0].second = 1; v[1].first = 1; v[1].second = 1;
   v[2].first = 0; v[2].second = 1;
   EXPECT_EQ(true, geo_mean.dyadic_weight_vector_test(v));
}

TEST(PowerGeoMean, sort) {
   GeoMeanIneq geo_mean;

   std::vector<double> num(6);
   num[0] = 12.5;
   num[1] = 11.2;
   num[2] = 13.9;
   num[3] = 5.6;
   num[4] = 6.3;
   num[5] = 7.5;

   std::vector<int> index(6);
   index = geo_mean.sort(num);

   EXPECT_EQ(3, index[0]);
   EXPECT_EQ(4, index[1]);
   EXPECT_EQ(5, index[2]);
   EXPECT_EQ(1, index[3]);
   EXPECT_EQ(0, index[4]);
   EXPECT_EQ(2, index[5]);
}

TEST(PowerGeoMean, fracify) {
   GeoMeanIneq geo_mean;
}

TEST(PowerGeoMean, make_frac) {
   GeoMeanIneq geo_mean;
   std::vector<std::pair<double, double>> frac;

   std::vector<double> a(3);
   a[0] = 0.123; a[1] = 0.345; a[2] = 0.532;

   int denominator = 10;
   frac = geo_mean.make_frac(a, denominator);
   EXPECT_EQ(1, frac[0].first);
   EXPECT_EQ(10, frac[0].second);
   EXPECT_EQ(4, frac[1].first);
   EXPECT_EQ(10, frac[1].second);
   EXPECT_EQ(5, frac[2].first);
   EXPECT_EQ(10, frac[2].second);

   denominator = 100;
   frac = geo_mean.make_frac(a, denominator);
   EXPECT_EQ(12, frac[0].first);
   EXPECT_EQ(100, frac[0].second);
   EXPECT_EQ(35, frac[1].first);
   EXPECT_EQ(100, frac[1].second);
   EXPECT_EQ(53, frac[2].first);
   EXPECT_EQ(100, frac[2].second);

   denominator = 1000;
   frac = geo_mean.make_frac(a, denominator);
   EXPECT_EQ(123, frac[0].first);
   EXPECT_EQ(1000, frac[0].second);
   EXPECT_EQ(345, frac[1].first);
   EXPECT_EQ(1000, frac[1].second);
   EXPECT_EQ(532, frac[2].first);
   EXPECT_EQ(1000, frac[2].second);
}

TEST(PowerGeoMean, dyad_completion){
   GeoMeanIneq geo_mean;

   std::vector<std::pair<double, double>> w(3);
   w[0].first = 1; w[0].second = 3; w[1].first = 1; w[1].second = 5;
   w[2].first = 7; w[2].second = 15;

   std::vector<std::pair<double, double>> dyad_completion;
   dyad_completion = geo_mean.dyad_completion(w);
   
   EXPECT_EQ(5, dyad_completion[0].first);
   EXPECT_EQ(16, dyad_completion[0].second);
   EXPECT_EQ(3, dyad_completion[1].first);
   EXPECT_EQ(16, dyad_completion[1].second);
   EXPECT_EQ(7, dyad_completion[2].first);
   EXPECT_EQ(16, dyad_completion[2].second);
   EXPECT_EQ(1, dyad_completion[3].first);
   EXPECT_EQ(16, dyad_completion[3].second);

   w[0].first = 1; w[0].second = 3; w[1].first = 1; w[1].second = 3;
   w[2].first = 1; w[2].second = 3;

   dyad_completion = geo_mean.dyad_completion(w);

   EXPECT_EQ(1, dyad_completion[0].first);
   EXPECT_EQ(4, dyad_completion[0].second);
   EXPECT_EQ(1, dyad_completion[1].first);
   EXPECT_EQ(4, dyad_completion[1].second);
   EXPECT_EQ(1, dyad_completion[2].first);
   EXPECT_EQ(4, dyad_completion[2].second);
   EXPECT_EQ(1, dyad_completion[3].first);
   EXPECT_EQ(4, dyad_completion[3].second);

   std::vector<std::pair<double, double>> v(4);
   v[0].first = 1; v[0].second = 1; v[1].first = 0; v[1].second = 1;
   v[2].first = 0; v[2].second = 1; v[3].first = 0; v[3].second = 1;

   dyad_completion = geo_mean.dyad_completion(v);

   EXPECT_EQ(1, dyad_completion[0].first);
   EXPECT_EQ(1, dyad_completion[0].second);
   EXPECT_EQ(0, dyad_completion[1].first);
   EXPECT_EQ(1, dyad_completion[1].second);
   EXPECT_EQ(0, dyad_completion[2].first);
   EXPECT_EQ(1, dyad_completion[2].second);
   EXPECT_EQ(0, dyad_completion[3].first);
   EXPECT_EQ(1, dyad_completion[3].second);
}

TEST(PowerGeoMean, approx_error) {
   GeoMeanIneq geo_mean;

   std::vector<std::pair<double, double>> w(2);
   w[0].first = 1; w[0].second = 3; w[1].first = 2; w[1].second = 3;
   std::vector<std::pair<double, double>> v(2);
   v[0].first = 1; v[0].second = 4; v[1].first = 3; v[1].second = 4;

   double max = geo_mean.approx_error(w, v);
   EXPECT_LE((w[0].first / w[0].second) - (v[0].first / v[0].second), max);
}

TEST(PowerGeoMean, next_power_two) {
   GeoMeanIneq geo_mean;

   EXPECT_EQ(4, geo_mean.next_power_two(3));
   EXPECT_EQ(8, geo_mean.next_power_two(8));
   EXPECT_EQ(1, geo_mean.next_power_two(0));
   EXPECT_EQ(1, geo_mean.next_power_two(1));
}

TEST(PowerGeoMean, check_dyad) {
   GeoMeanIneq geo_mean;

   std::vector<std::pair<double, double>> w(3);
   w[0].first = 1; w[0].second = 3; w[1].first = 1; w[1].second = 3;
   w[2].first = 1; w[2].second = 3;

   std::vector<std::pair<double, double>> w_dyad(4);
   w_dyad[0].first = 1; w_dyad[0].second = 4; w_dyad[1].first = 1; 
   w_dyad[1].second = 4; w_dyad[2].first = 1; w_dyad[2].second = 4;
   w_dyad[3].first = 1; w_dyad[3].second = 4;
   EXPECT_EQ(true, geo_mean.check_dyad(w, w_dyad));

   w[0].first = 1; w[0].second = 4; w[1].first = 0; w[1].second = 1;
   w[2].first = 3; w[2].second = 4;
   EXPECT_EQ(true, geo_mean.check_dyad(w, w));

   w[0].first = 1; w[0].second = 1; w[1].first = 0; w[1].second = 1;
   w[2].first = 0; w[2].second = 1;
   EXPECT_EQ(true, geo_mean.check_dyad(w, w));

   std::vector<std::pair<double, double>> v(2);
   v[0].first = 2; v[0].second = 3; v[1].first = 1; v[1].second = 1;
   EXPECT_EQ(false, geo_mean.check_dyad(v, v));

   v[0].first = 2; v[0].second = 5; v[1].first = 3; v[1].second = 5;

   std::vector<std::pair<double, double>> v_dyad(3);
   v_dyad[0].first = 3; v_dyad[0].second = 8; v_dyad[1].first = 4; 
   v_dyad[1].second = 8; v_dyad[2].first = 1; v_dyad[2].second = 8;
   EXPECT_EQ(false, geo_mean.check_dyad(v, v_dyad));
}

TEST(PowerGeoMean, split) {
   GeoMeanIneq geo_mean;

   std::vector<std::pair<double, double>> w_dyad(3);
   w_dyad[0].first = 3; w_dyad[0].second = 8; w_dyad[1].first = 4; 
   w_dyad[1].second = 8; w_dyad[2].first = 1; w_dyad[2].second = 8;

   std::vector<std::pair<double, double>> child_1(w_dyad.size());
   std::vector<std::pair<double, double>> child_2(w_dyad.size());
   std::pair<std::vector<std::pair<double, double>>, 
      std::vector<std::pair<double, double>>> split_w_dyad;

   split_w_dyad = geo_mean.split(w_dyad);
   child_1 = split_w_dyad.first;
   child_2 = split_w_dyad.second;

   EXPECT_EQ(0, child_1[0].first);
   EXPECT_EQ(1, child_1[0].second);
   EXPECT_EQ(1, child_1[1].first);
   EXPECT_EQ(1, child_1[1].second);
   EXPECT_EQ(0, child_1[2].first);
   EXPECT_EQ(1, child_1[2].second);
   EXPECT_EQ(6, child_2[0].first);
   EXPECT_EQ(8, child_2[0].second);
   EXPECT_EQ(0, child_2[1].first);
   EXPECT_EQ(1, child_2[1].second);
   EXPECT_EQ(2, child_2[2].first);
   EXPECT_EQ(8, child_2[2].second);

   w_dyad[0].first = 0; w_dyad[0].second = 1; w_dyad[1].first = 1; 
   w_dyad[1].second = 1; w_dyad[2].first = 0; w_dyad[2].second = 1;

   split_w_dyad = geo_mean.split(w_dyad);
   child_1 = split_w_dyad.first;
   child_2 = split_w_dyad.second;

   EXPECT_EQ(0, child_1[0].first);
   EXPECT_EQ(0, child_1[0].second);
   EXPECT_EQ(0, child_1[1].first);
   EXPECT_EQ(0, child_1[1].second);
   EXPECT_EQ(0, child_1[2].first);
   EXPECT_EQ(0, child_1[2].second);
   EXPECT_EQ(0, child_2[0].first);
   EXPECT_EQ(0, child_2[0].second);
   EXPECT_EQ(0, child_2[1].first);
   EXPECT_EQ(0, child_2[1].second);
   EXPECT_EQ(0, child_2[2].first);
   EXPECT_EQ(0, child_2[2].second);
}

TEST(PowerGeoMean, get_max_denom) {
   GeoMeanIneq geo_mean;

   std::vector<std::pair<double, double>> tup(3);
   tup[0].first = 3; tup[0].second = 8; tup[1].first = 4; 
   tup[1].second = 8; tup[2].first = 1; tup[2].second = 8;
   EXPECT_EQ(8, geo_mean.get_max_denom(tup));

   tup[0].first = 2; tup[0].second = 3; tup[1].first = 1; 
   tup[1].second = 3; tup[2].first = 1; tup[2].second = 5;
   EXPECT_EQ(5, geo_mean.get_max_denom(tup));
}

TEST(PowerGeoMean, get_number_of_digits) {
   GeoMeanIneq geo_mean;

   EXPECT_EQ(6, geo_mean.get_number_of_digits(153124));
   EXPECT_EQ(1, geo_mean.get_number_of_digits(1));
   EXPECT_EQ(1, geo_mean.get_number_of_digits(0));
}

TEST(PowerGeoMean, int_to_binary) {
   GeoMeanIneq geo_mean;

   EXPECT_EQ(1000, geo_mean.int_to_binary(8));
   EXPECT_EQ(11000, geo_mean.int_to_binary(24));
   EXPECT_EQ(0, geo_mean.int_to_binary(0));
   EXPECT_EQ(1, geo_mean.int_to_binary(1));
}

TEST(PowerGeoMean, lower_bound) {
   GeoMeanIneq geo_mean;

   std::vector<std::pair<double, double>> w_dyad(3);
   w_dyad[0].first = 0; w_dyad[0].second = 1; w_dyad[1].first = 1; 
   w_dyad[1].second = 1; w_dyad[2].first = 0; w_dyad[2].second = 1;
   EXPECT_EQ(0, geo_mean.lower_bound(w_dyad));

   std::vector<std::pair<double, double>> v_dyad(2);
   v_dyad[0].first = 1; v_dyad[0].second = 2; v_dyad[1].first = 1; 
   v_dyad[1].second = 2;
   EXPECT_EQ(1, geo_mean.lower_bound(v_dyad));

   v_dyad[0].first = 1; v_dyad[0].second = 8; v_dyad[1].first = 7; 
   v_dyad[1].second = 8;
   EXPECT_EQ(3, geo_mean.lower_bound(v_dyad));

   std::vector<std::pair<double, double>> z_dyad(4);
   z_dyad[0].first = 1; z_dyad[0].second = 4; z_dyad[1].first = 1; 
   z_dyad[1].second = 4; z_dyad[2].first = 1; z_dyad[2].second = 4;
   z_dyad[3].first = 1; z_dyad[3].second = 4;
   EXPECT_EQ(3, geo_mean.lower_bound(z_dyad));
}
