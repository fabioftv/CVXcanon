#include "cvxcanon/transform/PowerGeoMean.hpp"

#include <unordered_map>
#include <vector>
#include <math.h>
#include <iostream>
#include <assert.h>
#include <algorithm>

#include "cvxcanon/expression/Expression.hpp"
#include "cvxcanon/expression/ExpressionShape.hpp"
#include "cvxcanon/expression/ExpressionUtil.hpp"
#include "cvxcanon/expression/TextFormat.hpp"
#include "cvxcanon/util/MatrixUtil.hpp"
#include "glog/logging.h"

#define TOLERANCE 0.00000001
#define PRECISION 1000000000

GeoMeanIneq::GeoMeanIneq() {}
GeoMeanIneq::~GeoMeanIneq() {}

/*
typedef Expression(*TransformFunction)(
    const Expression& expr,
    std::vector<Expression>* constraints);

// t <= x0^p0 * x1^p1 * ... * xn^pn

Expression form_geo_mean_ineq(
           std::vector<double>* t, 
           std::vector<double>* variables, 
           std::vector<double>* p) {

// TODO(fabioftv): Build these functions
   weight_test(p);
   w = dyad_completion(p);
   tree = decompose(w);

   std::vector<double> d;
   d(w) = t;

   if (variables.size < w.size){
      variables += t;
   }
   variables.size = w.size;
   
   std::vector<int> tmp;

   for (int i = 0; i < p.size; i++){
      if (p[i] > 0){
         if (i = 0){
            tmp[i] = size.w;
         }
         else{
            tmp[i] = 1;
         }
      }
      d[tpm[i]] = v;
// TODO(fabioftv): Build v
   }

   std::vector<Expression>* constraints;

   return constraints

}

*/

long GeoMeanIneq::gcd(long a, long b) {
   if (a == 0){
      return b;
   }
   else if (b == 0){
      return a;
   }
   if (a < b){
      return gcd(a, b % a);
   }
   else{
      return gcd(b, a % b);
   }
}

std::pair<int, int> GeoMeanIneq::fraction(double number) {
   long gcd_aux;
   long denom;
   long num;
   int i;

   if (number >= 0){
      gcd_aux = gcd(round(number * PRECISION), PRECISION);
      denom = PRECISION / gcd_aux;
      num = round(number * PRECISION) / gcd_aux;
   }
   if (number < 0){
      gcd_aux = gcd(round(-number * PRECISION), PRECISION);
      denom = PRECISION / gcd_aux;
      num = -round(-number * PRECISION) / gcd_aux;
   }

   return std::make_pair(num, denom);
}

std::pair<std::pair<int, int>, std::pair<std::pair<int, int>, 
   std::pair<int, int>>> GeoMeanIneq::build_power_pos(double number) {
   
   assert(number > 1);   
   
   std::pair<int, int> frac;
   std::pair<int, int> frac_inv;
   std::pair<int, int> frac_minus;
   
   frac = fraction(number);
   frac_inv.first = frac.second;
   frac_inv.second = frac.first;
   frac_minus.first = frac_inv.second - frac_inv.first;
   frac_minus.second = frac_inv.second;

   return std::make_pair(frac, std::make_pair(frac_inv, frac_minus));
}

std::pair<std::pair<int, int>, std::pair<std::pair<int, int>, 
   std::pair<int, int>>> GeoMeanIneq::build_power_neg(double number) {
   
   assert(number < 0);
   
   std::pair<int, int> frac;
   std::pair<int, int> frac_div;
   std::pair<int, int> frac_minus;
   
   frac = fraction(number);
   frac_div.first = -frac.first;
   frac_div.second = (-frac.first) + frac.second;
   frac_minus.first = frac_div.second - frac_div.first;
   frac_minus.second = frac_div.second;

   return std::make_pair(frac, std::make_pair(frac_div, frac_minus));
}

std::pair<std::pair<int, int>, std::pair<std::pair<int, int>, 
   std::pair<int, int>>> GeoMeanIneq::build_power_middle(double number) {
   
   assert(number > 0 && number < 1);
   
   std::pair<int, int> frac;
   std::pair<int, int> frac_minus;
   
   frac = fraction(number);
   frac_minus.first = frac.second - frac.first;
   frac_minus.second = frac.second;

   return std::make_pair(frac, std::make_pair(frac, frac_minus));
}

bool GeoMeanIneq::is_integer(double number) {
   if (number > (floor(number)-TOLERANCE) && number < (floor(number)+TOLERANCE)){
      return true;
   }
   else{
      return false;
   }
}

bool GeoMeanIneq::power_two_test(double number) {
   if (is_integer(number) == true){
      if ((int) number > 0 && (! ((int) number & ((int) number -1))) == true){
         return true;
      }
      else{
         return false;
      }
   }
   else{
      return false;
   }  
}

bool GeoMeanIneq::dyadic_nonnegative_fraction_test(std::pair<double, double>
     fraction) {
   if ((fraction.first / fraction.second) == fraction.first && 
       is_integer(fraction.first) == true &&
       fraction.first >= 0){
      return true;
   }
   else if ((fraction.first / fraction.second) != fraction.first && 
            is_integer(fraction.first) == true && 
            is_integer(fraction.second) == true && 
            fraction.first / fraction.second  >= 0 && 
            power_two_test(fraction.second) == true){
      return true;
   }
   else {
      return false;
   }
}

bool GeoMeanIneq::weight_vector_test(std::vector<std::pair<double, double>> w){
   double sum = 0;   
   bool aux = true;
   int i = 0;

   while (i < w.size()){
      if ((w[i].first / w[i].second) == w[i].first){
         if (is_integer(w[i].first) == true && w[i].first >= 0){
            sum += w[i].first;
         }
         else{
            aux = false;
            i = w.size(); 
         }
      }
      else if ((w[i].first / w[i].second) != w[i].first){
         if (is_integer(w[i].first) == true && is_integer(w[i].second) == true
             && w[i].first / w[i].second >= 0){
            sum += w[i].first / w[i].second;
         }
         else{
            aux = false;
            i = w.size();          
         }
      }
      i++;
      if (sum > (1 + TOLERANCE) && i <= w.size()){
         aux = false;
         i = w.size();
      }
   }
   return aux;
}

bool GeoMeanIneq::dyadic_weight_vector_test(std::vector<std::pair<double, 
                                            double>> w){
   bool aux = true;
   int i = 0;

   while (i < w.size()){
      if (dyadic_nonnegative_fraction_test(w[i]) != true){
         aux = false;
         i = w.size();  
      }
      i++;
   }
   if (aux == true && weight_vector_test(w) == true){
      return true;
   }
   else{
      return false;
   }
}

std::vector<int> GeoMeanIneq::sort(std::vector<double> test){

   std::vector<int> y(test.size());
   std::size_t n(0);
   std::generate(std::begin(y), std::end(y), [&]{return n++;});
   std::sort (std::begin(y), std::end(y), [&](int i1, int i2) 
      {return test[i1] < test[i2];});

   return y;
}

std::vector<std::pair<int, int>> GeoMeanIneq::make_frac(std::vector<double>
                                 a, int denominator) {
   double sum_a = 0;
   int i = 0, sum_b = 0;
   std::vector<int> b(a.size());
   std::vector<double> error(a.size());
   std::vector<int> index_sort(a.size());

   for (i = 0; i < a.size(); i++){
      sum_a += a[i];
   }
   for (i = 0; i < a.size(); i++){
      a[i] = a[i] / sum_a;
   }
   for (i = 0; i < b.size(); i++){
      b[i] = (int) (a[i] * denominator);
      sum_b += b[i];
   }
   for (i = 0; i < error.size(); i++){
      error[i] = ((double) b[i] /  denominator) - a[i];
   }

   index_sort = sort(error);
   std::vector<int> inds(denominator - sum_b);
   for (i = 0; i < (denominator - sum_b); i++){
      inds[i] = index_sort[i];
      b[inds[i]] = b[inds[i]] + 1;
   }
   std::vector<std::pair<int, int>> frac(a.size());
   for (i = 0; i < frac.size(); i++){
      frac[i] = std::make_pair(b[i], denominator);
   }
   
   return frac;
}

std::vector<std::pair<int, int>> 
   GeoMeanIneq::dyad_completion(std::vector<std::pair<int, int>> w) {
   int i, d = 0, p;
   double div = 0;

   for (i = 0; i < w.size(); i++){
      if (w[i].second > d){
         d = w[i].second;
      }
   }

   p = next_power_two(d);

   if (p == d){
      return w;
   }
   else {
      std::vector<std::pair<int, int>> temp(w.size() + 1);
      for (i = 0; i < w.size(); i++){
         div = (((double) w[i].first) * ((double) d)) / 
               (((double) w[i].second) * ((double) p));
         std::pair<int, int> frac;
         frac = fraction(div);
         temp[i].first = frac.first;
         temp[i].second = frac.second;
      }
      temp[w.size()].first = (double) (p - d);
      temp[w.size()].second = p;
      return temp;
   }
}

double GeoMeanIneq::approx_error(std::vector<std::pair<double, double>> a_orig, 
                                 std::vector<std::pair<double, double>> 
                                 w_approx) {
   bool pos_vec = true;
   double sum_a_orig = 0, max = 0;
   int i = 0;

   for (i = 0; i < a_orig.size(); i++){
      if (a_orig[i].first / a_orig[i].second < 0){
         pos_vec = false;
      }
      else{
         sum_a_orig += a_orig[i].first / a_orig[i].second;
      }
   }
   
   assert (pos_vec == true);
   assert (weight_vector_test(w_approx) == true);
   assert (a_orig.size() == w_approx.size());

   std::vector<double> w_orig(a_orig.size());
   for (i = 0; i < w_orig.size(); i++){
      w_orig[i] = (a_orig[i].first / a_orig[i].second) / sum_a_orig;
   }

   for (i = 0; i < w_orig.size(); i++){
      if (fabs(w_orig[i] - (w_approx[i].first / w_approx[i].second)) > max){
         max = fabs(w_orig[i] - (w_approx[i].first / w_approx[i].second));
      }
   }

   return max;
}

int GeoMeanIneq::next_power_two(int number){
   if (number <= 0){
      return 1;
   }
   else{
      int next = pow(2, ceil(log(number) / log(2)));
      return next;
   }
}

bool GeoMeanIneq::check_dyad(std::vector<std::pair<double, double>> w, 
                             std::vector<std::pair<double, double>> w_dyad) {
   bool eq_vec = true;
   bool vec_one = true;
   double div = 0;
   int i;

   if (weight_vector_test(w) == false || 
       dyadic_weight_vector_test(w_dyad) == false){
      return false;
   }
   else if (w == w_dyad){
      return true;
   }
   else if(w_dyad.size() == (w.size() + 1)){
      std::vector<std::pair<double, double>> temp(w_dyad.size() - 1);
      std::vector<double> w_double(w.size());
      for (i = 0; i < temp.size(); i++){
         temp[i].first = w_dyad[i].first;
         temp[i].second = w_dyad[i].second;
         w_double[i] = w[i].first / w[i].second;
      }
      std::vector<std::pair<double, double>> new_dyad(temp.size());
      std::vector<double> new_dyad_double(new_dyad.size());
      for (i = 0; i < new_dyad.size(); i++){
         new_dyad[i].first = temp[i].first * w_dyad[w_dyad.size() - 1].second;
         new_dyad[i].second = temp[i].second * (w_dyad[w_dyad.size() - 1].second
                              - w_dyad[w_dyad.size() - 1].first);
        new_dyad_double[i] = new_dyad[i].first / new_dyad[i].second;
      }
      if (w_double == new_dyad_double){
         return true;
      }
      else{
         return false;
      }
   }
   else{
      return false;
   }
}

std::pair<std::vector<std::pair<double, double>>, std::vector<std::pair<double,
   double>>> GeoMeanIneq::split(std::vector<std::pair<double, double>> w_dyad) {
   bool condition = true, stat_a = true, stat_b = true;
   double sum = 0;
   int i, j;

   for (i = 0; i < w_dyad.size(); i++){
      if (w_dyad[i].first == w_dyad[i].second){
         condition = false;
         i = w_dyad.size();
      }
   }
//TODO(fabioftv): Need to return empty
/*
   if (condition == false){
      return ();
   }
*/
   std::vector<std::pair<double, double>> child_1(w_dyad.size());
   std::vector<std::pair<double, double>> child_2(w_dyad.size());
   for (i = 0; i < w_dyad.size(); i++){
      child_1[i].first = 0 * w_dyad.size();
      child_1[i].second = 1;
      child_2[i].first = w_dyad[i].first * 2;
      child_2[i].second = w_dyad[i].second;
   }
   std::vector<double> child_1_double(child_1.size());
   std::vector<double> child_2_double(child_2.size());
   for (i = 0; i < w_dyad.size(); i++){
      child_1_double[i] = child_1[i].first / child_1[i].second;
      child_2_double[i] = child_2[i].first / child_2[i].second;
   }
   std::pair<int, int> child_a;
   std::pair<int, int> child_b;
   double bit = 1;
   while (stat_a == true && stat_b == true){
      for (i = 0; i < w_dyad.size(); i++){
         if (child_2_double[i] >= bit){
            child_1_double[i] += bit;
            child_a = fraction(child_1_double[i]);
            child_1[i].first = (double) child_a.first;
            child_1[i].second = (double) child_a.second;
            child_2_double[i] -= bit;
            child_b = fraction(child_2_double[i]);
            child_2[i].first = (double) child_b.first;
            child_2[i].second = (double) child_b.second;
            stat_a = false;
         }
         sum = 0;
         for (j = 0; j < w_dyad.size(); j++){
            sum += child_1_double[j];
         }
         if (sum >= 1 - TOLERANCE && sum <= 1 + TOLERANCE){
            stat_b = false;
         }
         if ((stat_a != stat_b) || (stat_a == true && stat_b == true)){
            stat_a = stat_b = true;
         }
      }
   }

//TODO(fabioftv): Implement bit /= 2
//TODO(fabioftv): Implement error in case of infinite loop

   std::pair<std::vector<std::pair<double, double>>, 
      std::vector<std::pair<double, double>>> split_w_dyad;

   split_w_dyad.first = child_1;
   split_w_dyad.second = child_2;

   return split_w_dyad;
}

double GeoMeanIneq::get_max_denom(std::vector<std::pair<double, double>> tup) {
   double max = 0;
   int i;

   for (i = 0; i < tup.size(); i++){
      if (tup[i].second > max){
         max = tup[i].second;
      }
   }

   return max;
}

unsigned GeoMeanIneq::get_number_of_digits(unsigned number) {
   return number > 0 ? (int) log10 ((double) number) + 1 : 1;
}

int GeoMeanIneq::int_to_binary(int number) {
   int remainder, i = 1, binary = 0;

   while (number != 0){
      remainder = number % 2;
      number /= 2;
      binary += remainder * i;
      i *= 10;
   }

   return binary;
}

double GeoMeanIneq::lower_bound(std::vector<std::pair<double, double>> w_dyad) {
   int i, sum = 0;

   assert (dyadic_weight_vector_test(w_dyad) == true);
   double md = get_max_denom(w_dyad);
   int lb1 = get_number_of_digits(int_to_binary((int) md)) - 1;
   for (i = 0; i < w_dyad.size(); i++){
      if (w_dyad[i].first != 0){
         sum++;
      }
   }
   int lb2 = sum - 1;

   return std::max(lb1, lb2);
}
