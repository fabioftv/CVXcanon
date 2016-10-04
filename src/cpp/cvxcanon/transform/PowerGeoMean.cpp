#include "cvxcanon/transform/PowerGeoMean.hpp"

#include <vector>
#include <math.h>
#include <iostream>
#include <assert.h>
#include <algorithm>

#include "glog/logging.h"

#define TOLERANCE 0.00000001
#define PRECISION 1000000000

GeoMeanIneq::GeoMeanIneq() {}
GeoMeanIneq::~GeoMeanIneq() {}

//The gcd function calculates the greater common denominator between two numbers.
long GeoMeanIneq::gcd(long a, long b) {
   if (a == 0){
      return b;
   }
   else if (b == 0){
      return a;
   }
   else if (a < b){
      return gcd(a, b % a);
   }
   else{
      return gcd(b, a % b);
   }
}

//The function fraction converts a decimal number into a fraction. This fraction is represented as a pair where the first number of the pair is the numerator and the second number of the pair is the denominator.
std::pair<double, double> GeoMeanIneq::fraction(double number) {
   long gcd_aux, denom, num;
   int i;
   if (number >= 0){
      gcd_aux = gcd(round(number * PRECISION), PRECISION);
      denom = PRECISION / gcd_aux;
      num = round(number * PRECISION) / gcd_aux;
   }
   else {
      gcd_aux = gcd(round(-number * PRECISION), PRECISION);
      denom = PRECISION / gcd_aux;
      num = -round(-number * PRECISION) / gcd_aux;
   }
   return std::make_pair((double)num, (double)denom);
}

//The function build_power_pos considers p > 1. The function returns the power tuple (t,1,x) in which x <= t^(1/p)*1^(1-1/p).
std::pair<std::pair<double, double>, std::pair<std::pair<double, double>, 
   std::pair<double, double>>> GeoMeanIneq::build_power_pos(double number) {
   
   //If p <= 1, then exit.
   assert(number > 1);   
   
   std::pair<double, double> frac, frac_inv, frac_minus;
   
   frac = fraction(number);
   frac_inv.first = frac.second;
   frac_inv.second = frac.first;
   frac_minus.first = frac_inv.second - frac_inv.first;
   frac_minus.second = frac_inv.second;

   //Return (p, (1/p, 1-1/p)).
   return std::make_pair(frac, std::make_pair(frac_inv, frac_minus));
}

//The function build_power_neg considers p < 0. The function returns the power tuple (x,t,1) in which 1 <= x^(p/(p-1))*t^(-1/(p-1)).
std::pair<std::pair<double, double>, std::pair<std::pair<double, double>, 
   std::pair<double, double>>> GeoMeanIneq::build_power_neg(double number) {

   //If p >= 0, then exit.   
   assert(number < 0);
   
   std::pair<double, double> frac, frac_div, frac_minus;
   
   frac = fraction(number);
   frac_div.first = -frac.first;
   frac_div.second = (-frac.first) + frac.second;
   frac_minus.first = frac_div.second - frac_div.first;
   frac_minus.second = frac_div.second;

   //Return (p, (p/p-1, 1-p/p-1)).
   return std::make_pair(frac, std::make_pair(frac_div, frac_minus));
}

//The function build_power_middle considers 0 < p < 1. The function returns the power tuple (x,1,t) in which t <= x^p*1^(1-p).
std::pair<std::pair<double, double>, std::pair<std::pair<double, double>, 
   std::pair<double, double>>> GeoMeanIneq::build_power_middle(double number) {

   //If p <= 0 or p >= 1, then exit.     
   assert(number > 0 && number < 1);
   
   std::pair<double, double> frac, frac_minus;
   
   frac = fraction(number);
   frac_minus.first = frac.second - frac.first;
   frac_minus.second = frac.second;

   //Return (p, (p, 1-p)).
   return std::make_pair(frac, std::make_pair(frac, frac_minus));
}

//The function is_integer determines if a number is integer. In such a case, 4 is an integer number as the same way as 4.0. Note that this verification process is bounded by a tolerance so 4.00000001 is not an integer number but 4.000000001 is an integer.
bool GeoMeanIneq::is_integer(double number) {
   return (number > (floor(number)-TOLERANCE) && number < 
             (floor(number)+TOLERANCE));
}

//The function power_two_test determines if a number is a positive integer number power of 2.
bool GeoMeanIneq::power_two_test(double number) {
   return (is_integer(number) == true && number > 0 && (! ((int) number & 
             ((int) number -1))) == true);
}

//The function dyadic_nonnegative_fraction_test determines if a fraction is a nonnegative dyadic fraction or integer. The fraction is represented by a pair where the first number of the pair is the numerator and the second number of the pair is the denominator.
bool GeoMeanIneq::dyadic_nonnegative_fraction_test(std::pair<double, double>
     fraction) {
   //First, the function determines whether the fraction is an actual integer 
   //number >= 0. For instance, the number 4 must be represented by the 
   //pair (4,1). 
   if ((fraction.first / fraction.second) == fraction.first && 
       is_integer(fraction.first) == true &&
       fraction.first >= 0){
      return true;
   }
   //Second, the function determines whether the fraction is a real fraction
   //(e.g.: 1/3 represented by the pair (1,3)) where both the numerator and 
   //denominator are integer numbers and the fraction is >= 0.
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

//The function weight_vector_test determines if a vector is a valid weight vector. In such a case, all elements of the vector must be nonnegative integers or fractions, and they must sum to 1. Each element of the vector is represented by a pair where the first number of the pair is the numerator and the second number of the pair is the denominator.
bool GeoMeanIneq::weight_vector_test(std::vector<std::pair<double, double>> w){
   double sum = 0;   
   bool aux = true;
   int i = 0;
   //The while loop evaluates all the elements of the vector until one of the 
   //elements violates the conditions or it reaches the end of the vector.
   while (i < w.size()){
      //First, the loop checks whether the corresponding element is an actual
      //integer number >= 0 (e.g.: 4 represented by the pair (4,1)) or a
      //fraction of integers >= 0 (e.g.: 1/3 represented by the pair (1/3)).
      if ((w[i].first / w[i].second) == w[i].first){
         if (is_integer(w[i].first) == true && w[i].first >= 0){
            sum += w[i].first;
         }
         else{
            //The loop terminates because it violates the condition.
            aux = false;
            i = w.size(); 
         }
      }
      else {
         if (is_integer(w[i].first) == true && is_integer(w[i].second) == true
             && w[i].first / w[i].second >= 0){
            sum += w[i].first / w[i].second;
         }
         else{
            //The loop terminates because it violates the condition.
            aux = false;
            i = w.size();          
         }
      }
      i++;
      //The loop also terminates when the end of the vector is not reached but
      //the sum of the elements is >= 1.
      //is >= 1 or 
      if (sum > (1 + TOLERANCE) && i <= w.size()){
         aux = false;
         i = w.size();
      }
   }
   return aux;
}

// The function dyadic_weight_vector_test determines if a vector is a valid dyadic weight vector. In such a case, all elements of the vector must be nonnegative integers or dyadic fractions, and they must sum to 1. Each element of the vector is represented by a pair where the first number of the pair is the numerator and the second number of the pair is the denominator.
bool GeoMeanIneq::dyadic_weight_vector_test(std::vector<std::pair<double, 
                                            double>> w){
   bool aux = true;
   int i = 0;
   //The while loop evaluates all the elements of the vector until one of the 
   //elements violates the conditions or it reaches the end of the vector.
   while (i < w.size()){
      if (dyadic_nonnegative_fraction_test(w[i]) != true){
         aux = false;
         i = w.size();  
      }
      i++;
   }
   return (aux == true && weight_vector_test(w) == true);
}

//The function sort returns the indices of a sorted vector.
std::vector<int> GeoMeanIneq::sort(std::vector<double> test){
   std::vector<int> y(test.size());
   std::size_t n(0);
   std::generate(std::begin(y), std::end(y), [&]{return n++;});
   std::sort (std::begin(y), std::end(y), [&](int i1, int i2) 
      {return test[i1] < test[i2];});
   return y;
}

std::pair<std::vector<std::pair<double, double>>, std::vector<std::pair<double, 
   double>>> GeoMeanIneq::fracify(std::vector<std::pair<double, double>> a,
      double max_denom, bool force_dyad) {

   bool aux_a= true, aux_b = true;
   double sum = 0, sum_w_frac = 0;
   int i;

   std::vector<std::pair<double, double>> w_frac(a.size());
   std::vector<double> a_double(a.size());

   for (i = 0; i < a.size(); i++){
      sum += a[i].first / a[i].second;
      a_double[i] = a[i].first / a[i].second;
      if ((a[i].first < 0 || a[i].second < 0) && aux_a == true){
         aux_a = false;
      }
      if ((is_integer(a[i].first) != true || is_integer(a[i].second) != true) && 
         aux_b == true){
         aux_b = false;
      }
   }
   assert (aux_a);
   //TODO(fabioftv): Implement error message in case aux_a == false. "Input
   //powers must be nonnegative".
   assert (is_integer(max_denom) == true && max_denom > 0);
   //TODO(fabioftv): Implement error message in case max_denom <= 0 
   //and not integer. "Input senominator must be an integer".

   max_denom = next_power_two(max_denom);

   if (force_dyad == true){
      w_frac = make_frac(a_double, max_denom);
   }
   else if (aux_b == true){
      for (i = 0; i < w_frac.size(); i++){
         w_frac[i].first = a[i].first;
         w_frac[i].second = a[i].second * sum;
      }
      double d = get_max_denom(w_frac);
      if (d > max_denom) {
      //TODO(fabioftv): Implement error meassage in case d > max_denom. "Cannot
      //reliably represent the input weight vector".
      }
   }
   else{
      for (i = 0; i < w_frac.size(); i++){
         w_frac[i].first = a[i].first;
         w_frac[i].second = a[i].second * sum;
         sum_w_frac += w_frac[i].first / w_frac[i].second;
      }
      if (sum_w_frac > 1+TOLERANCE || sum_w_frac < 1-TOLERANCE){
         w_frac = make_frac(a_double, max_denom);
      }
   }

   std::pair<std::vector<std::pair<double, double>>, 
      std::vector<std::pair<double, double>>> tup;
   tup.first = w_frac;
   tup.second = dyad_completion(w_frac);

   return tup;
}

//The function make_frac approximates a/sum(a) with a tuple of fractions with an exact denominator. In such a case, a is a vector with decimal elements. The function returns a vector of pairs where the first number of the pair is the numerator and the second number of the pair is the denominator.
std::vector<std::pair<double, double>> GeoMeanIneq::make_frac(std::vector<double>
                                 a, int denominator) {
   double sum_a = 0;
   int i = 0, sum_b = 0;
   std::vector<int> b(a.size()), index_sort(a.size());
   std::vector<double> error(a.size());

   for (i = 0; i < a.size(); i++){
      sum_a += a[i];
   }
   for (i = 0; i < a.size(); i++){
      a[i] = a[i] / sum_a;
      b[i] = (int) (a[i] * denominator);
      sum_b += b[i];
      error[i] = ((double) b[i] /  denominator) - a[i];
   }

   index_sort = sort(error);
   std::vector<int> inds(denominator - sum_b);
   for (i = 0; i < (denominator - sum_b); i++){
      inds[i] = index_sort[i];
      b[inds[i]] = b[inds[i]] + 1;
   }
   std::vector<std::pair<double, double>> frac(a.size());
   for (i = 0; i < frac.size(); i++){
      frac[i] = std::make_pair((double) b[i], (double) denominator);
   }
   return frac;
}

//The function dyad_completion returns the dyadic completion of a vector. Every element of the input vector has nonnegative fractions or integers that sum to 1.
std::vector<std::pair<double, double>> 
   GeoMeanIneq::dyad_completion(std::vector<std::pair<double, double>> w) {
   double div = 0;
   int i;

   double d = get_max_denom(w);
   int p = next_power_two((int) d);

   if ((double) p == d){
      return w;
   }
   else {
      std::vector<std::pair<double, double>> temp(w.size() + 1);
      for (i = 0; i < w.size(); i++){
         div = ((w[i].first) * ((double) d)) / 
               ((w[i].second) * ((double) p));
         std::pair<double, double> frac;
         frac = fraction(div);
         temp[i].first = frac.first;
         temp[i].second = frac.second;
      }
      temp[w.size()].first = (double) p - d;
      temp[w.size()].second = (double) p;
      return temp;
   }
}

//The function approx_error returns the norm error from approximating the vector a_orig/sum(a_orig) with the weight vector w_approx.
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

//The function next_power_two returns the next power of two from a given number.
int GeoMeanIneq::next_power_two(int number){
   if (number <= 0){
      return 1;
   }
   else{
      int next = pow(2, ceil(log(number) / log(2)));
      return next;
   }
}

//The function check_dyad checks that w_dyad is a valid dyadic completion of w.
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

//The function split, splits a tuple of dyadic rationals into two children so that the returned function represents 1/2*(child_1 + child_2).
std::pair<std::vector<std::pair<double, double>>, std::vector<std::pair<double,
   double>>> GeoMeanIneq::split(std::vector<std::pair<double, double>> w_dyad) {
   bool condition = true, stat_a = true, stat_b = true;
   double sum = 0;
   int i, j;

   std::pair<std::vector<std::pair<double, double>>, 
      std::vector<std::pair<double, double>>> split_w_dyad;
   std::vector<std::pair<double, double>> child_1(w_dyad.size());
   std::vector<std::pair<double, double>> child_2(w_dyad.size());

   for (i = 0; i < w_dyad.size(); i++){
      if (w_dyad[i].first == w_dyad[i].second){
         condition = false;
         i = w_dyad.size();
      }
   }
   //TODO(fabioftv): Need to return empty.
   if (condition == false){
      for (i = 0; i < w_dyad.size(); i++){
         child_1[i].first = 0;
         child_1[i].second = 0;
         child_2[i].first = 0;
         child_2[i].second = 0;
      }
      split_w_dyad.first = child_1;
      split_w_dyad.second = child_2;
      return split_w_dyad;
   }
   else {
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
      //TODO(fabioftv): Implement bit /= 2.
      //TODO(fabioftv): Implement error message in case of infinite loop.
      split_w_dyad.first = child_1;
      split_w_dyad.second = child_2;
      return split_w_dyad;
   }
}

std::vector<std::pair<std::vector<std::pair<double, double>>, 
   std::pair<std::vector<std::pair<double, double>>, 
      std::vector<std::pair<double, double>>>>> GeoMeanIneq::decompose
         (std::vector<std::pair<double, double>> w_dyad) {

}

//The function get_max_denom returns the max denominator of a vector of fractions.
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

//The function get_number_of_digits determines the number of digits of a given number.
unsigned GeoMeanIneq::get_number_of_digits(unsigned number) {
   return number > 0 ? (int) log10 ((double) number) + 1 : 1;
}

//The function int_to_binary converts a decimal number to a binary number.
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

//The function lower_bound returns a lower bound on the number of cones needed to represent the tuple based on two simple lower bounds.
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

