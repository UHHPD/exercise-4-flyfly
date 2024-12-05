#include "Data.hh"

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

using namespace std;

double f_x(double x){
  double alpha = 0.005;
  double beta = -0.00001;
  double gamma = 0.08;
  double delta = 0.015;

  return alpha + beta * x + gamma * exp(-delta*x);
}

bool Data::compare(Data i, Data j,int bin, int n) const {
  double diff = i.m_data[bin] - j.m_data[bin];
    if (diff ==0) return true;
  double error =  sqrt(pow((i.m_data[bin] - j.m_data[bin])/(abs(diff))* i.m_error[bin],2) + pow((j.m_data[bin] - i.m_data[bin])/(abs(diff))* j.m_error[bin],2));
  return abs(diff) < n*error;
}
Data::Data(const std::string& filename) {
  ifstream file(filename);

  if (!file.is_open()) {
    throw std::runtime_error("Cannot open file: " + filename);
  }

  int size = -1;
  file >> size;

  for (int i = 0; i < size + 1; ++i) {
    double edge;
    file >> edge;
    m_bins.push_back(edge);
  }

  for (int i = 0; i < size; ++i) {
    double entries;
    file >> entries;
    m_data.push_back(entries);
  }

    for (int i = 0; i < size; ++i) {
      double error;
      file >> error;
      m_error.push_back(error);
    }
m_name = filename;

  file.close();

  assertSizes();
};



int Data::checkCompatibility(const Data& in, int n) {
    int count = 0;
    for (int i = 0; i < m_data.size(); ++i) {
        if (!compare(*this,in,i, n)) {
        count++;
        }
    }
    return count;


}
void Data::assertSizes() { assert(m_data.size() + 1 == m_bins.size()); }
Data Data::operator+(Data& in){
Data y; 
   if (this->size() != in.size())
   throw std::invalid_argument("Datasets have different sizes");
   y.m_data.resize(this->size());
   y.m_bins.resize(this->size());
  y.m_bins = this->m_bins;
   y.m_error.resize(this->size());
        for (int i = 0; i < this->size(); ++i){
          double w_1 =1/pow(this->error(i),2);
          double w_2 =1/pow(in.error(i),2);
       y.m_data[i] = (this->measurement(i)*(w_1 )+ in.measurement(i)* w_2)/(w_1 + w_2);
            y.m_error[i] = sqrt(1/(w_1+w_2));
        }
        return y;
    }

double Data::chi_2(){
  double sum = 0;
  for (int i = 0; i < this->size(); ++i){
    sum += pow(this->measurement(i) - f_x(this->binCenter(i)),2)/pow(this->error(i),2);
  }
  return sum/56;
}