#ifndef DATA_HH
#define DATA_HH

#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>

class Data {
 public:
  explicit Data(const std::string& filename);

  unsigned int size() const { return m_data.size(); }
  double measurement(int i) const { return m_data[i]; }
  double binCenter(int i) const { return 0.5 * (m_bins[i] + m_bins[i + 1]); }
  double binLow(int i) const { return m_bins.front();}
  double binHigh(int i) const { return m_bins.back();}
  double error(int i) const { return m_error[i]; }

std::string name () {return m_name;}
  int checkCompatibility(const Data& in, int n);
 Data operator+( Data& in) ;
double chi_2();
 private:
            Data() {}  // disallow empty data set
  void assertSizes();
  bool compare(Data i, Data j,int bin, int n) const;
  std::vector<double> m_data;
  std::vector<double> m_error;
  std::vector<double> m_bins;
std::string m_name;
};

#endif
