#include <iostream>
#include <vector>
#include <functional>
#include <string>

#include "Data.hh"

// generic function comparing two values of some type T used later for int and
// double
template <class T>
bool testEqual(const std::string& name, T expected, T real) {
  if (expected != real) {
    std::cout << "(" << name << ": <" << expected << ">"
              << ", " << real << ") ";
    return false;
  }
  return true;
}

bool testReadingSize() {
  std::cout << "testReadingSize...";
  Data datA("testA");
  return testEqual<int>("size", 1, datA.size());
}

bool testReadingMeasurement() {
  std::cout << "testReadingMeasurement...";
  Data datA("testA");
  return testEqual("measurement", 10., datA.measurement(0));
}

bool testReadingBinEdges() {
  std::cout << "testReadingBinEdges...";
  Data datA("testA");
  return testEqual("bin low", 0., datA.binLow(0)) &&
         testEqual("bin high", 1., datA.binHigh(0));
}

bool testReadingErrors() {
  std::cout << "testReadingErrors...";
  Data datA("testA");
  return testEqual("error", 2., datA.error(0));
}

bool testCopyConstructor() {
  std::cout << "testCopyConstructor...";
  Data a("testA");
  Data b(a);
  Data c = a;  // equivalent to Data c(a)
  return testEqual("measurement", 10., b.measurement(0)) &&
         testEqual("measurement", 10., c.measurement(0));
}

void runTests() {
  std::cout << "running tests...\n";
  std::vector<std::function<bool()> > tests(
      {testReadingSize, testReadingMeasurement, testReadingBinEdges,
       testReadingErrors, testCopyConstructor});
  for (auto test : tests)
    std::cout << (test() ? " ok" : " FAILED!") << std::endl;
}

int main() {
  using namespace std;

  cout << "******************************************************" << endl;
  runTests();
  cout << "******************************************************" << endl;
  // create an object which holds data of experiment A
    std::vector<Data> data;
  Data datA("exp_A");
  Data datB("exp_B");
  Data datC("exp_C");
  Data datD("exp_D");
data.push_back(datA);
data.push_back(datB);
data.push_back(datC);
data.push_back(datD);

  // here is the data from experiment A
  cout << "bin 27: from " << datA.binLow(27) << " to " << datA.binHigh(27)
       << endl;
  cout << "measurement of experiment A in bin 27: " << datA.measurement(27)
       << endl;

  for (Data i : {datA, datB, datC, datD}){

    cout << "numb of incompatible bins: " << i.name() << "and A with 2 sigma:" << i.checkCompatibility(datA, 2) << "" << i.name() << "and B with 2 sigma:" << i.checkCompatibility(datB, 2) << "" << i.name() << "and C with 2 sigma:" << i.checkCompatibility(datC, 2) << "" << i.name() << "and D with 2 sigma:" << i.checkCompatibility(datD, 2) << "" << endl;
    cout << "numb of incompatible bins: " << i.name() << "and A with 3 sigma:" << i.checkCompatibility(datA, 3) << "" << i.name() << "and B with 3 sigma:" << i.checkCompatibility(datB, 3) << "" << i.name() << "and C with 3 sigma:" << i.checkCompatibility(datC, 3) << "" << i.name() << "and D with 3 sigma:" << i.checkCompatibility(datD, 3) << "" << endl;
  }

  cout << "chi_2 of experiment A: " << datA.chi_2() << endl;
  cout << "chi_2 of experiment B: " << datB.chi_2() << endl;
  cout << "chi_2 of experiment C: " << datC.chi_2() << endl;
  cout << "chi_2 of experiment D: " << datD.chi_2() << endl;

  Data sumABCD = datA + datB + datC + datD;
  cout << sumABCD.measurement(27) << endl;
  cout << "chi_2 of experiment A+B+C+D: " << sumABCD.chi_2() << endl;

    Data sumAB = datA + datB;
    cout << "meausurement of A+B: " << sumAB.measurement(27) << endl;

  return 0;
}

