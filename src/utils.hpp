#pragma once
#include <cassert>
#include <iostream>
#include <string>
#include <sys/stat.h>

using namespace std;

string time_to_str(time_t time_taken);

int binomial(int n, int k);

double smart_round(double d);
double smart_round(double d, double precision);

bool create_dir(string path);

bool starts_with_flag_expression(string filename);

//  The following is by James Kanze
// http://gabisoft.free.fr/articles/fltrsbf1.html
// https://openclassrooms.com/forum/sujet/fluxtampon-gestion-31091?page=1
// https://lists.boost.org/Archives/boost/att-49459/fltrsbf1.htm
struct FilteringInputStreambuf : public std::streambuf {
  FilteringInputStreambuf(streambuf *source, bool deleteWhenFinished = false)
      : mySource(source), myDeleteWhenFinished(deleteWhenFinished) {}

  virtual ~FilteringInputStreambuf() { resetSource(NULL); }

  void resetSource(streambuf *newSource, bool deleteWhenFinished = false) {
    sync();

    if (myDeleteWhenFinished)
      delete mySource;

    mySource = newSource;
    myDeleteWhenFinished = deleteWhenFinished;
    setg(NULL, NULL, NULL);
  }

  virtual int underflow() {
    int result(EOF);

    if (gptr() < egptr())
      result = static_cast<unsigned char>(*gptr());
    else if (mySource != NULL) {
      result = extract(*mySource);

      if (result != EOF) {
        assert(
            result >=
            0); // && result <= UCHAR_MAX); // does not work on older machines
        myPushbackBuffer = result;
        setg(&myPushbackBuffer, &myPushbackBuffer, &myPushbackBuffer + 1);
      }
    }

    return result;
  }
  virtual int sync() {
    int result(EOF);

    if (mySource != NULL) {
      if (gptr() == egptr() || mySource->sputbackc(*gptr()) != EOF)
        result = mySource->pubsync();

      setg(NULL, NULL, NULL);
    }

    return result;
  }
  virtual std::streambuf *setbuf(char *buffer, std::streamsize length) {
    return mySource == NULL ? static_cast<std::streambuf *>(NULL)
                            : mySource->pubsetbuf(buffer, length);
  }

  virtual int extract(std::streambuf &source) {
    int ch(source.sbumpc());
    if (ch == '#') {
      while (ch != EOF && ch != '\n' && ch != '\r') {
        ch = source.sbumpc();
      }
    }
    if (ch >= 126) {
      cerr << "WARNING: Reading from input: " << ch << " which is '" << char(ch)
           << "'" << endl;
    }
    return ch;
  }

  std::streambuf *mySource;
  char myPushbackBuffer;
  bool myDeleteWhenFinished;
};

class FilteringIstream : private FilteringInputStreambuf, public istream {
public:
  FilteringIstream(istream &source)
      : FilteringInputStreambuf(source.rdbuf()), istream(this) {}
  virtual ~FilteringIstream() {};
  virtual FilteringInputStreambuf *rdbuf() { return this; }
  void changeSource(streambuf *newSource) { resetSource(newSource); }
};

// This is a macro that creates  istream *istr variable that opened the filename
// in a smart and nice way.
//
//     istream *istr = &std::cin;
#define OPEN_FILE_SMARTLY_BASE(istr, filename, fail_operation)                 \
  FilteringIstream fscin_XXX(std::cin);                                        \
  istream *istr = &fscin_XXX;                                                  \
  stringstream ss_filename_XXX(filename);                                      \
  ifstream infile_XXX;                                                         \
  bool file_exists_XXX = false;                                                \
  if (filename != "cin") {                                                     \
    infile_XXX.open(filename.c_str(), ifstream::in);                           \
    if (!infile_XXX.good()) {                                                  \
      std::cerr << "Failed opening file " << filename << endl;                 \
      if (starts_with_flag_expression(filename)) {                             \
        std::cerr << "Using it as the string itself" << endl;                  \
        istr = &ss_filename_XXX;                                               \
      } else {                                                                 \
        std::cerr << "Failed interpreting as flags without comments" << endl;  \
        fail_operation;                                                        \
      }                                                                        \
    } else {                                                                   \
      istr = &infile_XXX;                                                      \
      file_exists_XXX = true;                                                  \
    }                                                                          \
  }                                                                            \
  FilteringIstream infileFilter_XXX(infile_XXX);                               \
  if (filename != "cin" && file_exists_XXX) {                                  \
    istr = &infileFilter_XXX;                                                  \
  }

#define OPEN_FILE_SMARTLY_RETURN_FALSE_ON_FAIL(istr, filename)                 \
  OPEN_FILE_SMARTLY_BASE(istr, filename, return false)
#define OPEN_FILE_SMARTLY(istr, filename)                                      \
  OPEN_FILE_SMARTLY_BASE(istr, filename, assert(0))
