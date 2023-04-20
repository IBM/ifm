#ifndef __WATCH_H
#define __WATCH_H

#include <string>
#include <list>
#include "ifm_common.h"

using namespace std;

class WatchPoint {
  public:
    WatchPoint() : x(0), y(0), _name(""), _reported(false), _maxSize(10) {
      _history.clear();
    }
    ~WatchPoint() { }

    void init(int xx, int yy, string nn) {
      x = xx;
      y = yy;
      _name = nn;
    }

    void report(int kk) {
      if (! _reported) 
        printf("Watch: %s is draining after %d seconds\n", _name.c_str(), kk);
      _reported = true;
    }

    void push(real_t value) {
      _history.push_back(value);
      if (_history.size() == _maxSize)
        _history.pop_front();
    }

    bool isDraining() {
      if (_history.size() == _maxSize-1) {
        real_t last = _history.front() + 1;
        list<real_t>::iterator it;
        for (it=_history.begin(); it!=_history.end(); it++) {
          if (*it > last)
            return false;
          last = *it;
        }
        return true;
      }
      return false;
    }

  // public variables
  uint64_t x;
  uint64_t y;

  private:
    string _name;
    bool _reported;
    unsigned int _maxSize;
    list<real_t> _history;
};

#endif
