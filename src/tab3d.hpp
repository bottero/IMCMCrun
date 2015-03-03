// tab3d.hpp
// Author : fabien ors
// ENSMP / Armines

#ifndef _TAB3D_H_
#define _TAB3D_H_
#include <iostream>

template <class T> class tab3d;

template <class T> tab3d<T> &operator+(const tab3d<T>& );
template <class T> void max(const tab3d<T>& );
template <class T> void min(const tab3d<T>& );

template<class T>
class tab3d
{
public:
    tab3d();
    tab3d(const tab3d& ref);
    tab3d(int nz, int nx, int ny);
    tab3d(int nz, int nx, int ny, const T& val);
    tab3d(int nz, int nx, int ny, const T *val);
    virtual ~tab3d();

    void reset(int nx, int ny, int nz);
    void fill(const T& val);
    void fill(const T *val);
    tab3d<T> &operator+=(const tab3d<T>& tab);
    tab3d<T> &operator-=(const tab3d<T>& tab);
    tab3d<T> &operator*=(const tab3d<T>& tab);
    tab3d<T> &operator=(const tab3d<T>& tab);
    tab3d<T> &operator+(const tab3d<T> & tab);

    void sqrtt();
    void max( tab3d<T> & tab);
    void min( tab3d<T> & tab);

    void  inverse();
    void  multiply(double val);
    bool set(int iz, int ix, int iy, const T& val);
    const T& get(int iz, int ix, int iy) const;
    int get_nx() const {return _nx;}
    int get_ny() const {return _ny;}
    int get_nz() const {return _nz;}
    T* get_values() const {return _values;}
private:
  void alloc();
  void free();

private:
    T* _values;
    T _udf;
    int _nz;
    int _nx;
    int _ny;
};


template<class T>
tab3d<T>::tab3d() :
  _values(NULL), _udf(), _nz(0) ,_nx(0), _ny(0)
{
}

template<class T>
tab3d<T>& tab3d<T>::operator+=(const tab3d<T>& tab) {
	for(int ix=0; ix < _nx; ix++) {
	    for(int iy=0; iy < _ny; iy++) {
	      for(int iz=0; iz < _nz; iz++) {
	        _values[ix*_ny*_nz+iy*_nz+iz] += tab._values[ix*_ny*_nz+iy*_nz+iz];

	      }
	    }
	  }

	return *this;

};

template<class T>
tab3d<T>& tab3d<T>::operator-=(const tab3d<T>& tab) {
	for(int ix=0; ix < _nx; ix++) {
	    for(int iy=0; iy < _ny; iy++) {
	      for(int iz=0; iz < _nz; iz++) {
	        _values[ix*_ny*_nz+iy*_nz+iz] -= tab._values[ix*_ny*_nz+iy*_nz+iz];

	      }
	    }
	  }

	return *this;

};


template<class T>
tab3d<T>& tab3d<T>::operator*=(const tab3d<T>& tab) {
	for(int ix=0; ix < _nx; ix++) {
	    for(int iy=0; iy < _ny; iy++) {
	      for(int iz=0; iz < _nz; iz++) {
	        _values[ix*_ny*_nz+iy*_nz+iz] *= tab._values[ix*_ny*_nz+iy*_nz+iz];

	      }
	    }
	  }

	return *this;

};

template<class T>
tab3d<T> &tab3d<T>::operator+(const tab3d<T>& tab)
    {tab3d<T> tabres=*this;
    tabres+=tab;
    return tabres;
    }


template<class T>
tab3d<T>&tab3d<T>::operator=(const tab3d<T>& tab)
{
  if(this!=&tab)
  {
    free();
    _nx=tab._nx;
    _ny=tab._ny;
    _nz=tab._nz;
    alloc();
    for(int ix=0; ix < _nx; ix++) {
      for(int iy=0; iy < _ny; iy++) {
        for(int iz=0; iz < _nz; iz++) {
          _values[ix*_ny*_nz+iy*_nz+iz] =tab._values[ix*_ny*_nz+iy*_nz+iz];
        }
      }
    }
  }
  return *this;
}



template<class T>
void tab3d<T>::sqrtt()
{for(int ix=0; ix < _nx; ix++) {
    for(int iy=0; iy < _ny; iy++) {
      for(int iz=0; iz < _nz; iz++) {

       _values[ix*_ny*_nz+iy*_nz+iz] =sqrt(_values[ix*_ny*_nz+iy*_nz+iz]);

      }
    }
  }


}

template<class T>
void tab3d<T>::max( tab3d<T> & tab)
{
	for(int ix=0; ix < _nx; ix++) {
	    for(int iy=0; iy < _ny; iy++) {
	      for(int iz=0; iz < _nz; iz++) {
	       if(tab._values[ix*_ny*_nz+iy*_nz+iz]>_values[ix*_ny*_nz+iy*_nz+iz])
	    	   {_values[ix*_ny*_nz+iy*_nz+iz]= tab._values[ix*_ny*_nz+iy*_nz+iz];}

	      }
	    }
	  }



}

template<class T>
void tab3d<T>::min( tab3d<T> & tab)
{
	for(int ix=0; ix < _nx; ix++) {
	    for(int iy=0; iy < _ny; iy++) {
	      for(int iz=0; iz < _nz; iz++) {
	       if(tab._values[ix*_ny*_nz+iy*_nz+iz]<_values[ix*_ny*_nz+iy*_nz+iz])
	    	   {_values[ix*_ny*_nz+iy*_nz+iz]= tab._values[ix*_ny*_nz+iy*_nz+iz];}

	      }
	    }
	  }



}


template<class T>
void tab3d<T>::inverse() {
	for(int ix=0; ix < _nx; ix++) {
	    for(int iy=0; iy < _ny; iy++) {
	      for(int iz=0; iz < _nz; iz++) {
	       _values[ix*_ny*_nz+iy*_nz+iz] =1/ _values[ix*_ny*_nz+iy*_nz+iz];

	      }
	    }
	  }
}


template<class T>
void tab3d<T>::multiply(double val) {
	for(int ix=0; ix < _nx; ix++) {
	    for(int iy=0; iy < _ny; iy++) {
	      for(int iz=0; iz < _nz; iz++) {
	       _values[ix*_ny*_nz+iy*_nz+iz] *=val;

	      }
	    }
	  }
}

template<class T>
tab3d<T>::tab3d(const tab3d& ref) :
  _values(NULL), _udf(), _nz(ref._nz), _nx(ref._nx), _ny(ref._ny) {
  alloc();
  fill(ref._values);
}

template<class T>
tab3d<T>::tab3d(int nz, int nx, int ny) :
  _values(NULL), _udf(), _nz(nz), _nx(nx), _ny(ny)
{
  alloc();
}

template<class T>
tab3d<T>::tab3d(int nz, int nx, int ny, const T& val) :
  _values(NULL), _udf(), _nz(nz), _nx(nx), _ny(ny)
{
  alloc();
  fill(val);
}


template<class T>
tab3d<T>::tab3d(int nz, int nx, int ny, const T *val) :
  _values(NULL), _udf(), _nz(nz), _nx(nx), _ny(ny)
{
  alloc();
  fill(val);
}


template<class T>
tab3d<T>::~tab3d()
{
  free();
}

template<class T>
void tab3d<T>::reset(int nz, int nx, int ny)
{
  _nx = nx;
  _ny = ny;
  _nz = nz;
  alloc();
}

template<class T>
void tab3d<T>::fill(const T& val)
{
  T* pval = _values;
  for (int ii=0; ii < _nx*_ny*_nz; ii++) {
    *pval = val;
    pval++;
  }
}

template<class T>
void tab3d<T>::fill(const T *val)
{
  T* pval = _values;
  for (int ii=0; ii < _nx*_ny*_nz; ii++) {
    *pval = val[ii];
    pval++;
  }
}

template<class T>
bool tab3d<T>::set(int iz, int ix, int iy, const T& val)
{
  if(ix < 0 || ix >= _nx ||
     iy < 0 || iy >= _ny ||
     iz < 0 || iz >= _nz) {
    return false;
  }

  _values[iy*_nx*_nz+ix*_nz+iz] = val;
  return true;
}

template<class T>
const T& tab3d<T>::get(int iz, int ix, int iy) const
{
  if(ix < 0 || ix >= _nx ||
     iy < 0 || iy >= _ny ||
     iz < 0 || iz >= _nz) {
    return _udf;
  }
  return _values[iy*_nx*_nz+ix*_nz+iz];
}

template<class T>
void tab3d<T>::alloc()
{

  if(_values != NULL)
    free();

  _values = new T [_nx*_ny*_nz];

}

template<class T>
void tab3d<T>::free()
{
  if(_values != NULL) {
    delete[] _values;
  }
  _values = NULL;
}

#endif

