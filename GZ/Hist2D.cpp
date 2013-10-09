/*
 *  Hist2D.cpp
 *  GZ
 *
 *  Created by Sven Lafèbre on 29-11-07.
 *
 */

#include "Hist2D.h"

using namespace std;

int Hist2D::bin2mem(int x, int y) const {
	return x*_nbinsy + y;
}

Hist2D::Hist2D() {
	_nbinsx  = 0;
	_nbinsy  = 0;
	_values  = NULL;
	_weights = NULL;
}

Hist2D::Hist2D(const Hist2D &that) {
	_nbinsx  = 0;
	_nbinsy  = 0;

	SetBins(that.GetNbinsX(), that.GetXmin(), that.GetXmax(),
			that.GetNbinsY(), that.GetYmin(), that.GetYmax());
}

Hist2D::Hist2D(int nx, double xmin, double xmax,
			   int ny, double ymin, double ymax) {
	_nbinsx  = 0;
	_nbinsy  = 0;

	SetBins(nx, xmin, xmax,
			ny, ymin, ymax);
}

Hist2D::~Hist2D() {
	Clear();
}

void Hist2D::SetBins(int nx, double xmin, double xmax,
					 int ny, double ymin, double ymax) {
	if (_nbinsx > 0 || _nbinsy > 0) {
		Clear();
	}
	
	if (nx < 1 || ny < 1) {
		return;
	}

	_nbinsx = nx;
	_xmin   = xmin;
	_xmax   = xmax;
	_nbinsy = ny;
	_ymin   = ymin;
	_ymax   = ymax;
	
	_values  = new double[_nbinsx*_nbinsy];
	_weights = new double[_nbinsx*_nbinsy];
	
	for (int x = 0; x < _nbinsx; x ++) {
		for (int y = 0; y < _nbinsy; y ++) {
			_values [bin2mem(x,y)] = 0.;
			_weights[bin2mem(x,y)] = 0.;
	}	}
}

void Hist2D::Clear() {
	delete [] _values;
	delete [] _weights;
	_nbinsx = 0;
	_nbinsy = 0;
}

void Hist2D::Copy(const Hist2D& that) {
	if (this == &that)
		return;

	SetBins(that.GetNbinsX(), that.GetXmin(), that.GetXmax(),
			that.GetNbinsY(), that.GetYmin(), that.GetYmax());

	for (int x = 0; x < _nbinsx; x ++) {
		for (int y = 0; y < _nbinsy; y ++) {
			SetBinValue (x, y, that.GetBinValue (x,y));
			SetBinValue (x, y, that.GetBinValue (x,y));
			SetBinWeight(x, y, that.GetBinWeight(x,y));
			SetBinWeight(x, y, that.GetBinWeight(x,y));
	}	}
}

void Hist2D::operator=(const Hist2D& that) {
	Copy(that);
}

double Hist2D::GetBinValue(int x, int y) const {
	if (x < 0 || y < 0 || x >= _nbinsx || y >= _nbinsy)
		return 0;
	
	return _values[bin2mem(x,y)];
}

double Hist2D::GetValue(double x, double y) const {
	if (x < _xmin || y < _ymin || x > _xmax || y > _ymax)
		return 0;
	
	return _values[bin2mem(GetXbin(x),GetYbin(y))];
}

double Hist2D::GetBinAverage(int x, int y) const {
	if (x < 0 || y < 0 || x >= _nbinsx || y >= _nbinsy)
		return 0;
	
	return _values[bin2mem(x,y)]/
		  _weights[bin2mem(x,y)];
}

double Hist2D::GetAverage(double x, double y) const {
	if (x < _xmin || y < _ymin || x > _xmax || y > _ymax)
		return 0;
	
	return _values[bin2mem(GetXbin(x),GetYbin(y))]/
		  _weights[bin2mem(GetXbin(x),GetYbin(y))];
}

double Hist2D::GetBinWeight(int x, int y) const {
	if (x < 0 || y < 0 || x >= _nbinsx || y >= _nbinsy)
		return 0;

	return _weights[bin2mem(x,y)];
}

double Hist2D::GetWeight(double x, double y) const {
	if (x < _xmin || y < _ymin || x > _xmax || y > _ymax)
		return 0;

	return _weights[bin2mem(GetXbin(x),GetYbin(y))];
}

void Hist2D::SetBinValue(int x, int y, double val) {
	if (x < 0 || y < 0 || x >= _nbinsx || y >= _nbinsy)
		return;
	
	_values[bin2mem(x,y)] = val;
}

void Hist2D::SetValue(double x, double y, double val) {
	if (x < _xmin || y < _ymin || x > _xmax || y > _ymax)
		return;
	
	_values[bin2mem(GetXbin(x),GetYbin(y))] = val;
}

void Hist2D::SetBinWeight(int x, int y, double wght) {
	if (x < 0 || y < 0 || x >= _nbinsx || y >= _nbinsy)
		return;

	_weights[bin2mem(x,y)] = wght;
}

void Hist2D::SetWeight(double x, double y, double wght) {
	if (x < _xmin || y < _ymin || x > _xmax || y >= _ymax)
		return;

	_weights[bin2mem(GetXbin(x),GetYbin(y))] = wght;
}

void Hist2D::AddBin(int x, int y, double val) {
	AddBin(x, y, val, 1.);
}

void Hist2D::AddBin(int x, int y, double val, double wght) {
	if (x < 0 || y < 0 || x >= _nbinsx || y >= _nbinsy)
		return;

	_values [bin2mem(x,y)] += val;
	_weights[bin2mem(x,y)] += wght;
}

void Hist2D::Add(double x, double y, double val) {
	Add(x, y, val, 1.);
}

void Hist2D::Add(double x, double y, double val, double wght) {
	if (x < _xmin || y < _ymin || x > _xmax || y >= _ymax)
		return;

	_values [bin2mem(GetXbin(x),GetYbin(y))] += val;
	_weights[bin2mem(GetXbin(x),GetYbin(y))] += wght;
}

void Hist2D::PrintValues(std::ostream* s) const {
	for (int y = 0; y < _nbinsy; y ++) {
		for (int x = 0; x < _nbinsx; x ++) {
			*s << _values [bin2mem(x,y)] << "\t";
		}
		*s << endl;
	}
	*s << endl << endl;
}

void Hist2D::PrintAverages(std::ostream* s) const {
	for (int y = 0; y < _nbinsy; y ++) {
		for (int x = 0; x < _nbinsx; x ++) {
			if (_weights[bin2mem(x,y)] > 0) {
				*s << _values [bin2mem(x,y)]/
					  _weights[bin2mem(x,y)] << "\t";
			} else {
				*s << 0 << "\t";
			}
		}
		*s << endl;
	}
	*s << endl << endl;
}

void Hist2D::PrintAverageX() {
	for (int x = 0; x < _nbinsx; x ++) {
		double sum = 0;
		for (int y = 0; y < _nbinsy; y ++) {
			sum += _values [bin2mem(x,y)]/
				   _weights[bin2mem(x,y)];
		}
		cout << x << "\t" << sum/_nbinsy << "\t";
	}
	cout << endl;
}

void Hist2D::PrintAverageY() {
	for (int y = 0; y < _nbinsy; y ++) {
		double sum = 0;
		for (int x = 0; x < _nbinsx; x ++) {
			sum += _values [bin2mem(x,y)]/
				   _weights[bin2mem(x,y)];
		}
		cout << y << "\t" << sum/_nbinsx << endl;
	}
	cout << endl;
}
