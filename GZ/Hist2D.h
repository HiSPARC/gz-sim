/*
 *  Hist2D.h
 *  GZ
 *
 *  Created by Sven Lafèbre on 29-11-07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __HIST2D_H__
#define __HIST2D_H__ 1

#include <iostream>

///
/// \brief	Two-dimensional histogram class
///
/// This class describes a two-dimensional histogram of data type <tt>double</tt>
class Hist2D {
	int _nbinsx;		///< Number of bins in histogram in \a x direction
	double _xmin;		///< Lower edge of first bin in \a x direction
	double _xmax;		///< Upper edge of last bin in \a x direction
	
	int _nbinsy;		///< Number of bins in histogram in \a y direction
	double _ymin;		///< Lower edge of first bin in \a y direction
	double _ymax;		///< Upper edge of last bin in \a y direction
	
	double* _values;    ///< Histogram data values
	double* _weights;   ///< Histogram data weights

	/// \brief	Convert bin coordinates to array element number
	///
	/// Convert a set of \a x and \a y bin number to an element number for use in the
	/// ::_values and ::_weights arrays. It is private, as GetXbin(), GetYbin()
	/// and GetBinValue() supply all the access the user needs.
	/// \param	x	\a x bin number [0,<tt>_nbinsx</tt>)
	/// \param	y	\a y bin number [0,<tt>_nbinsy</tt>)
	/// \return	Element to be used in #_values and #_weights array
	int bin2mem(int x, int y) const;

public:
	/// \brief	Default constructor
	///
	/// Create an empty instance of Hist2D. No memory is allocated; use SetBins() or Copy()
	/// to do this at a later stage.
	Hist2D();
	
	/// \brief	Copy constructor
	///
	/// Create a new instance of Hist2D. Histogram contents and boundaries are equal to those
	/// of the histogram copied.
	/// \param	src	Source histogram to copy
	Hist2D(const Hist2D& src);
	
	/// \brief	Constructor
	///
	/// Create a new instance of Hist2D. Histogram boundaries are set to those
	/// supplied to the function and memory is allocated.
	/// \param	nx		Desired number of bins in \a x direction
	/// \param	xmin	Desired value of lower edge of first bin
	/// \param	xmax	Desired value of upper edge of last bin
	/// \param	ny		Desired number of bins in \a y direction
	/// \param	ymin	Desired value of lower edge of first bin
	/// \param	ymax	Desired value of upper edge of last bin
	Hist2D(int nx, double xmin, double xmax,
		   int ny, double ymin, double ymax);
	
	/// \brief	Destructor
	~Hist2D();
	
	/// \brief	Initialize histogram
	///
	/// Initialize a histogram. Histogram boundaries are set to those
	/// supplied to the function and memory is allocated.
	/// \param	nx		Desired number of bins in \a x direction
	/// \param	xmin	Desired value of lower edge of first bin
	/// \param	xmax	Desired value of upper edge of last bin
	/// \param	ny		Desired number of bins in \a y direction
	/// \param	ymin	Desired value of lower edge of first bin
	/// \param	ymax	Desired value of upper edge of last bin
	void SetBins(int nx, double xmin, double xmax,
				 int ny, double ymin, double ymax);
				 
	/// \brief	Clear histogram contents and free memory
	///
	/// Clear histograms #_values and #_weights by releasing the memory. The number of bins
	/// #_nbinsx and #_nbinsy are set to zero.
	void Clear();
	
	/// \brief	Copy a histogram's contents
	///
	/// Copy the contents and geometry of a histogram to the current histogram. The current
	/// histogram is ::Clear()ed first.
	/// \param	src	Source histogram to copy from
	void Copy(const Hist2D& src);
	
	/// \brief	Copy a histogram's contents by assignment
	///
	/// Copy the contents and geometry of a histogram to the current histogram by using the
	/// assignment operator. The current histogram is ::Clear()ed first.
	/// \param	src	Source histogram to copy from
	void operator= (const Hist2D& src);
	
	/// \brief	Get number of bins of \a x direction
	///
	/// \return #_nbinsx
	int GetNbinsX()  const {return _nbinsx;}

	/// \brief	Get minimum value of \a x direction
	///
	/// \return #_xmin
	double GetXmin() const {return _xmin;}

	/// \brief	Get maximum value of \a x direction
	///
	/// \return #_xmax
	double GetXmax() const {return _xmax;}

	/// \brief	Get number of bins of \a y direction
	///
	/// \return #_nbinsy
	int GetNbinsY()  const {return _nbinsy;}

	/// \brief	Get minimum value of \a y direction
	///
	/// \return #_ymin
	double GetYmin() const {return _ymin;}

	/// \brief	Get maximum value of \a y direction
	///
	/// \return #_ymax
	double GetYmax() const {return _ymax;}

	/// \brief	Convert coordinate in \a x direction to bin in \a x direction
	///
	/// \param	x	Coordinate to convert
	/// \return	Bin number in \a x direction
	int GetXbin(double x) const {return int((x-_xmin)/(_xmax-_xmin)*_nbinsx); }

	/// \brief	Convert coordinate in \a y direction to bin in \a y direction
	///
	/// \param	y	Coordinate to convert
	/// \return	Bin number in \a y direction
	int GetYbin(double y) const {return int((y-_ymin)/(_ymax-_ymin)*_nbinsy); }

	/// \brief	Get value of a bin
	double GetBinValue(int, int) const;

	/// \brief	Get value of a bin at a coordinate
	double GetValue(double, double) const;

	/// \brief	Get average of a bin
	double GetBinAverage(int, int) const;

	/// \brief	Get average of a bin at a coordinate
	double GetAverage(double, double) const;

	/// \brief	Get weight of a bin
	double GetBinWeight(int, int) const;

	/// \brief	Get weight of a bin at a coordinate
	double GetWeight(double, double) const;

	/// \brief	Set value of a bin
	void SetBinValue(int, int, double);

	/// \brief	Set value of a bin at a coordinate
	void SetValue(double, double, double);

	/// \brief	Set weight of a bin
	void SetBinWeight(int, int, double);

	/// \brief	Set weight of a bin at coordinate
	void SetWeight(double, double, double);
	
	/// \brief	Add a value to a bin
	void AddBin(int, int, double);

	/// \brief	Add a value with a certain weight to a bin
	void AddBin(int, int, double, double);
	
	/// \brief	Add a value to a bin at a coordinate
	void Add(double, double, double);

	/// \brief	Add a value with a certain weight to a bin at a coordinate
	void Add(double, double, double, double);
	
	/// \brief	Stream out a matrix of histogram values
	///
	/// Stream out an ascii matrix of the histogram values. Each line contains one row of
	/// comma separated values for a single \a y value.
	///
	/// The output of this function can, for instance, be used by the gnuplot program with the
	/// \c matrix option.
	void PrintValues  (std::ostream*) const;

	/// \brief	Stream out a matrix of histogram average values
	///
	/// Stream out an ascii matrix of the histogram values devided by their weights. Each line
	/// contains one row of comma separated \a values for a single \a y value.
	///
	/// The output of this function can, for instance, be used by the gnuplot program with the
	/// \c matrix option.
	void PrintAverages(std::ostream*) const;
	
	void PrintAverageX();
	void PrintAverageY();
};

#endif // __HIST2D_H__
