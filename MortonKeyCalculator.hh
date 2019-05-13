#ifndef __MORTONKEYCALCULATOR_H__
#define __MORTONKEYCALCULATOR_H__

#include <ostream>
#include "Body.hh"

// compare x coordinates
bool compX(const Body &a, const Body &b);

// compare y coordinates
bool compY(const Body &a, const Body &b);

class MortonKeyCalculator{
	public:
		MortonKeyCalculator(Body *bodies, int N);

		uint32_t x(double xval) const{
			return (uint32_t)(0xffffffff * (xval - xmin) / max_range);
		};

		uint32_t y(double yval) const{
			return (uint32_t)(0xffffffff * (yval - ymin) / max_range);
		}

		// helps calculate max value of width and length for a given level's
		// sub-squares
 		double cellWidth(int level) const{
			return max_range / (1 << (level -1));
		}

		std::ostream& printKey(std::ostream &os, const Body &b) const;

		//returns true if a is less than b in a MortonKey ordering.
		//Note that the y coordinate is more significant.
		bool operator() (const Body &a, const Body &b) const;

private:
	double xmin;
	double ymin;
	double xmax;
	double ymax;
	double max_range;
};

#endif
