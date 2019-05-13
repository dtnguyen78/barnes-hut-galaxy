#ifdef _OPENMP
#include <parallel/algorithm>
#include <parallel/numeric>
#else
#include <algorithm>
#endif

#include "MortonKeyCalculator.hh"

bool compX(const Body &a, const Body &b)
{
	return a.x() < b.x();
}

bool compY(const Body &a, const Body &b)
{
	return a.y() < b.y();
}

MortonKeyCalculator::MortonKeyCalculator(Body *bodies, int N){
	// Find in parallel the minimum x and y values among bodies and
	// the range of values (i.e. max(xmax - xmin, ymax - ymin)).

	// Use min_element and max_element functions.
	xmin = __gnu_parallel::min_element(bodies, bodies+N, compX)->x();
	xmax = __gnu_parallel::max_element(bodies, bodies+N, compX)->x();
	ymin = __gnu_parallel::min_element(bodies, bodies+N, compY)->y();
	ymax = __gnu_parallel::max_element(bodies, bodies+N, compY)->y();

	this->max_range = std::max(xmax - xmin, ymax - ymin);
};

std::ostream& MortonKeyCalculator::printKey(std::ostream &os, const Body &b) const{
	uint32_t bx = x(b.x());
	uint32_t by = y(b.y());

	for(int i = 31; i >= 0; i--){
		os<<((by >> i) % 2);
		os<<((bx >> i) % 2);
	}
	os<<std::endl;

	return os;
}

bool MortonKeyCalculator::operator() (const Body &a, const Body &b) const{
	//Return true is a is less than b in a Morton Key ordering; 
        //otherwise, return false.
	uint32_t ax = x(a.x());
        uint32_t ay = y(a.y());
        uint32_t bx = x(b.x());
        uint32_t by = y(b.y());

        for (int i = 31; i >= 0; i--)
        {
                // Keep in mind the y coordinate is more significant, so
                // its evaluations should go first
                if ( ((ay >> i) % 2) < ((by >> i) % 2) ) return true;
                if ( ((ay >> i) % 2) > ((by >> i) % 2) ) return false;

                if ( ((ax >> i) % 2) < ((bx >> i) % 2) ) return true;
                if ( ((ax >> i) % 2) > ((bx >> i) % 2) ) return false;
        }

        return true;
}
