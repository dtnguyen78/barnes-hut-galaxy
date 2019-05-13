#include "ForceBarnesHut.hh"
#include "Timer.hh"

#ifdef _OPENMP
#include <parallel/algorithm>
#include <parallel/numeric>
#else
#include <algorithm>
#endif

using namespace std;

/* Taken from MortonKeyCalculator.cc */
// It's bad practice to include a .cc file

bool compX(const Body &a, const Body &b) { return a.x() < b.x(); }
bool compY(const Body &a, const Body &b) { return a.y() < b.y(); }

MortonKeyCalculator::MortonKeyCalculator(Body *bodies, int N){
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

/* Definitions for Quadtree data structure */

Quadtree::Quadtree(double x, double y, double d, double m, Quadtree **qt_bodies, Body **bodies) :
	Body(x, y, 0.0, 0.0, m), // body isn't moving yet. initialize (x,y) coords and mass
	qt_bodies_(qt_bodies),
	bodies_(bodies),
	d(d) {};

Quadtree::~Quadtree(){
	for(int i = 0; i < 4; i++)
		delete qt_bodies_[i];
	delete[] qt_bodies_;
		//qt_bodies_[i] = NULL;
	//qt_bodies_ = NULL;
};

// used to determine bounds when creating Quadtrees and recalculating center of mass
CompBody::CompBody(int level, MortonKeyCalculator *mkc) :
	level(level),
	mkc(mkc) {};

bool CompBody::operator()(Body &body, int n)
{
	int x = ((this->mkc->x(body.x()) & this->level) > 0) * 1;
	int y = ((this->mkc->y(body.y()) & this->level) > 0) * 2;
	bool result = (x + y) < n;

	return result;
}

void Quadtree::updateForce(Body *pulled, double theta) {
	if (pulled == NULL) return;

	// s = width of region represented by node
	// d = distance between body and node's center of mass;
	//	sqrt((x1-x2)^2 + (y1-y2)^2)
	double xx = pulled->x() - this->x();
	double yy = pulled->y() - this->y();
//	double dist = sqrt(xx*xx + yy*yy);
	xx *= xx;
	yy *= yy;
	double dist = sqrt(xx+yy);

	// If s/d < theta, treat node as a single body, and calculate
	// the force it exerts on body b, and add this to b's net
	// force.
	if (this->d / dist < theta)
		pulled->accGravityFrom(*this); //*(this->body)
	else
	{
		for (int i = 0; i < 4; i++)
		{
			if (qt_bodies_[i])
				qt_bodies_[i]->updateForce(pulled, theta);
			else if (bodies_[i])
				pulled->accGravityFrom(*bodies_[i]);
		}
	}
};

/* Start implementing ForceBarnesHut algorithm */

Quadtree* ForceBarnesHut::recalculate(Body *left, Body *right, uint32_t level, double d)
{
	auto comp = CompBody(level, &mkc);

	// Define bounds to arrange bodies to update (x,y), mass, and force later
	Body **bounds = new Body*[5];
	bounds[0] = left;
	# pragma omp parallel for
	for(int i = 1; i < 4; i++)
		bounds[i] = lower_bound(left, right, i, comp);
	bounds[4] = right;

	// always 4 children
	Quadtree **qt_bodies = new Quadtree*[4]; // internal nodes
	Body **bodies = new Body*[4]; // external nodes

	level = level >> 1;
	d = d / 2.0;

	// for (int i = 0; i < 4; i++)
	//	if (bounds[i+1] - bounds[0])
	//		qt_bodies[i] = recalculate(bounds[i], bounds[i+1], level, d);

	// This was 3x slower. :(
	// Use loop unrolling instead
	qt_bodies[0] = NULL;
	qt_bodies[1] = NULL;
	qt_bodies[2] = NULL;
	qt_bodies[3] = NULL;

	if(bounds[1] - bounds[0] > 1)
		qt_bodies[0] = recalculate(bounds[0], bounds[1], level, d);
	if(bounds[2] - bounds[1] > 1)
		qt_bodies[1] = recalculate(bounds[1], bounds[2], level, d);
	if(bounds[3] - bounds[2] > 1)
		qt_bodies[2] = recalculate(bounds[2], bounds[3], level, d);
	if(bounds[4] - bounds[3] > 1)
		qt_bodies[3] = recalculate(bounds[3], bounds[4], level, d);

	for(int i = 0; i < 4; i++)
	{
		if (bounds[i+1] - bounds[i] == 1)
			bodies[i] = bounds[i];
		else
			bodies[i] = NULL;
	}

	// Update center of mass (x,y) and create new body
	double x = 0, y = 0, m = 0;
//	#pragma omp parallel for default(shared) reduction(+:x,y,m)
	for (int i = 0; i < 4; i++)
	{
		if (qt_bodies[i])
		{
			Quadtree* &body = qt_bodies[i]; // remember to localize for speed
			x += body->x() * body->m();
			y += body->y() * body->m();
			m += body->m(); // Update total mass
		}
		else if (bodies[i])
		{
			Body* &body = bodies[i]; // remember to localize for speed
			x += body->x() * body->m();
			y += body->y() * body->m();
			m += body->m(); // Update total mass	
		}
	}
	x = x/m; // x = (x1*m1 + ... + xn*mn) / m
	y = y/m; // y = (y1*m1 + ... + yn*yn) / m

	return new Quadtree(x, y, d, m, qt_bodies, bodies);
}

ForceBarnesHut::ForceBarnesHut(Body *bodies, int N, double theta):
	bodies_(bodies),
	N_(N),
	theta_(theta),
	mkc(bodies_, N)
{

	//#ifdef _OPENMP
		__gnu_parallel::sort(bodies, bodies+N, mkc);
	//#endif

	//cache?
	auto level = 1u << 31u;
	auto d = mkc.cellWidth(1);

	#pragma omp parallel
	#pragma omp single nowait
	tree = recalculate(bodies, bodies + N, level, d);
};

ForceBarnesHut::~ForceBarnesHut()
{
	//tree = NULL;
	delete tree;
}

void ForceBarnesHut::operator()(Body *pulled){
	#pragma omp parallel
	#pragma omp single
	tree->updateForce(pulled, theta_);
};
