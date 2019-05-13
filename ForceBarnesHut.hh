#ifndef _FORCEBARNESHUT_HH_
#define _FORCEBARNESHUT_HH_

#include "ForceCalculator.hh"
#include "MortonKeyCalculator.hh"
#include "Body.hh"

// Leverage Body class because Quadtree contains Bodies
class Quadtree : public Body
{
	public:
		Quadtree(double x, double y, double d, double m, Quadtree **qt_bodies, Body **bodies);
		~Quadtree();
		void updateForce(Body* body, double d); //update gravitational force

	private:
		Body *body;
		Quadtree **qt_bodies_;
		Body **bodies_;
		double d;
};

// helper class for ForceBarnesHut::recalculate
class CompBody {
	public:
		CompBody(int level, MortonKeyCalculator *mkc); //comparator to determine bounds
		bool operator() (Body &body, int n);

	private:
		int level;
		MortonKeyCalculator* mkc;
};

class ForceBarnesHut : public ForceCalculator {
	public:
		ForceBarnesHut(Body *body, int N, double theta);

		virtual void operator() (Body *pulled);

		virtual ~ForceBarnesHut();

		Quadtree* recalculate(Body*, Body*, uint32_t, double); //recalculate center of mass

	private:
		Body *bodies_;
		const int N_;
		double theta_;
		MortonKeyCalculator mkc;
		Quadtree* tree;
};

#endif
