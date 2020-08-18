#include <fstream>
#include <iostream>
#include <vector>
#include <stdlib.h>	
#include <time.h>						
#include <math.h>							
#include <chrono>
#include <random>								
#include <algorithm>	
#include <deque>
#include <iomanip>
#include <cmath>
#include <string>
#include <cstring>	
#include <cassert>		
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <boost/geometry/algorithms/is_convex.hpp>

using namespace std;
namespace boost_geo = boost::geometry;
namespace trans = boost::geometry::strategy::transform;

typedef boost_geo::model::point<long double, 2, boost_geo::cs::cartesian> Point;
typedef boost::geometry::model::linestring<Point> Line;
typedef boost_geo::model::polygon<Point, false> Polygon;
typedef boost_geo::model::multi_polygon<Polygon> MultiPolygon;
typedef boost_geo::model::box<Point> Box;
typedef Point Vector;

const long double EPS = 1e-8;
const long double INF = 4e18;
const long double PI = acos(-1);
const long double RUN_TIME_DURATION = 20;
const long double FEASIBILTY = 1e-4;
const int NUMBER_OF_HOST_NESTS = 15;
const long double NFP_TRANSLATOR_COF = 1e4;
const long double INITIAL_STOCK_LENGTH = 100000000;
const int MAXIMUM_ITERATIONS_FOR_LOCAL_MINIMA = 50;
const long double LENGTH_INCREASING_RATE = 0.01;
const long double LENGTH_DECREASING_RATE = 0.04;
const vector<long double> ALLOWABLE_ROTATIONS = {0, 90, 180, 270};

static int frameno;

# define M_PI				3.14159265358979323846
# define N_ITERACIONES		30


namespace cs
{
	long double testFuncion(vector<long double>& x);
	long double function3(vector<long double>& x);
	long double function3(vector<long double>& x);
	long double Schwefel(vector<long double>& x);
	long double Ackley(vector<long double>& x, long double a = 20.0, long double b = 0.2, long double c = 2 * M_PI);
	void saveData(vector<long double>& bestValues);
	void initializeParams(int& n, long double& Pa);
	long double fRand(long double fMin, long double fMax);
	void normalDistribution(vector<long double>& normalDistr, long double& sigma, bool multiplicate);
	void mantegna(vector<long double>& step, vector<long double>& u, vector<long double>& v, long double& beta);
	vector<long double> simpleBounds(vector<long double>& s, vector<long double>& lowerBound, vector<long double>& upperBound);
	vector<vector<long double>> getCuckoos(vector<vector<long double>>& nest, vector<long double>& bestNest, vector<long double>& lowerBound, vector<long double>& upperBound);
	long double getBestNest(vector<long double>& bestNest, vector<vector<long double>>& nest, vector<vector<long double>>& newNest, vector<vector<long double>>& fitness, 
		MultiPolygon packing, vector<vector<long double>> &penalty, int polygon_id, long double rotationAngle, long double width, long double length);
	vector<vector<long double>> randPerm(vector<vector<long double>> nest);
	vector<vector<long double>> emptyNests(vector<vector<long double>>& nest, vector<long double>& lowerBound, vector<long double>& upperBound, long double& Pa);
	Point cuckooSearch(MultiPolygon &packing, vector<vector<long double>> &penalty, int polygon_id, long double rotationAngle, long double width, long double length);
};

namespace cp
{
	int dblcmp(long double d, long double eps);
	long double getLength(MultiPolygon &multiPolygon);
	Polygon translate(Polygon &polygon, Point translationPoint);
	Polygon rotateCW(Polygon &polygon, long double rotationAngleInDegree, Point referencePoint);
	long double randomInRange(long double a, long double b);
	long double getPackingDensity(MultiPolygon &packing);
	long double getPenetrationDepth(Polygon &polygonA, Polygon &polygonB);
	long double getTotalPenetrationDepth(MultiPolygon &packing);
	Polygon getInnerFitRectangle(MultiPolygon cluster, long double length, long double width);
	void increasePenalty(MultiPolygon &packing, vector<vector<long double>> &penalty);
	void pushDown(MultiPolygon &packing, long double length);
	void visualize(MultiPolygon multiPolygon, string outputLocation, string datasetName);
	long double polygonPolygonIntersectionArea(Polygon &polygon1, Polygon &polygon2);
	bool isFeasible(MultiPolygon &packing, long double totalAreaOfInputPolygons);
	long double getOverlapPenalty(
	MultiPolygon packing, vector<vector<long double>> &penalty, int polygon_id,
	long double rotationAngle, Point translationPoint);
	void readWKTPolygon(Polygon &polygon, string filename);
	void readWKTMultiPolygon(MultiPolygon &multiPolygon, string filename);
	MultiPolygon minimizeOverlap( MultiPolygon packing, vector<long double> allowableRoatations, long double width, long double length, long double totalAreaOfInputPolygons);
	void cuckooPacking(long double width, long double runTimeDuration);
};



int main()
{
	cp::cuckooPacking(100, 30); // dighe1
	return 0;
}
namespace cs 
{

	/* ============================================ Objective Functions ================================================ */

	/* A d-dimensional Test function 
	 * @param x: a d-dimensional vector 
	 * @return the total sum of each (x_i - 1)^2 */
	long double testFuncion(vector<long double>& x)
	{
		long double totalSum = 0.0;
		for (size_t i = 0; i < x.size(); ++i) {
			totalSum += pow(x[i] - 1.0, 2);
		}

		return totalSum;
	}

	/* This function is characterized by forming waves of local maximums and minimums 
	 * @param x: a d-dimensional vector 
	 * @return the function3 evaluated in this x vector */
	long double function3(vector<long double>& x){
		long double summation = 0.0;

		for (size_t i = 0; i < x.size(); ++i) {
			summation += pow(x[i], 2.0);
		}

		long double numerator = pow(sin(sqrt(summation)), 2.0) - 0.5;
		long double denominator = pow(1.0 + 0.001 * summation, 2.0);

		return 0.5 - numerator / denominator;
	}

	/* This function is very complex with many local minimums
	* @param x: a d-dimensional vector
	* @return the Schwefel evaluated in this x vector */
	long double Schwefel(vector<long double>& x){
		int d = x.size();

		long double sumatoria1 = 0.0;
		for (size_t i = 0; i < x.size(); ++i) {
			sumatoria1 += x[i] * sin(sqrt(abs(x[i])));
		}

		return 418.9829 * d - sumatoria1;
	}

	/* This function is characterized by an almost flat region and a large hole in the center.
	* @param x: a d-dimensional vector
	* @return the Ackley evaluated in this x vector */
	long double Ackley(vector<long double>& x, long double a, long double b, long double c){
		int d = x.size();

		// First Summation
		long double firstSummation = 0.0;
		for (size_t i = 0; i < x.size(); ++i) {
			firstSummation += pow(x[i], 2);
		}
		firstSummation = sqrt(firstSummation / d);

		// Second Summation
		long double secondSummation = 0.0;
		for (size_t i = 0; i < x.size(); ++i) {
			secondSummation += cos(c * x[i]);
		}
		secondSummation /= d;

		long double totalValue = (-1 * a * exp(-1 * b * firstSummation)) - exp(secondSummation) + a + exp(1);
		return totalValue;
	}


	/*================================================ Util Functions ==================================================== */

	/* Save in a *.txt format the best values */
	void saveData(vector<long double>& bestValues) 
	{
		ofstream myFile;
		myFile.open("data.txt");

		if (myFile.fail()) {
			cerr << "The file was not found" << endl;
			return;
		}

		for (size_t i = 0; i < bestValues.size(); i++) {
			myFile << bestValues[i];
			myFile << "\n";
		}
	}

	/* Params are initialized as indicated in the paper */
	void initializeParams(int& n, long double& Pa)
	{
		n = 25;
		Pa = 0.25;
	}

	/* Only generate random values between two numbers */
	long double fRand(long double fMin, long double fMax)
	{
		long double f = (long double)rand() / RAND_MAX;
		return fMin + f * (fMax - fMin);
	}

	/* Function to calculate normal distribution */
	void normalDistribution(vector<long double>& normalDistr, long double& sigma, bool multiplicate)
	{
		// construct a trivial random generator engine from a time-based seed:
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		default_random_engine generator(seed);
		normal_distribution<long double> distribution(0.0, 1.0);

		if (multiplicate == true) {
	#pragma omp parallel for		
			for (size_t i = 0; i < normalDistr.size(); ++i) {
				normalDistr[i] = distribution(generator) * sigma;
			}
		}
		else {
	#pragma omp parallel for		
			for (size_t i = 0; i < normalDistr.size(); ++i) {
				normalDistr[i] = distribution(generator);
			}
		}

	}

	/* Auxiliar algorithm for Mantegna */
	void mantegna(vector<long double>& step, vector<long double>& u, vector<long double>& v, long double& beta)
	{
		for (size_t i = 0; i < step.size(); ++i) {
			step[i] = u[i] / pow(abs(v[i]), 1 / beta);
		}
	}


	/* Application of simple constraints */
	vector<long double> simpleBounds(vector<long double>& s, vector<long double>& lowerBound, vector<long double>& upperBound)
	{
		// Apply the lower bound
		vector<long double> nsTemp = s;
	#pragma omp parallel for	
		for (size_t i = 0; i < s.size(); ++i) {
			if (nsTemp[i] < lowerBound[i]) nsTemp[i] = lowerBound[i];
		}

		// Apply the upper bound
	#pragma omp parallel for	
		for (size_t i = 0; i < s.size(); ++i) {
			if (nsTemp[i] > upperBound[i]) nsTemp[i] = upperBound[i];
		}

		// Update this new move 
		return nsTemp;
	}

	/* Get cuckoos by ramdom walk - Levy implementation */
	vector<vector<long double>> getCuckoos(vector<vector<long double>>& nest, vector<long double>& bestNest, vector<long double>& lowerBound, vector<long double>& upperBound)
	{
		int n = nest.size();
		vector<vector<long double>> newNest(n, vector<long double>(nest[0].size()));

		// Levy exponent and coefficient
		long double beta = 1.5;
		long double sigma = (tgamma(1 + beta)* sin(M_PI * beta / 2) / (tgamma((1 + beta) / 2) * beta * pow(2.0, (beta - 1) / 2)));
		long double sigmaPow = pow(sigma, 1 / beta);

		for (int i = 0; i < n; ++i) {
			vector<long double> s = nest[i];

			// Levy flights by Mantegna's algorithm
			vector<long double> u(s.size()), v(s.size()), step(s.size());
			normalDistribution(u, sigmaPow, true);
			normalDistribution(v, sigmaPow, false);
			mantegna(step, u, v, beta);													// fill step variable

			vector<long double> stepsize(s.size());
			for (size_t j = 0; j < stepsize.size(); ++j) {
				stepsize[j] = 0.01 * step[j] * (s[j] - bestNest[j]);
			}

			// Actual random walks or flights
			vector<long double> auxNormalD(s.size());
			normalDistribution(auxNormalD, sigmaPow, false);

			for (size_t j = 0; j < s.size(); ++j) {
				s[j] = s[j] + stepsize[j] * auxNormalD[j];
			}

			// Apply simple bounds/limits
			nest[i] = simpleBounds(s, lowerBound, upperBound);			// newNest[i] = simpleBounds(s, lowerBound, upperBound);	
		}

		newNest = nest;

		return newNest;
	}

	/* Find the current best nest and fill bestNest param */
	long double getBestNest(vector<long double>& bestNest, vector<vector<long double>>& nest, vector<vector<long double>>& newNest, vector<vector<long double>>& fitness, 
		MultiPolygon packing, vector<vector<long double>> &penalty, int id, long double rotationAngle, long double width, long double length)
	{
		// Evaluating all new solutions
		for (size_t i = 0; i < nest.size(); ++i) {
			Point translationPoint(newNest[i][0], newNest[i][1]);


			Polygon poly = packing[id];
			packing[id] = cp::rotateCW(packing[id], rotationAngle, packing[id].outer()[0]);
			packing[id] = cp::translate(packing[id], Point(-packing[id].outer()[0].get<0>(), -packing[id].outer()[0].get<1>()));
			packing[id] = cp::translate(packing[id], translationPoint);

			long double fNew = cp::getTotalPenetrationDepth(packing);		
		//	std::cout << " fitness@ " << fNew << " : nest@ " << newNest[i][0] << " " << newNest[i][1] << "\n";
			if (fNew <= fitness[i][0]) {
				fitness[i][0] = fNew;
				nest[i] = newNest[i];
			}

			packing[id]=poly;
		}

		// Find the current best
		long double fitnessMin = fitness[0][0];
		int k = 0;
	#pragma omp parallel	
		for (size_t i = 1; i < fitness.size(); ++i) {
			if (fitness[i][0] < fitnessMin) {
				fitnessMin = fitness[i][0];
				k = i;
			}
		}

		bestNest = nest[k];

		// std::cout << "best value: " << fitnessMin << "\n"; 
		// std::cout << "best nest: " << bestNest[0] << " " << bestNest[1] << "\n"; 

		return fitnessMin;
	}

	/* rearrange the matrix from rows */
	vector<vector<long double>> randPerm(vector<vector<long double>> nest)
	{
		int n = nest.size();
		vector<vector<long double>> nestCopy = nest;
		vector<int> index(n);
		for (int i = 0; i < n; ++i) {
			index[i] = i;
		}

		// Random permutation the order
		for (int i = 0; i < n; ++i) {
			int j = rand() % (n - i) + i;
			// Swap i and j
			int t = index[j];
			index[j] = index[i];
			index[i] = t;
		}
	#pragma omp parallel
		for (int i = 0; i < n; i++) {
			nestCopy[i] = nest[index[i]];
		}

		return nestCopy;
	}

	/* Replace some nests by constructing new solutions/nests */
	vector<vector<long double>> emptyNests(vector<vector<long double>>& nest, vector<long double>& lowerBound, vector<long double>& upperBound, long double& Pa)
	{
		int n = nest.size();

		// Discovered or not -- a status vector (K contains 0 or 1)
		unsigned seed2 = std::chrono::system_clock::now().time_since_epoch().count();
		srand(seed2);
		vector<vector<int>> K(n, vector<int>(nest[0].size()));
		for (int i = 0; i < n; ++i) {
			for (size_t j = 0; j < nest[0].size(); ++j) {
				if (fRand(-1.0, 1.0) > Pa) K[i][j] = 1.0;
				else  K[i][j] = 0.0;
			}
		}

		// New solution by biased/selective random walks
		unsigned seed3 = std::chrono::system_clock::now().time_since_epoch().count();
		srand(seed3);
		vector<vector<long double>> nestRand1 = randPerm(nest);
		vector<vector<long double>> nestRand2 = randPerm(nest);
		vector<vector<long double>> stepsize(n, vector<long double>(nest[0].size()));

		for (int i = 0; i < n; ++i) {
			for (size_t j = 0; j < nest[0].size(); ++j) {
				stepsize[i][j] = fRand(0, 1.0) * (nestRand1[i][j] - nestRand2[i][j]);
			}
		}

		vector<vector<long double>> newNest(n, vector<long double>(nest[0].size()));
		// Fill newNest
		for (int i = 0; i < n; ++i) {
			for (size_t j = 0; j < nest[0].size(); ++j) {
				newNest[i][j] = nest[i][j] + stepsize[i][j] * K[i][j];
			}

		}

	#pragma omp parallel
		for (int i = 0; i < n; ++i) {
			vector<long double> s = newNest[i];
			newNest[i] = simpleBounds(s, lowerBound, upperBound);
		}

		return newNest;
	}

	/* Principal function to Cuckoo Search */
	Point cuckooSearch(MultiPolygon &packing, vector<vector<long double>> &penalty, int polygon_id, long double rotationAngle, long double width, long double length)
	{
		int n = 15;
		long double Pa = 0.25;
		Polygon polygon = packing[polygon_id];
		polygon = cp::rotateCW(polygon, rotationAngle, polygon.outer()[0]);
		Polygon innerFitRectangle = cp::getInnerFitRectangle({polygon}, length, width);

		long double max_x = -INF, min_x = INF, max_y = -INF, min_y = INF;
		for (Point &point : innerFitRectangle.outer())
		{
			max_x = std::max(max_x, point.get<0>());
			min_x = std::min(min_x, point.get<0>());
			max_y = std::max(max_y, point.get<1>());
			min_y = std::min(min_y, point.get<1>());
		}

		long double tolerance = 1.0e-5;
		int numberDomain = 2;											// x = (x_1, x_2, ... , x_n)
		int nIterations = 0;											// number of iterations

		// Bounds
		vector<long double> lowerBound(numberDomain, -5.0);					// limits of domain
		vector<long double> upperBound(numberDomain, 5.0);					// limits of domain

		lowerBound[0] = min_x;
		lowerBound[1] = min_y;
		upperBound[0] = max_x;
		upperBound[1] = max_y;

		// Random initial solutions
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		srand(seed);

		// std::cout << "Min_x: " << min_x << " Max_x: " << max_x << "\n";
		// std::cout << "Min_y: " << min_x << " Max_y: " << max_y << "\n";
		vector<vector<long double>> nest(n, vector<long double>(numberDomain));
		for (int i = 0; i < n; ++i) {
			nest[i][0] = min_x + (max_x - min_x) * fRand(0, 1.0);
			nest[i][1] = min_y + (max_y - min_y) * fRand(0, 1.0);
			// std::cout << "Nest[i] = " << nest[i][0] << " " << nest[i][1] << "\n";
		}
		// getchar();

		// Get the current best
		vector<vector<long double>> fitness(n, vector<long double>(1));
		for (size_t i = 0; i < fitness.size(); ++i) {
			fitness[i][0] = INF;
		}

		vector<long double> bestNest;
		long double fMin = getBestNest(bestNest, nest, nest, fitness, packing, penalty, polygon_id, rotationAngle, length, width);

		//getchar();
		cout << "Current minima@" << fMin << " ";
		// To save data
		vector<long double> bestValues(N_ITERACIONES);

		while (fMin > tolerance && nIterations < N_ITERACIONES) {

			// Generate new solutions (but keep the current best)
			vector<vector<long double>> newNest = getCuckoos(nest, bestNest, lowerBound, upperBound);
			vector<long double> best;
			long double fNew = getBestNest(best, nest, newNest, fitness, packing, penalty, polygon_id, rotationAngle, length, width);

			// Discovery and randomization
			newNest = emptyNests(nest, lowerBound, upperBound, Pa);

			// Evaluate this set of solutions
			fNew = getBestNest(best, nest, newNest, fitness, packing, penalty, polygon_id, rotationAngle, length, width);
			//cout << "fmin N " << nIterations << ": " << fNew << endl;

			// Find the best objective so far
			if (fNew < fMin) {
				fMin = fNew;
				bestNest = best;
			}
			
			// Update the counter
			nIterations++;
		}
		cout << "best minima@" << fMin << "\n";
		return Point(bestNest[0], bestNest[1]);
	}
};
namespace cp
{
	int dblcmp(long double d, long double eps)
	{
		if (std::fabs(d) < eps)
		{
			return 0;
		}
		return d > eps ? 1 : -1;
	}

	long double getLength(MultiPolygon &multiPolygon)
	{
		long double min_y = INF, max_y = -INF;
		for (auto polygon : multiPolygon)
		{
			for (auto point : polygon.outer())
			{
				min_y = std::min(point.get<1>(), min_y);
				max_y = std::max(point.get<1>(), max_y);
			}
		}
		return max_y - min_y;
	}

	Polygon translate(Polygon &polygon, Point translationPoint)
	{
		Polygon translatedPolygon;
		boost_geo::transform(polygon, translatedPolygon, trans::translate_transformer<long double, 2, 2>(translationPoint.get<0>(), translationPoint.get<1>()));
		return translatedPolygon;
	}

	Polygon rotateCW(Polygon &polygon, long double rotationAngleInDegree, Point referencePoint)
	{
		boost_geo::multiply_value(referencePoint, -1);
		Polygon translatedPolygon = translate(polygon, referencePoint), rotatedPolygon;
		boost_geo::transform(translatedPolygon, rotatedPolygon, trans::rotate_transformer<boost_geo::degree, long double, 2, 2>(rotationAngleInDegree));
		return rotatedPolygon;
	}
	long double randomInRange(long double a, long double b)
	{
		long double random = ((long double)std::rand()) / (long double)RAND_MAX;
		long double diff = b - a;
		long double r = random * diff;
		return a + r;
	}

	long double getPackingDensity(MultiPolygon &packing)
	{
		Box stock;
		boost_geo::envelope(packing, stock);
		long double stockArea = std::fabs(boost_geo::area(stock));
		long double polygonsArea = std::fabs(boost_geo::area(packing));
		long double packingDensity = polygonsArea / stockArea;
		return packingDensity * 100.0;
	}

	long double getPenetrationDepth(Polygon &polygonA, Polygon &polygonB)
	{
		MultiPolygon intersections;
		boost_geo::intersection(polygonA, polygonB, intersections);
		return std::fabs(boost_geo::area(intersections));
	}

	long double getTotalPenetrationDepth(MultiPolygon &packing)
	{
		long double totalPenetrationDepth = 0.0;
		int numberOfPolygons = packing.size();
		for (int i = 0; i < numberOfPolygons; i += 1)
		{
			for (int j = i + 1; j < numberOfPolygons; j += 1)
			{
				totalPenetrationDepth += getPenetrationDepth(packing[i], packing[j]);
			}
		}
		return totalPenetrationDepth;
	}
	Polygon getInnerFitRectangle(MultiPolygon cluster, long double length, long double width)
	{
		Polygon innerFitRectangle;
		Point reference = cluster[0].outer()[0];

		long double max_x = -INF, min_x = INF, max_y = -INF, min_y = INF;
		for (Polygon &polygon : cluster)
		{
			for (Point &point : polygon.outer())
			{
				point = Point(point.get<0>() - reference.get<0>(), point.get<1>() - reference.get<1>());
				max_x = std::max(max_x, point.get<0>());
				min_x = std::min(min_x, point.get<0>());
				max_y = std::max(max_y, point.get<1>());
				min_y = std::min(min_y, point.get<1>());
			}
		}
		innerFitRectangle.outer().push_back(Point(std::abs(min_x), std::abs(min_y)));
		innerFitRectangle.outer().push_back(Point(std::abs(min_x), length - max_y));
		innerFitRectangle.outer().push_back(Point(width - max_x, length - max_y));
		innerFitRectangle.outer().push_back(Point(width - max_x, std::abs(min_y)));
		innerFitRectangle.outer().push_back(Point(std::abs(min_x), std::abs(min_y)));
		return innerFitRectangle;
	}

	void increasePenalty(MultiPolygon &packing, vector<vector<long double>> &penalty)
	{
		int n = penalty.size();
		long double maximumPenetrationDepth = -INF;
		vector<vector<long double>> penetrationDepths(n, vector<long double>(n, 0));
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				if (i == j)
				{
					continue;
				}
				penetrationDepths[i][j] = getPenetrationDepth(packing[i], packing[j]);
				if (i < j)
				{
					maximumPenetrationDepth = std::max(maximumPenetrationDepth, penetrationDepths[i][j]);
				}
			}
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				penalty[i][j] += (penetrationDepths[i][j] / maximumPenetrationDepth);
			}
		}
	}

	void pushDown(MultiPolygon &packing, long double length)
	{
		for (auto &polygon : packing)
		{
			long double pushDownBy = 0;
			for (auto &point : polygon.outer())
			{
				long double dif = point.get<1>() - length;
				if (dblcmp(dif, EPS) > 0)
				{
					pushDownBy = std::max(pushDownBy, dif);
				}
			}
			for (auto &point : polygon.outer())
			{
				point = Point(point.get<0>(), point.get<1>() - pushDownBy);
			}
			for (auto &inner : polygon.inners())
			{
				for (auto &point : inner)
				{
					point = Point(point.get<0>(), point.get<1>() - pushDownBy);
				}
			}
		}
	}

	void visualize(MultiPolygon multiPolygon, string outputLocation, string datasetName)
	{
		Box box;
		boost::geometry::envelope(multiPolygon, box);
		// std::cout << "make_envelope..............: " << boost::geometry::dsv(box) << std::endl;
		std::ostringstream name;
		name << "frame_" << std::setw(4) << std::setfill('0') << frameno++ << "_" << datasetName << ".svg";
		std::ofstream svg(outputLocation + "/" + name.str());
		boost_geo::svg_mapper<Point> mapper(svg, 400, 500);
		mapper.add(multiPolygon);
		mapper.map(multiPolygon, "fill-opacity:0.5;fill:rgb(169,169,169);stroke:rgb(169,169,169);stroke-width:1");
		mapper.map(box, "opacity:0.8;fill:none;stroke:rgb(0,0,0);stroke-width:1;stroke-linecap:round");
	}

	long double polygonPolygonIntersectionArea(Polygon &polygon1, Polygon &polygon2)
	{
		std::deque<Polygon> polygons;
		boost_geo::intersection(polygon1, polygon2, polygons);
		long double intersectedArea = 0.0;
		for (auto polygon : polygons)
		{
			intersectedArea += boost_geo::area(polygon);
		}
		return std::fabs(intersectedArea);
	}

	bool isFeasible(MultiPolygon &packing, long double totalAreaOfInputPolygons)
	{
		long double overlappingArea = 0.0;
		long double feasibilityRatio = 0.0;
		int n = packing.size();
		for (int i = 0; i < n; i += 1)
		{
			for (int j = i + 1; j < n; j += 1)
			{
				overlappingArea += polygonPolygonIntersectionArea(packing[i], packing[j]);
			}
		}
		feasibilityRatio = overlappingArea / totalAreaOfInputPolygons;
		return dblcmp(feasibilityRatio - FEASIBILTY, EPS) <= 0 ? true : false;
	}

	long double getOverlapPenalty(
	MultiPolygon packing, vector<vector<long double>> &penalty, int polygon_id,
	long double rotationAngle, Point translationPoint)
	{
		int numberOfPolygons = packing.size();

		// rotating
		packing[polygon_id] =
			rotateCW(packing[polygon_id], rotationAngle, packing[polygon_id].outer().front());
		Point referecePoint = packing[polygon_id].outer().front();
		// translating
		Point newPosition = translationPoint;
		boost_geo::subtract_point(newPosition, referecePoint);
		packing[polygon_id] = translate(packing[polygon_id], newPosition);

		long double overlapPenalty = 0;
		for (int j = 0; j < numberOfPolygons; j++)
		{
			if (j != polygon_id)
			{
				overlapPenalty +=
					(penalty[polygon_id][j] * getPenetrationDepth(packing[polygon_id], packing[j]));
			}
		}
		return overlapPenalty;
	}

	void readWKTPolygon(Polygon &polygon, string filename)
	{
		std::ifstream polygonWKTFile(filename);
		assert(polygonWKTFile.is_open());

		string wktStr;
		polygonWKTFile.seekg(0, std::ios::end);
		wktStr.reserve(polygonWKTFile.tellg());
		polygonWKTFile.seekg(0, std::ios::beg);
		wktStr.assign((std::istreambuf_iterator<char>(polygonWKTFile)), std::istreambuf_iterator<char>());
		wktStr.pop_back();
		boost_geo::read_wkt(wktStr, polygon);
		boost_geo::correct(polygon);
		polygonWKTFile.close();
	}

	void readWKTMultiPolygon(MultiPolygon &multiPolygon, string filename)
	{
		std::ifstream multiPolygonWKTFile(filename);
		assert(multiPolygonWKTFile.is_open());

		string wktStr;
		multiPolygonWKTFile.seekg(0, std::ios::end);
		wktStr.reserve(multiPolygonWKTFile.tellg());
		multiPolygonWKTFile.seekg(0, std::ios::beg);
		wktStr.assign((std::istreambuf_iterator<char>(multiPolygonWKTFile)), std::istreambuf_iterator<char>());
		wktStr.pop_back();
		boost_geo::read_wkt(wktStr, multiPolygon);
		boost_geo::correct(multiPolygon);

		multiPolygonWKTFile.close();
	}

	MultiPolygon minimizeOverlap( MultiPolygon packing, vector<long double> allowableRoatations, long double width, long double length, long double totalAreaOfInputPolygons)
    {
		int numberOfPolygons = packing.size();
		vector<vector<long double>> penalty(numberOfPolygons, vector<long double>(numberOfPolygons, 1.0));
		int it = 0;
		long double fitness = INF;
		vector<int> Q(numberOfPolygons);
		for (int i = 0; i < numberOfPolygons; i++)
		{
			Q[i] = i;
		}
		MultiPolygon bestPacking = packing;
		long double bestTotalOverlapPenalty = getTotalPenetrationDepth(bestPacking);
		int noProgressOfPenalyCount = 0;
		while (it < MAXIMUM_ITERATIONS_FOR_LOCAL_MINIMA)
		{
			std::random_shuffle(Q.begin(), Q.end());
			for (int i = 0; i < numberOfPolygons; i++)
			{
				long double overlapPenalty = getTotalPenetrationDepth(packing);

				Point bestLocation(INF, INF);
				long double bestRotationAngle = 0.0;
				std::cout << "iteration@" << it << " Polygon@" << i + 1 << "/" << numberOfPolygons << " penalty@before@" << overlapPenalty << " "; 
				for (long double rotationAngle : allowableRoatations)
				{
					Point translationPoint = cs::cuckooSearch(packing, penalty, Q[i], rotationAngle, width, length);

					int id = Q[i];
					Polygon poly = packing[id];
					packing[id] = cp::rotateCW(packing[id], rotationAngle, packing[id].outer()[0]);
					packing[id] = cp::translate(packing[id], Point(-packing[id].outer()[0].get<0>(), -packing[id].outer()[0].get<1>()));
					packing[id] = cp::translate(packing[id], translationPoint);

					long double currentOverlapPenalty = getTotalPenetrationDepth(packing);

					if (currentOverlapPenalty < overlapPenalty)
					{
						overlapPenalty = currentOverlapPenalty;
						bestLocation = translationPoint;
						bestRotationAngle = rotationAngle;
					}

					packing[id] = poly;
				}
				if (dblcmp(INF - bestLocation.get<0>(), EPS) <= 0 or
					dblcmp(INF - bestLocation.get<1>(), EPS) <= 0)
				{
					continue;
				}
				cout << "penalty@after@" << overlapPenalty << "\n";
				packing[Q[i]] = rotateCW(packing[Q[i]], bestRotationAngle, packing[Q[i]].outer().front());
				Point newPosition = bestLocation;
				boost_geo::subtract_point(newPosition, packing[Q[i]].outer().front());
				packing[Q[i]] = translate(packing[Q[i]], newPosition);
			}
			long double currentTotalOverlapPenalty = getTotalPenetrationDepth(packing);
			if (currentTotalOverlapPenalty < bestTotalOverlapPenalty)
			{
				noProgressOfPenalyCount = 0;
				bestTotalOverlapPenalty = currentTotalOverlapPenalty;
				bestPacking = packing;
				if (isFeasible(bestPacking, totalAreaOfInputPolygons))
				{
					return bestPacking;
				}
			}
			long double totalPenetrationDepth = getTotalPenetrationDepth(packing);
			if (dblcmp(totalPenetrationDepth, EPS) == 0)
			{
				return packing;
			}
			else if (totalPenetrationDepth < fitness)
			{
				fitness = totalPenetrationDepth;
				it = 0;
			}
			increasePenalty(packing, penalty); // increase penalty
			it += 1;
		}
		return bestPacking;
	}

	void cuckooPacking(long double width, long double runTimeDuration)
	{
		string initialSolutionFileName = "./initial_packing.wkt";
		MultiPolygon initialPacking;
		readWKTMultiPolygon(initialPacking, initialSolutionFileName);

		std::cout << "reading input initial packing...\n";
		std::cout << boost_geo::wkt(initialPacking) << "\n";
		std::cout << "reading successful!\n";

		long double totalAreaOfInitialPackingPolygons = 0;
		for (Polygon &polygon : initialPacking)
		{
			totalAreaOfInitialPackingPolygons += fabs(boost_geo::area(polygon));
		}

		MultiPolygon bestPacking = initialPacking;
		MultiPolygon currentPacking = initialPacking;

		long double bestLenght = getLength(initialPacking);
		long double decreasingRate = LENGTH_DECREASING_RATE;
		long double increasingRate = LENGTH_INCREASING_RATE;
		long double currentLength = bestLenght;
		int feasiblePackingId = 0;

		std::cout << "Initial packing density: " << getPackingDensity(initialPacking) << "\n";
		auto start = std::chrono::high_resolution_clock::now();
		while (true)
		{
			auto stop = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start);
			if (duration.count() >= runTimeDuration)
			{
				break;
			}
			visualize(currentPacking, "./dighe1/", "minimized");
			if (isFeasible(currentPacking, totalAreaOfInitialPackingPolygons))
			{
				long double currentAccuracy = getPackingDensity(currentPacking);
				bestPacking = currentPacking;
				std::cout << "Current Accuray@"  << frameno << "@" << currentAccuracy << "\n";
				bestLenght = currentLength;
				currentLength = (1.0 - decreasingRate) * currentLength;
				pushDown(currentPacking, currentLength);
			}
			else
			{
				currentLength = (1.0 + increasingRate) * currentLength;
			}
			currentPacking = minimizeOverlap(currentPacking, ALLOWABLE_ROTATIONS, width, currentLength, totalAreaOfInitialPackingPolygons);
		}
	}
};