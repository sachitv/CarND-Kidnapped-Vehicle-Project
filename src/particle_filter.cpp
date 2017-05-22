/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <limits>
#include <map>

#include "particle_filter.h"

static int const NUM_PARTICLES = 500;
static double const INIT_WEIGHT = 1.f;

double NormalizeAngleRadians( double angleInRadians )
{
	return remainder( angleInRadians, 2.0 * M_PI );
}

void ParticleFilter::init( double x, double y, double theta, double std[] )
{
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	using namespace std;

	//Set the number of particles
	num_particles = NUM_PARTICLES;

	//Generate Normal Distribution aka Gaussian Distributions
	double const &std_x = std[ 0 ];
	normal_distribution<double> dist_x( x, std_x );

	double const &std_y = std[ 1 ];
	normal_distribution<double> dist_y( y, std_y );

	double const &std_theta = std[ 2 ];
	//Normalize theta before using it
	theta = NormalizeAngleRadians( theta );
	normal_distribution<double> dist_theta( theta, std_theta );

	//Initialize particles and weights
	default_random_engine gen;
	particles.reserve( NUM_PARTICLES );

	for ( size_t i = 0; i < NUM_PARTICLES; ++i )
	{
		Particle p;
		p.id = (int) i;
		p.weight = INIT_WEIGHT;
		p.x = dist_x( gen );
		p.y = dist_y( gen );
		p.theta = NormalizeAngleRadians( dist_theta( gen ));

		particles.emplace_back( p );
		weights.emplace_back( INIT_WEIGHT );

		//cout<<"Particle "<< p.id <<" initialized "<<p.x<<", "<<p.y<<endl;
	}
	is_initialized = true;
}

void ParticleFilter::prediction( double delta_t, double std_pos[], double velocity, double yaw_rate )
{
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	//Generate Normal Distribution aka Gaussian Distributions
	using namespace std;
	double const &std_x = std_pos[ 0 ];
	normal_distribution<double> dist_x( 0, std_x );

	double const &std_y = std_pos[ 1 ];
	normal_distribution<double> dist_y( 0, std_y );

	double const &std_theta = std_pos[ 2 ];
	normal_distribution<double> dist_theta( 0, std_theta );

	default_random_engine gen;

	//Compute the prediction for each of the particles and add the noise sampled from the above distribution
	for ( auto &particle : particles )
	{
		double const yawdelta = yaw_rate * delta_t;
		double const angle_delta = NormalizeAngleRadians( particle.theta + (yawdelta));
		double const vel_by_yaw_rate = (velocity / yaw_rate);

		particle.x = particle.x + (vel_by_yaw_rate * (sin( angle_delta ) - sin( particle.theta )));
		particle.y = particle.y + (vel_by_yaw_rate * (cos( particle.theta ) - cos( angle_delta )));
		particle.theta = particle.theta + yawdelta;

		particle.x += dist_x( gen );
		particle.y += dist_y( gen );
		particle.theta += NormalizeAngleRadians( dist_theta( gen ));
		NormalizeAngleRadians( particle.theta );

		//cout<<"Particle "<< particle.id <<" prediction "<<particle.x<<", "<<particle.y<<endl;
	}
}

//This function is not necessary
void ParticleFilter::dataAssociation( std::vector<LandmarkObs> &predicted, std::vector<LandmarkObs> &observations ) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

	//UNUSED PARAMETERS
	(void)predicted;
	(void)observations;

}

void ParticleFilter::updateWeights( double sensor_range, double std_landmark[],
									std::vector<LandmarkObs> observations, Map map_landmarks )
{
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html

	//For each particle
	using namespace std;

	int weightcounter = 0;
	for ( auto &particle : particles )
	{
		double resultweight = 1.f;
		particle.weight = 1.f;

		map<Map::single_landmark_s const*, int> landmarkToDist;
		map<Map::single_landmark_s const*, LandmarkObs> landmarkToObs;

		//Check each observation
		for ( auto const &observation : observations )
		{
			//Transform this observation into the space of the map
			LandmarkObs txObs;

			txObs.id = observation.id;
			txObs.x = (observation.x * cos( particle.theta )) - (observation.y * sin( particle.theta )) + particle.x;
			txObs.y = (observation.x * sin( particle.theta )) + (observation.y * cos( particle.theta )) + particle.y;

			//Look through the map and see which ones this could possibly be?
			double dist2min = std::numeric_limits<double>::max();
			Map::single_landmark_s const *bestLandmark = &(map_landmarks.landmark_list[ 0 ]);

			for ( auto const &landmark : map_landmarks.landmark_list )
			{
				double const dist2landmark = dist2( txObs.x, txObs.y, landmark.x_f, landmark.y_f );
				if ( dist2landmark < dist2min )
				{
					dist2min = dist2landmark;
					bestLandmark = &landmark;
				}
			}

			if(landmarkToDist.find(bestLandmark) == landmarkToDist.end())
			{
				//Not Found, insert
				landmarkToDist[bestLandmark] = dist2min;
				landmarkToObs[bestLandmark] = txObs;
			}
			else
			{
				//Insert
				if(landmarkToDist[bestLandmark] > dist2min)
				{
					landmarkToDist[bestLandmark] = dist2min;
					landmarkToObs[bestLandmark] = txObs;
				}
			}
		}

		for(auto const& iter : landmarkToObs)
		{
			Map::single_landmark_s const * const bestLandmark = iter.first;
			LandmarkObs const& txObs = iter.second;

			//Calculate multivariate gaussian prob
			double const p1 = 1.f / (2.f * M_PI * std_landmark[ 0 ] * std_landmark[ 1 ]);
			double const p2 = -(
					(pow((txObs.x - bestLandmark->x_f), 2.f ) / (2.f * std_landmark[ 0 ] * std_landmark[ 0 ])) +
					(pow((txObs.y - bestLandmark->y_f), 2.f ) / (2.f * std_landmark[ 1 ] * std_landmark[ 1 ])));
			double const prob = p1 * exp( p2 );

			resultweight *= prob;
		}

		particle.weight = resultweight;

		//cout<<"Particle "<< particle.id <<" final weight "<<particle.weight<<endl;

		weights[ weightcounter ] = resultweight;
		weightcounter++;
	}

}

void ParticleFilter::resample()
{
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	using namespace std;
	discrete_distribution<> d( weights.begin(), weights.end());
	default_random_engine gen;

	std::vector<Particle> newParticles;
	newParticles.reserve( NUM_PARTICLES );

	for ( auto i = 0; i < NUM_PARTICLES; ++i )
	{
		newParticles.emplace_back( particles[ d( gen ) ] );
	}

	swap( particles, newParticles );
}

void ParticleFilter::write( std::string filename )
{
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open( filename, std::ios::app );
	for ( int i = 0; i < NUM_PARTICLES; ++i )
	{
		dataFile << particles[ i ].x << " " << particles[ i ].y << " " << particles[ i ].theta << "\n";
	}
	dataFile.close();
}
