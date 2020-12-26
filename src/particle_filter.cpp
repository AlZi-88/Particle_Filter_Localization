/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1.
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method
   *   (and others in this file).
   */
  num_particles = 1000;  // TODO: Set the number of particles
  particles = vector<Particle>(num_particles);
  std::default_random_engine gen;
  std::normal_distribution<double> dist_x(x,std[0]);
  std::normal_distribution<double> dist_y(y,std[1]);
  std::normal_distribution<double> dist_theta(theta,std[2]);

  for (int i = 0; i < num_particles; i++){
    particles[i].id = i;
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);
    particles[i].weight = 0.0;

  }

  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[],
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
   double theta_old, new_theta, new_x, new_y;

   std::default_random_engine gen;

   if (yaw_rate != 0){
     for (int i = 0; i<num_particles; i++){
       theta_old = particles[i].theta;
       new_theta = theta_old + yaw_rate*delta_t;
       new_x = particles[i].x + velocity*(std::sin(new_theta) - std::sin(theta_old))/yaw_rate;
       new_y = particles[i].y + velocity*(std::cos(theta_old) - std::cos(new_theta))/yaw_rate;

       std::normal_distribution<double> dist_x(new_x,std_pos[0]);
       std::normal_distribution<double> dist_y(new_y,std_pos[1]);
       std::normal_distribution<double> dist_theta(new_theta,std_pos[2]);

       particles[i].x = dist_x(gen);
       particles[i].y = dist_y(gen);
       particles[i].theta = dist_theta(gen);
     }
   }


}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted,
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each
   *   observed measurement and assign the observed measurement to this
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will
   *   probably find it useful to implement this method and use it as a helper
   *   during the updateWeights phase.
   */
   double error, error_old;

   for (int o=0; o< observations.size(); o++){
     LandmarkObs Z;
     Z = observations[o];
     error_old = 9999999999;
     error = 0;
     for (int p = 0; p< predicted.size(); p++){
       LandmarkObs prediction;
       prediction = predicted[p];
       error = dist(prediction.x, prediction.y, Z.x, Z.y);

       if (error < error_old){
         error_old  = error;
         observations[o].id = prediction.id;
         //std::cout << "observations[o].id: " << observations[o].id << "error: " << error << std::endl;
       }

     }
   }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                  vector<LandmarkObs> &observations,
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian
   *   distribution. You can read more about this distribution here:
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system.
   *   Your particles are located according to the MAP'S coordinate system.
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
   double norm, factor, exponent;
   int closest;

   factor = 1/(2*M_PI*std_landmark[0]*std_landmark[1]);

   for (int i = 0; i < num_particles; i++){
     vector<LandmarkObs> predicted;
     Particle p = particles[i];
     for (int m = 0; m < map_landmarks.landmark_list.size(); m++){
       LandmarkObs pred;


       Map::single_landmark_s landmark_temp;
       landmark_temp = map_landmarks.landmark_list[m];
       pred.id = landmark_temp.id_i;
       float l_x = landmark_temp.x_f;
       float l_y = landmark_temp.y_f;

       pred.y = ((l_y - p.y) - (l_x - p.x)*std::tan(p.theta))/
       (std::cos(p.theta) + std::sin(p.theta)*std::tan(p.theta));

       pred.x = ((l_x - p.x) + pred.y*std::sin(p.theta))/std::cos(p.theta);

       predicted.push_back(pred);
     }
     ParticleFilter::dataAssociation(predicted, observations);
     vector<int> associations;
     vector<double> sense_x;
     vector<double> sense_y;

     for (int o=0; o< observations.size(); o++){
       //std::cout << "assigned Id: " << observations[o].id << std::endl;
       associations.push_back(observations[o].id);
       sense_x.push_back(observations[o].x);
       sense_y.push_back(observations[o].y);
     }

     ParticleFilter::SetAssociations(p, associations, sense_x, sense_y);
     //exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
               //+ (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));
    double weight = 1.0;
    for (int k = 0; k < p.associations.size(); k++){
      int id = p.associations[k];
      std::cout << "partcle association: " << id << std::endl;
      double x_map;
      x_map = p.x + (std::cos(p.theta))*(p.sense_x[k]) - (std::sin(p.theta))*(p.sense_y[k]);
      double y_map;
      y_map = p.y + std::sin(p.theta)*p.sense_x[k] + std::cos(p.theta)*p.sense_y[k];
      for (int m = 0; m < map_landmarks.landmark_list.size(); m++){
          Map::single_landmark_s landmark_temp;
          landmark_temp = map_landmarks.landmark_list[m];
          int l_id = landmark_temp.id_i;
          float l_x;
          float l_y;
          if (l_id == id){
            l_x = landmark_temp.x_f;
            l_y = landmark_temp.y_f;
            exponent = -((pow(p.sense_x[k] - l_x, 2) / (2 * pow(std_landmark[0], 2)))
                      + (pow(p.sense_y[k] - l_y, 2) / (2 * pow(std_landmark[1], 2))));
            weight *= factor*std::exp(exponent);
          }
    }
    std::cout << "particle weight = " << weight << std::endl;
    particles[i].weight = weight;

   }


}
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional
   *   to their weight.
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

   vector<Particle> p_new;
   double beta;
   std::default_random_engine generator;
   std::discrete_distribution<int> distribution (0,num_particles-1);
   int index = distribution(generator);
   double max_weight = 0.0;
   for (int i = 0; i < num_particles; i++){
     max_weight = std::max(max_weight, particles[i].weight);
     std::cout << "max weight = " << particles[i].weight << std::endl;
   }

   std::normal_distribution<double> distribution2(0,2*max_weight);
   for (int i = 0; i < num_particles; i++){
     beta = beta + distribution2(generator);
     while (particles[index].weight < beta){
       beta = beta - particles[index].weight;
       index++;
       index = index % num_particles;
     }
     p_new.push_back(particles[index]);
   }
particles = p_new;
}

void ParticleFilter::SetAssociations(Particle& particle,
                                     const vector<int>& associations,
                                     const vector<double>& sense_x,
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association,
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
