#ifndef __FANGOOSTERLEE_H_INCLUDED__
#define __FANGOOSTERLEE_H_INCLUDED__
#define _USE_MATH_DEFINES
#include <cmath>
#include "Complex.h"
#include <vector>
#include <unordered_map>
#include <iostream>
#include <string>
#include <functional>
//typedef Complex (*cf)(Complex, std::map<std::string, double>); //defines cf as a pointer to a function which takes complex and outputs complex as arguments...as of now all arguments must be doubles...

class FangOosterlee {
	private:
		int k;
		int h;
		double exloss;
		double vloss;
		//double M_PI;
	public:
		FangOosterlee(int, int);
		//template< typename FN, typename... ARGS>
		//std::map<std::string, std::vector<double> > computeDistribution(double, double, FN&& fn, ARGS&&... args);
		double getEL();
		double getVariance();
		template< typename FN, typename... ARGS>
		std::unordered_map<std::string, std::vector<double> > computeDistribution(double xmin, double xmax, FN&& fn, ARGS&&... args) {
			double xRange=xmax-xmin;
			double du=M_PI/xRange;
			double dx=xRange/(double)(h-1);
			double cp=2.0/xRange;
			std::vector<double> f=std::vector<double> (k);
			std::vector<double> y=std::vector<double> (h);

			std::vector<double> x=std::vector<double> (h);
			std::vector<double> VaR=std::vector<double> (h);

			#pragma omp parallel//multithread using openmp
			{
				#pragma omp for //multithread using openmp
				for(int j=0; j<k; j++){
					Complex u=Complex(0, du*j);

					f[j]=fn(u, args...).multiply(u.multiply(-xmin).exp()).getReal()*cp;

				}
			}
			f[0]=.5*f[0];
			exloss=0;
			vloss=0;
			double cdf=0;
			for(int i=0;  i<h; i++){
				y[i]=0;
				x[i]=xmin+dx*i;
				for(int j=0; j<k; j++){
					y[i]=y[i]+f[j]*std::cos(du*j*dx*i);
				}
				vloss=vloss+y[i]*i*dx*dx*i;
				exloss=exloss+y[i]*(xmin+i*dx)*dx;
				cdf=cdf+dx*y[i];
				VaR[i]=cdf;
			}
			std::unordered_map<std::string, std::vector<double> > distribution;
			distribution["x"]=x;
			distribution["y"]=y;
			distribution["VaR"]=VaR;
			return distribution;
		}
		template< typename FN, typename... ARGS>
		void computeDistributionJSON(double xmin, double xmax, FN&& fn, ARGS&&... args) {
			computeDistributionJSON(false, xmin, xmax, fn, args...);
		}
		
		template< typename FN, typename... ARGS>
		void computeDistributionJSON(bool showProgress, double xmin, double xmax, FN&& fn, ARGS&&... args) {
			double xRange=xmax-xmin;
			double du=M_PI/xRange;
			double dx=xRange/(double)(h-1);
			double cp=2.0/xRange;
			std::vector<double> f=std::vector<double> (k);
			int trackProgress=0;
			int numProgress=20;//relatively smooth if 20 updates
			int interval=k/numProgress;
			#pragma omp parallel//multithread using openmp
			{
				#pragma omp for //multithread using openmp
				for(int j=0; j<k; j++){
					Complex u=Complex(0, du*j);
					f[j]=fn(u, args...).multiply(u.multiply(-xmin).exp()).getReal()*cp;
					if(showProgress){
						if(trackProgress%interval==0){
							std::cerr<<"{\"progress\":"<<(double)trackProgress/(double)k<<"}\\n"<<std::endl;
						}
					}
					trackProgress++;
				}
			}
			f[0]=.5*f[0];
			exloss=0;
			vloss=0;
			double cdf=0;
			std::cout<<"{\"y\":[";
			double y=0;

			for(int i=0;  i<(h-1); i++){
				y=0;
				for(int j=0; j<k; j++){
					y=y+f[j]*std::cos(du*j*dx*i);
				}

				std::cout<<y<<", ";
				vloss=vloss+y*i*dx*dx*i;
				exloss=exloss+y*(xmin+i*dx)*dx;
			}
			y=0;
			for(int j=0; j<k; j++){
				y=y+f[h-1]*std::cos(du*j*dx*(h-1));
			}
			std::cout<<y<<"], \"xmin\":"<<xmin<<", \"dx\":"<<dx<<"}\\n"<<std::endl;
			vloss=vloss+y*(h-1)*dx*dx*(h-1);
			exloss=exloss+y*(xmin+(h-1)*dx)*dx;
		}
};
#endif
