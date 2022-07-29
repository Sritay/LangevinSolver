/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *    This program simulates a point particle infront of
 *    an LJ wall in a fluid using the Langevin equations
 *    of motion in 1-D.
 *    % to compile: g++ lngevin_solver.cpp
 *    % to run: ./a.out
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <random>
using namespace std;

/*units

mass = amu (grams/mol)
length = angstrom
time = femtosecond
temperature = K
energy = (angstrom/femtosecond)^2 * (gram/mol)
force = (angstrom/femtosecond^2) * (gram/mol)
pressure = gram/(mol*angstrom*femtosecond^2)
viscosity = gram/(mol*angstrom*femtosecond)

*/
// conversion factors
double energy_real2native = 0.00041868; // lammps real units to native units
double viscosity_si2native = 60.221418;
double pressure_si2native = 6.0221418e-14;
double force_native2real = 2388.459;
double force_real2native = 0.00041868;

// global variables
double runtime = 5e7;
double dt = 2;
double eps = 15 * energy_real2native;
double sig = 17.81;
double wall_pos = 0;
double T = 600;
double eta = 0.855e-3 * viscosity_si2native; // viscosity
const double pi = 3.14159265359;
const double boltzmann = 8.3145422e-7;

default_random_engine generator(time(NULL)); // random number generator (seed)
normal_distribution<double> distribution(0.0, 0.341);
double damp = 1e3;

// particle properties
double m = 790 * 12.0107; // mass in amu
double r = 20.34;		  // radius, angstrom

vector<double> ppx; // vector of interpolation positions
vector<double> ppf; // vector of corresponding forces
/*********************************************************************/
double Evaluate_force(double, double);
double interpolate(double);

int main()
{

	// read-in particle wall force interpolation file
	double tx, tf; // temporary value holders
	ifstream pwf("wallf.txt", ios::in);
	while (true)
	{
		pwf >> tx >> tf;
		if (pwf.eof())
			break;
		ppx.push_back(tx);
		ppf.push_back(tf);
	}

	// output file
	ofstream dumpwrite("output.txt", ios::out);
	int c = 1000; // write every 1000 timesteps

	// initialise point particle
	double x = 11, y = 0, z = 0; // position, angstrom
	double v = 0;				 // velocity
	double f = 0;				 // force on particle

	// integrate particle position - velocity verlet
	double vhalf, vnext, anext, xnext;

	for (long n = 0; n < runtime; n++)
	{
		vhalf = v + 0.5 * f / m * dt;
		xnext = x + vhalf * dt;
		anext = Evaluate_force(xnext, v) / m;
		vnext = vhalf + 0.5 * anext * dt;

		// x=x-0.1;
		x = xnext;
		v = vnext;
		f = m * anext;
		if (x > 39)
			break;

		// output parameters every c timesteps
		if (n % c == 0)
		{
			dumpwrite<<n<<'\t'<<x<<'\t'<<v<<'\t'<<f*force_native2real<<endl;
		}
	}

	dumpwrite.close();
}

double Evaluate_force(double x, double v)
{
	// calculate LJ force
	double dx = x - wall_pos;
	double flj1, flj2, flj;
	if (dx < 34) // artificial cutoff at 14 angstroms
	{
		flj2 = -6 * pow(sig, 6) * pow(dx, -7);	 // attractive, towards -ve x direction
		flj1 = 12 * pow(sig, 12) * pow(dx, -13); // repulsive, towards +ve x direction
		flj = 4 * eps * (flj1 + flj2);
	}
	else
		flj = 0;

	// calculate interpolated force
	double finterp;
	if (dx < 38)
	{
		finterp = interpolate(dx);
	}
	else
		finterp = 0;

	// calculate drag force
	double fdrag = -6 * pi * eta * r * v;

	// calculate brownian force
	double random_n = distribution(generator);
	double fbr = random_n * sqrt(boltzmann * T * m / (dt * damp)); // random number

	// added force
	double fadd = 1.60 * force_real2native;
	double ftotal = finterp + fdrag + fadd + fbr;
	// double ftotal = flj + fdrag + fbr;
	return ftotal;
}

double interpolate(double x)
{
	int i = 0; // find left end of interval for interpolation
	while (x < ppx[i + 1])
		i++;
	double xR = ppx[i], yR = ppf[i], xL = ppx[i + 1], yL = ppf[i + 1]; // points on either side (unless beyond ends)

	double dydx = (yR - yL) / (xR - xL); // gradient
	return yL + dydx * (x - xL);		 // linear interpolation
}
