#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <list>
#include <iterator>
#include "vec.h"
#include "mol.h"
using namespace std;

/* Configruation */
float rCut = pow (2., 1./6.);
float detaT = 0.005;
int stepAvg = 100;
float density = 0.8;
float temperature = 1;

/* Dimensions of the system*/
int height = 10;
int width = 10;
int length = 10;


/* Number of molecules (height*width*length)*/
int total = 1000;
/* Array for number of Molecules */
Mol mole[1000];
/* List for each Molecule for neighbour molecules position*/
list <int> neig[1000];


/* Global parameters */
float regionX, regionY, regionZ;
int steps ;
float uSum, virSum;
float kinEnergy, totEnergy, pressure;

/* declaration of functions */
void setup();
void initVel();
void simulate();
void checkBoundary();
void computeForce();
void update();
void evaluate();
void PrintSummary();
void findNeighbours();
float timeNow;


/* Main function */
int main(int argc, char const *argv[])
{
	cout << "steps | timeNow | kinEnergy | totEnergy | pressure" << endl;

	srand (time(NULL));
	steps = 0;
	setup();
	int loop = 1;

	while(loop)
	{
		simulate();
		if (steps > 1000)
		{		loop = 0; }
	}

	return 0;
}


void setup()
{
	/* Inititalization of positions ans scaling it to region*/
	int index;
	for(int i =  0; i< width; i++){
        for(int j =0; j< height; j++){
            for(int k=0; k< length; k++){
                index = k+j*length+i*length*height;
                mole[index] = Mol(i+0.5, j+0.5, k+0.5);
                mole[index].pos.x *= sqrt(1/density);
				mole[index].pos.y *= sqrt(1/density);
				mole[index].pos.z *= sqrt(1/density);
            }
        }
    }

	/* creating region from height width length */
	regionX = sqrt(1/density)*width;
	regionY = sqrt(1/density)*height;
	regionZ = sqrt(1/density)*length;
	initVel();
}


/* Inititalization of velocities */
void initVel()
{
	for (int j = 0; j < total; ++j)
	{
		float a, b, c, mag;
		a = rand()%200-100;
		b = rand()%200-100;
		c = rand()%200-100;
		mag = sqrt(a*a + b*b + c*c);
		float velMag = sqrt (2 * (1. - 1. / total) * temperature);
		// cout << mag << endl;
		a = (a/mag)*velMag;
		b = (b/mag)*velMag;
		c = (c/mag)*velMag;
		mole[j].vel.setXYZ(a,b,c);
	}
}


void simulate()
{
	++steps;
	timeNow = steps * detaT;
	update();
	evaluate();
	if (steps % stepAvg == 0) {
        PrintSummary ();
      }
}

void checkBoundary()
{
	for (int j = 0; j < total; ++j)
	{
		createVector m = mole[j].getPos();
		if(m.getX() < 0)
		{
			mole[j].pos.x += regionX ;
		} else
		if(m.getX() > regionX)
		{
			mole[j].pos.x -= regionX;
		}

		if(m.getY() < 0)
		{
			mole[j].pos.y += regionY ;
		} else
		if(m.getY() > regionY)
		{
			mole[j].pos.y -= regionY ;
		}

		if(m.getZ() < 0)
		{
			mole[j].pos.z += regionZ ;
		} else
		if(m.getZ() > regionZ)
		{
			mole[j].pos.z -= regionZ ;
		}

	}
}

void findNeighbours()
{
	createVector dr;

	for (int j1 = 0; j1 < total - 1; j1 ++) {
      for (int j2 = j1 + 1; j2 < total; j2 ++) {
      	dr.setX(mole[j1].pos.getX()-mole[j2].pos.getX());
      	dr.setY(mole[j1].pos.getY()-mole[j2].pos.getY());
      	dr.setZ(mole[j1].pos.getZ()-mole[j2].pos.getZ());

      	if(fabs(dr.getX()) > 0.5*regionX)
      	{	if (dr.getX() > 0)
      		{
      			dr.setX(-(regionX - dr.getX()));
      		} else if (dr.getX() < 0)
      		{
      			dr.setX(regionX + dr.getX());
      		}

      	}


      	if(fabs(dr.getY()) > 0.5*regionY)
      	{	if (dr.getY() > 0)
      		{
      			dr.setY(-(regionY - dr.getY()));
      		} else if (dr.getY() < 0)
      		{
      			dr.setY(regionY + dr.getY());
      		}
      	}

      	if(fabs(dr.getZ()) > 0.5*regionZ)
      	{	if (dr.getZ() > 0)
      		{
      			dr.setZ(-(regionZ - dr.getZ()));
      		} else if (dr.getZ() < 0)
      		{
      			dr.setZ(regionZ+ dr.getZ());
      		}
      	}

      	if (dr.getMagnitude() < 3)
      	{
      		neig[j1].push_back(j2);
      	}

      }
    }
}



void computeForce()
{
	createVector dr;
	uSum = 0;
	virSum = 0;

	for (int j = 0; j < total; ++j)
	{
		mole[j].acc.setX(0);
		mole[j].acc.setY(0);
		mole[j].acc.setZ(0);
	}

	for (int j1 = 0; j1 < total - 1; j1 ++) {
		list <int> :: iterator k;
      for (k = neig[j1].begin(); k != neig[j1].end(); ++k) {
      	int j2 = *k;

      	dr.setX(mole[j1].pos.getX()-mole[j2].pos.getX());
      	dr.setY(mole[j1].pos.getY()-mole[j2].pos.getY());
      	dr.setZ(mole[j1].pos.getZ()-mole[j2].pos.getZ());

      	if(fabs(dr.getX()) > 0.5*regionX)
      	{	if (dr.getX() > 0)
      		{
      			dr.setX(-(regionX - dr.getX()));
      		} else if (dr.getX() < 0)
      		{
      			dr.setX(regionX + dr.getX());
      		}

      	}


      	if(fabs(dr.getY()) > 0.5*regionY)
      	{	if (dr.getY() > 0)
      		{
      			dr.setY(-(regionY - dr.getY()));
      		} else if (dr.getY() < 0)
      		{
      			dr.setY(regionY + dr.getY());
      		}
      	}

      	if(fabs(dr.getZ()) > 0.5*regionZ)
      	{	if (dr.getZ() > 0)
      		{
      			dr.setZ(-(regionZ - dr.getZ()));
      		} else if (dr.getZ() < 0)
      		{
      			dr.setZ(regionZ + dr.getZ());
      		}
      	}



      	float rr = pow(dr.getMagnitude(),2);


      	if (rr < pow(rCut,2))
      	{
          float rri = 1. / rr;
	        float rri3 = pow(rri,3);
	        float fcVal = 48. * rri3 * (rri3 - 0.5) * rri;
	        rr =sqrt(rr);
	        mole[j1].acc.setX(mole[j1].acc.getX() + fcVal*dr.getX()/rr);
	        mole[j1].acc.setY(mole[j1].acc.getY() + fcVal*dr.getY()/rr);
	        mole[j1].acc.setZ(mole[j1].acc.getZ() + fcVal*dr.getZ()/rr);
	        mole[j2].acc.setX(mole[j2].acc.getX() - fcVal*dr.getX()/rr);
	        mole[j2].acc.setY(mole[j2].acc.getY() - fcVal*dr.getY()/rr);
	        mole[j2].acc.setZ(mole[j2].acc.getZ() - fcVal*dr.getZ()/rr);

	        uSum += 4. * rri3 * (rri3 - 1.) + 1.; //Potential Energy
        	virSum += fcVal * rr; // Viral sum
      	}
      }
  }

}

/* leapfrogStep */
void update()
{
	for (int i = 0; i < total; ++i)
	{
		mole[i].pos.setX(mole[i].pos.getX() + detaT*mole[i].vel.getX() + 0.5*detaT*detaT*mole[i].acc.getX());
		mole[i].pos.setY(mole[i].pos.getY() + detaT*mole[i].vel.getY() + 0.5*detaT*detaT*mole[i].acc.getY());
		mole[i].pos.setZ(mole[i].pos.getZ() + detaT*mole[i].vel.getZ() + 0.5*detaT*detaT*mole[i].acc.getZ());

		mole[i].vel.setX(mole[i].vel.getX() + 0.5*detaT*mole[i].acc.getX());
		mole[i].vel.setY(mole[i].vel.getY() + 0.5*detaT*mole[i].acc.getY());
		mole[i].vel.setZ(mole[i].vel.getZ() + 0.5*detaT*mole[i].acc.getZ());
	}
	checkBoundary();
	if(steps%10 == 0)
	{
		for(int q =0; q < total; q++)
		{
			neig[q].clear();
		}
		findNeighbours();
	}

	computeForce();

	for (int i = 0; i < total; ++i)
	{
		mole[i].vel.setX(mole[i].vel.getX() + 0.5*detaT*mole[i].acc.getX());
		mole[i].vel.setY(mole[i].vel.getY() + 0.5*detaT*mole[i].acc.getY());
		mole[i].vel.setZ(mole[i].vel.getZ() + 0.5*detaT*mole[i].acc.getZ());

	}


}

void evaluate()
{	float sqVel = 0;
	for (int j = 0; j < total; ++j)
	{
		sqVel += pow(mole[j].vel.getMagnitude(), 2);
	}

	kinEnergy = 0.5*sqVel/total;
	totEnergy = kinEnergy + uSum/total;
	pressure = density*(sqVel + virSum)/(total*2);
}

void PrintSummary()
{
	printf("%4d     %4.2f       %4.3f       %4.3f      %4.3f \n",steps, timeNow, kinEnergy, totEnergy, pressure);
}
