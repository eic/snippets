// This code calculates the input phi angle needed to shoot the charged particle such that it hits the middle of a given DIRC bar in the magnetic field
// This prints the phi values for theta range 30-150 deg in 5 deg steps
// Author: Nilanga Wickramaarachchi
// Based on code by William Llope

#include "TMath.h"

void phi_code_dirc()
{
  int kTargetBar = -1; 
  double fTargetBarY =  0;
  double fBarsGap =  0.15; // gap between bars in y-direction (mm)
  double momentum = 6.0; // momentum of charged particle (GeV/c)
  double fRadiatorW = 34.865; // width of bar (mm)
  
  double theta_value[25];
  double phi_value[25];
  
  for(int i=0; i < 25; i++)
      {
	theta_value[i] = 30 + 5*i;
       	double theta = theta_value[i];
    
	kTargetBar = 4; // use the bar number from 0-9
	fTargetBarY = ((double)(kTargetBar-5))*(fRadiatorW + fBarsGap) + fRadiatorW/2. + fBarsGap/2.;
    
	//std::cout << "kTargetBar = " << kTargetBar << std::endl; 
	//std::cout << "fTargetBarY = "<< fTargetBarY <<std::endl;
    
	//double partheta = (180. - fRun->getTheta());   // deg
	double partheta = 180. - theta; // deg
	double distheta = fabs(partheta-90.);// deg
	//std::cout << "distheta = " << distheta << std::endl;
	double Bval = (-1.75 + (1.75-1.6)*(distheta/60.)); // MARCO 1.7T, theta scalings to handle radial component of field
	//std::cout << "Bval = " << Bval << std::endl;
    
	double pt = momentum*TMath::Sin(theta*TMath::DegToRad());
	//std::cout << "pt = " << pt << std::endl;

	double R = 1000.*pt/0.3/fabs(Bval);
	double x1 = 0;// primary vertex
	double y1 = 0;// primary vertex
	double x2 = 708.5;
	//double x2 = 770.5; // dirc radius (mm)
	double y2 = fTargetBarY;// target point !!!! HERE ASSUMES TEN BARS PER BOX
	double xa = (x2 - x1)/2.;
	double ya = (y2 - y1)/2.;
	double x0 =  x1 + xa;
	double y0 =  y1 + ya;
	double a  = sqrt(xa*xa + ya*ya);
	double b  = sqrt(R*R - a*a);
	double x3 = x0 + b*ya/a;
	double y3 = y0 - b*xa/a;
	double x4 = x0 - b*ya/a;
	double y4 = y0 + b*xa/a;
	double ang3 = atan2(x3,y3);
	// reverse order to get phi angle of perpendicular to ray from PV to (x3,y3), NEG
	double ang4 = atan2(x4,y4);
	// reverse order to get phi angle of perpendicular to ray from PV to (x4,y4), POS
	//if (acharge>0.){ phi = -ang4; } else // pos
	//if (acharge<0.){ phi = -ang3; } // neg

	phi_value[i] = 360 - ang4*TMath::RadToDeg();

	std::cout << phi_value[i] << "\t";
	if(i==24) std::cout << "\n" << std::endl;

	//phi = ang3;
	//std::cout<<"ang3 deg = "<< -ang3*TMath::RadToDeg() <<std::endl;
      }

  for(int i=0; i < 25; i++)
    {
      std::cout << "theta = " << theta_value[i] << "\t" << "phi = " << phi_value[i] << std::endl;
    }

}
    


