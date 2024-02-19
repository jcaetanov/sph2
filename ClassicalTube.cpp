#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;
typedef unsigned int uint;

const bool barbosa=true;

const double t0 = 0.02;
const double t_max = 1.5;
const double total_mass=15; 		//Tube configuration
const double tube_length=3.;
const double eos_K = 0.1;		//Eos parameters	
const double eos_gamma = 2.;		

const double alpha_visc = 0;

const double time_step=0.01;
const int particle_number = 1000;
const double h0_SPH = 0.04;
const double h_plot = 0.02/sqrt(particle_number/1000.);
const double c0 = sqrt(eos_K*eos_gamma*pow(total_mass/tube_length,eos_gamma-1));


//set<double> plot_times({0.08,0.16,0.28,0.48,0.6,0.76,1.2}); 			//times to plot
//vector<float> plot_times{0.06,0.2,0.52,0.8,1.2}; 			//times to plot
vector<float> plot_times{0.08,0.16,0.28,0.48,0.6,0.76}; 			//times to plot
vector<int> plot_steps;


const double RK_arr[5]={1.,1.,0.5,0.5,1.};  				//convenient aux array for RK4 integration

double Kernel(double r, double h);
double GradKernel(double r, double h);
double EoS(double rho,double K=eos_K,double gamma=eos_gamma){ 	//add a table based-EoS
	return K*pow(rho,gamma);};    				//later	for generality			
double SpeedOfSound(double rho,double K=eos_K,double gamma=eos_gamma){
	return sqrt(gamma*K*pow(rho,gamma-1));
};
double Density_Analytic(double,double);
double Velocity_Analytic(double,double);

template <class T>
int sgn(T x){ return (x > 0) - (x < 0);};

template <typename T>
bool contains(vector<T> vec, const T & elem){
    bool result = false;
    if( find(vec.begin(), vec.end(), elem) != vec.end() )
    {
        result = true;
    }
    return result;
}
int myround(double d){return (int)floor(d + 0.5);};


//____________________________Class declarations________________________//
class Particle1D{
public:
	double pos,vel,rho,mass,press,h,cs;
	double kx[5],kv[5],kh[5];
	uint index;
	
	Particle1D(){};
	Particle1D(double x,uint i,double m0, double h_ini){
		pos=x;
		vel=0;
		index=i;
		mass=m0;
		kx[0]=0;
		kv[0]=0;
		kh[0]=0;
		h = h_ini; //is changed in SPHSystem constructor
	};

	~Particle1D(){};
};

class SPHSystem{
public:
//objects
	vector<Particle1D> SPH;
	uint Npart;
	double length,fluid_mass,dt;
//--tructors	
	SPHSystem(uint,double,double);
	~SPHSystem(){};
//methods
	void   ComputeDensities();            		//compute density for all particles at once,
	void   ComputeDensities(Particle1D*);		//for a single particle,
	double ComputeDensities(double,double);		//for a position in space with another h_plot;
	double ComputeVelSmoothed(double,double);	//compute smoothed velocities
	void   ComputeVel(uint);       			//compute velocity in a RK iteration (the uint input);
	void   ComputeGradPForce(uint);			//compute pressure force for all particles in a RK iteration,
	double ComputeGradPForce(Particle1D*,uint);	//for a single particle
	void   ComputeSmoothing(uint);
	double ComputeSmoothing(Particle1D*,uint);
	double ComputeArtificialVisc(Particle1D*,Particle1D*,uint); //Goffin (3.22)
	
	void RK4_step(double dt);			//implementation of a RK4 with time step dt 
};

//____________________________MAIN________________________________//
int main(int argc,char** argv){	

	for(float t:plot_times){
		int step = myround((t-t0)/time_step);
		plot_steps.push_back(step);
		cout << step << " ";
	}

	SPHSystem fluid(particle_number, tube_length, h0_SPH/sqrt(particle_number/1000.));

	string out_rho_name(argv[1]);
	ofstream out_rho (out_rho_name+"_rho_varh.dat");
	string out_vel_name(argv[1]);
	ofstream out_vel (out_vel_name+"_vel_varh.dat");

	out_rho << "# Npart dt length total_mass h0_SPH eos_K eos_gamma alpha c0\n# "
		<< particle_number<<" "<<time_step<<" "<<tube_length<<" "<< total_mass<<" "
		<<" "<<h0_SPH<<" "<<eos_K<<" "<<eos_gamma<<" "<<alpha_visc<<" " << c0 << endl;

	out_vel << "# Npart dt length total_mass h0_SPH eos_K eos_gamma c0\n# "
		<< particle_number<<" "<<time_step<<" "<<tube_length<<" "<< total_mass<<" "
		<<" "<<h0_SPH<<" "<<eos_K<<" "<<eos_gamma<<" " << c0 << endl;

	out_rho << "# plot_times\n#";
	out_vel << "# plot_times\n#";
	for(float t:plot_times){
		out_rho << " " << t;
		out_vel << " " << t;
	}

	for(float t=t0; t<t_max; t+=time_step){

		if(contains(plot_steps,myround(((t-t0)/time_step)))){
			cout << t << endl;
			out_rho << endl; 
			out_vel << endl;
//plot particles
			for(uint i=0;i<fluid.Npart;i++){
//				out_rho << fluid.SPH[i].pos << " " <<  fluid.SPH[i].rho<< " ";
				out_vel <<fluid.SPH[i].pos << " "  << fluid.SPH[i].vel << " ";
			}
//plot in space
			for(double x=-tube_length/2;x<tube_length;x+=0.005){
				out_rho << x << " " << fluid.ComputeDensities(x,h_plot) << " ";
//				out_vel << x << " " << fluid.ComputeVelSmoothed(x,h_plot) << " ";
			}
		}
		fluid.RK4_step(time_step);
	}

	out_rho.close();

	return 0;
}

//___________________Constructors_______________________________//
SPHSystem::SPHSystem(uint N=1000,double len=2.,double h0=0.04){
	Npart=N;
	length = len;
	dt = time_step;
	fluid_mass = total_mass;

	double dx = (len + 2*c0*t0/(eos_gamma-1))/N;
	double m0 = fluid_mass/N;
	double h_ini = h0/sqrt(N/1000.);

	if(barbosa==false)
		for(uint i=0;i<Npart;i++){
			Particle1D* part = new Particle1D(-len+i*dx,i,m0,h_ini);
			SPH.push_back(*part);
			delete part;
		}
	else
		for(uint i=0;i<Npart;i++){
			double z = -len+i*dx;
			double m_i = Density_Analytic(z,t0)*dx;
			Particle1D* part = new Particle1D(z,i,m_i,h_ini);
			SPH.push_back(*part);
			delete part;
		}
	cout << endl;

	ComputeDensities();
	for(uint i=0;i<Npart;i++)
		if(SPH[i].rho != 0){
			SPH[i].vel = Velocity_Analytic(SPH[i].pos,t0);
			SPH[i].h = 25*SPH[i].mass/(4.*SPH[i].rho);
		}
}

//____________________SPH method______________________________//
void SPHSystem::ComputeDensities(){
	for(uint i=0;i<Npart;i++){
		ComputeDensities(&SPH[i]);
		SPH[i].press = EoS(SPH[i].rho);
		SPH[i].cs = SpeedOfSound(SPH[i].rho);
	}
}

void SPHSystem::ComputeDensities(Particle1D* p){
	p->rho=0;
	for(uint j=0;j<Npart;j++){
		double hab = (p->h+SPH[j].h)/2;
        	p->rho += SPH[j].mass*Kernel(p->pos-SPH[j].pos,hab);
	}
}

double SPHSystem::ComputeDensities(double x,double h2){
        double rho=0;
	for(uint j=0;j<Npart;j++){
		//double hab = max(h2,SPH[j].h);
		double hab = SPH[j].h;     	
		rho += SPH[j].mass*Kernel(x-SPH[j].pos,hab);
	}
	return rho;
}

double SPHSystem::ComputeVelSmoothed(double x,double h2){		//compute smoothed velocity
        double vel=0;
        for(uint j=0;j<Npart;j++){
                //double hab = max(h2,SPH[j].h); 
                double hab = SPH[j].h; 
		vel += (SPH[j].mass/SPH[j].rho)*SPH[j].vel*Kernel(x-SPH[j].pos,hab);
	}
        return vel;
}


void SPHSystem::ComputeVel(uint RK_it=1){
	for(uint i=0;i<Npart;i++)
		SPH[i].kx[RK_it]=SPH[i].vel + SPH[i].kv[RK_it-1]*RK_arr[RK_it]*dt;
}

void SPHSystem::ComputeGradPForce(uint RK_it=1){
        for(uint i=0;i<Npart;i++)
                SPH[i].kv[RK_it]=ComputeGradPForce(&SPH[i],RK_it); 
}

double SPHSystem::ComputeGradPForce(Particle1D* pa, uint RK_it){ 
	double force=0;
	double xa = pa->pos + pa->kx[RK_it-1]*RK_arr[RK_it]*dt;
	for(uint j=0;j<Npart;j++){
        	Particle1D* pb = &SPH[j];
		double xb = pb->pos + pb->kx[RK_it-1]*RK_arr[RK_it]*dt;
		double hab = (pa->h + pa->kh[RK_it-1]*RK_arr[RK_it]*dt + pb->h + pb->kh[RK_it-1]*RK_arr[RK_it]*dt)/2;
		if (pa->rho!=0 && pb->rho!=0){
			double pressure_ab = pa->press/pow(pa->rho,2.) + pb->press/pow(pb->rho,2.);
			//double visc_ab = ComputeArtificialVisc(pa,pb,RK_it);
			//if(visc_ab > 0.01) cout << xa << " " << xb << " " << pressure_ab << " " << visc_ab << endl;
			force -= pb->mass*pressure_ab*GradKernel(xa-xb,hab);
		
		}
	}
	return force;
}

double SPHSystem::ComputeArtificialVisc(Particle1D* pa,Particle1D* pb,uint RK_it){
	double Pi_ab=0;
	double xa = pa->pos + pa->kx[RK_it-1]*RK_arr[RK_it]*dt, va = pa->vel + pa->kv[RK_it-1]*RK_arr[RK_it]*dt;
	double xb = pb->pos + pb->kx[RK_it-1]*RK_arr[RK_it]*dt, vb = pb->vel + pb->kv[RK_it-1]*RK_arr[RK_it]*dt;
	if((xa-xb)*(va-vb)<0){
		double hab = (pa->h + pa->kh[RK_it-1]*RK_arr[RK_it]*dt + pb->h + pb->kh[RK_it-1]*RK_arr[RK_it]*dt)/2;
		double mu_ab = hab*(va-vb)*(xa-xb)/( (xa-xb)*(xa-xb) + 0.01*hab*hab );
		double rho_ab = (pa->rho + pb->rho)/2., cs_ab = (pa->cs + pb->cs)/2.;
		Pi_ab -= alpha_visc*cs_ab*mu_ab/rho_ab;

	}
	return Pi_ab;
}

void SPHSystem::ComputeSmoothing(uint RK_it=1){
        for(uint i=0;i<Npart;i++)
                SPH[i].kh[RK_it]=ComputeSmoothing(&SPH[i],RK_it); 
}

double SPHSystem::ComputeSmoothing(Particle1D* pa, uint RK_it){
	if(0 == pa->rho) return 0;

	double sum=0;
        double xa = pa->pos + pa->kx[RK_it-1]*RK_arr[RK_it]*dt;
        double va = pa->vel + pa->kv[RK_it-1]*RK_arr[RK_it]*dt;
	double ha = pa->h + pa->kh[RK_it-1]*RK_arr[RK_it]*dt;
	for(uint j=0;j<Npart;j++){
		Particle1D& pb = SPH[j];
                double xb = pb.pos + pb.kx[RK_it-1]*RK_arr[RK_it]*dt;
                double vb = pb.vel + pb.kv[RK_it-1]*RK_arr[RK_it]*dt;
                double hab = (ha+pb.h + pb.kh[RK_it-1]*RK_arr[RK_it]*dt)/2;
		sum += pb.mass*(va-vb)*GradKernel(xa-xb,hab);
	}
	return -1*ha*sum/pa->rho;
}


//____________________________________RK4 step___________________________________//

void SPHSystem::RK4_step(double dt){
	for(uint rk=1;rk<=4;rk++){
		ComputeVel(rk);
		ComputeGradPForce(rk);
		ComputeSmoothing(rk);
	}
	
	for(uint i=0;i<Npart;i++){
		double dx=0, dv=0, dh=0;
		for(uint rk=1;rk<=4;rk++){
			dx += (dt/6.)*SPH[i].kx[rk]/RK_arr[rk];
			dv += (dt/6.)*SPH[i].kv[rk]/RK_arr[rk];
			dh += (dt/6.)*SPH[i].kh[rk]/RK_arr[rk];
		}
		SPH[i].pos += dx;
		SPH[i].vel += dv;
		SPH[i].h += dh;
	}

	ComputeDensities();
}

//_____________________________Kernel functions_________________________________//

double Kernel(double r, double h){
	double q=abs(r)/h;
	if(q<=1)
		return (pow(2-q,3.) - 4*pow(1-q,3.))/(6.*h);
	else if(q<=2)
		return pow(2-q,3.)/(6.*h);
	else 
		return 0;
}


double GradKernel(double r, double h){ 
        double q=abs(r)/h;
	if(q<=1)    
        	return (sgn(r)/(6.*h*h))*(-3*pow(2-q,2.) + 12*pow(1-q,2.));
        else if(q<=2)
                return (sgn(r)/(6.*h*h))*(-3*pow(2-q,2.));
        else
                return 0;
}

//____________________________Barbosa method_________________________________//


double Density_Analytic(double z, double t){	
	double rho=0;
	if(t==0){
		if(z<0)
			rho=total_mass/tube_length;
	}
	else if(z<=-c0*t){
		rho=total_mass/tube_length;
	}

	else if(z<2*c0*t/(eos_gamma-1)){ 
		double a = 2*c0/(eos_gamma+1) - (z/t)*(eos_gamma-1)/(eos_gamma+1);
		double b = a*a/(eos_K*eos_gamma);
		rho = pow(b,1/(eos_gamma-1));
	}
	return rho;
}

double Velocity_Analytic(double z, double t){
        double vel=0;
	
	if(t>0)
              if (-c0*t<z)  
		      vel=2*(c0+z/t)/(eos_gamma+1);
        return vel;
}

