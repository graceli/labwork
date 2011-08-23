#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>

//#define DEBUG
using namespace std;
//     0   1   2    3   4   5
enum {O1, HAA, C1, H1, C2, H2, O2, HAB, C3, H3, O3, HAC};

void readGroFile(ifstream& gro, double** dipepCoords){
	string title, numAtoms;
	string line;

	//get a line and parse each line using the stringstream
	int lineNum=0;
	string resname, atomname;	
	int atomnum;
	double x,y,z;	
	while(getline(gro,line)){
		//ignore comments the number of atoms, and the box numbers
		if(lineNum  < 2 || lineNum == 26){
#ifdef DEBUG
			cerr<<lineNum<<" ignored: "<<line<<endl;
#endif
			lineNum++;
			continue;
		}
		//read body: read in dipeptide coordinates
		istringstream ist(line);
#ifdef DEBUG
		cerr<<line<<endl;
#endif
		ist>>resname>>atomname>>atomnum>>x>>y>>z;
		dipepCoords[lineNum-2][0]=x;
		dipepCoords[lineNum-2][1]=y;
		dipepCoords[lineNum-2][2]=z;
		lineNum++;
	}
}

void output(double** dipepCoords){
	for(int i=0; i<24; i++){
		for(int j=0; j<3; j++){
			cout<<dipepCoords[i][j]<<" ";
		}
		cout<<endl;
	}
}
void normalize(double vec[3]){
	double norm = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
	vec[0] = vec[0]/norm;
	vec[1] = vec[1]/norm;
	vec[2] = vec[1]/norm;
}

double dot(double v1[3], double v2[3]){
	return (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]);
}

double dihedral(double* u, double* v, double* w){
	double* u_cross_v = new double[3];
	double* v_cross_w = new double[3];

	//compute u x v 
	u_cross_v[0] = u[1]*v[2]-u[2]*v[1];
	u_cross_v[1] = u[2]*v[0]-u[0]*v[2];
	u_cross_v[2] = u[0]*v[1]-u[1]*v[0];

	//normalize u x v
	double norm = sqrt(dot(u_cross_v,u_cross_v));
	for(int i=0; i<3; i++){
		u_cross_v[i]/=norm;
	}
	
	//compute v x w
	v_cross_w[0] = v[1]*w[2] - w[1]*v[2];
	v_cross_w[1] = w[0]*v[2] - v[0]*w[2];
	v_cross_w[2] = v[0]*w[1] - w[0]*v[1];

	//normalize v x w
	double normp = sqrt(dot(v_cross_w, v_cross_w));
	for(int i=0; i<3; i++){
		v_cross_w[i]/=normp;
	}
	
	//calculate dihedral angle
	double dot_prod = dot(u_cross_v, v_cross_w);
	double angle = acos(dot(u_cross_v, v_cross_w));

#ifdef DEBUG
	cerr<<"dot(u x v , v x w)="<<dot_prod<<endl;
	cerr<<"angle in rad = "<<angle<<endl;
	cerr<<"angle in deg = "<<angle*180/3.14159<<endl;
	//determine sign of angle
	double sign = dot(u_cross_v, w);
	cerr<<"sign="<<sign<<endl;
#endif

	return angle*180/3.14159;
}

int main(int argc, char* argv[]){
	double** dipepCoords = new double*[24];
	ifstream gro("/home/grace/work/inositol/gro/epi_em.gro");
	for(int i=0; i<24; i++){
		dipepCoords[i]=new double[3];
	}
	readGroFile(gro, dipepCoords);

	//define the 3 vectors for dihedral angle calculation
	double* u = new double[3];
	double* v = new double[3];
	double* w = new double[3];

	double* uax = new double[3];
	double* vax = new double[3];
	double* wax = new double[3];

	u[0] = dipepCoords[C1][0] - dipepCoords[O1][0];
	u[1] = dipepCoords[C1][1] - dipepCoords[O1][1];
	u[2] = dipepCoords[C1][2] - dipepCoords[O1][2];

	v[0] = dipepCoords[C2][0] - dipepCoords[C1][0];
	v[1] = dipepCoords[C2][1] - dipepCoords[C1][1];
	v[2] = dipepCoords[C2][2] - dipepCoords[C1][2];

	w[0] = dipepCoords[H2][0] - dipepCoords[C2][0];
	w[1] = dipepCoords[H2][1] - dipepCoords[C2][1];
	w[2] = dipepCoords[H2][2] - dipepCoords[C2][2];
	double equa_angle = dihedral(u,v,w);
	
	uax[0] = dipepCoords[C2][0] - dipepCoords[O2][0];
	uax[1] = dipepCoords[C2][1] - dipepCoords[O2][1];
	uax[2] = dipepCoords[C2][2] - dipepCoords[O2][2];

	vax[0] = dipepCoords[C3][0] - dipepCoords[C2][0];
	vax[1] = dipepCoords[C3][1] - dipepCoords[C2][1];
	vax[2] = dipepCoords[C3][2] - dipepCoords[C2][2];

	wax[0] = dipepCoords[H3][0] - dipepCoords[C3][0];
	wax[1] = dipepCoords[H3][1] - dipepCoords[C3][1];
	wax[2] = dipepCoords[H3][2] - dipepCoords[C3][2];
	double axial_angle = dihedral(uax, vax, wax);

	return 0;
}




















