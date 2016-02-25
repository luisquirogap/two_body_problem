
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<malloc.h>

typedef struct
{
  double pos[3];
  double vel[3];
  double mass;
  double pot;
}Particles;

typedef struct
{
  double pos[3];
  double vel[3];
  double accel[3];
}Aux_particles;

typedef struct
{
  double accel;
  double timestep;
  int ntimestep;
}Adaptative_timestep;

Particles *galaxies=NULL,*particles_set=NULL;
Aux_particles *aux_part=NULL;
Adaptative_timestep *adap_timestep = NULL;
int N_particles_set, nsteps, Npart;
double etha, dt, G = 43007.1;
double initLinealMom,initAngularMom,initEnergy;

#include"common_routines.c"
#include"leapfrog.c"

int counterLines(char *infile);

int main(int argc, char *argv[])
{

  double T_end, Poscm[3],Velcm[3];
  
  Npart = 2;
  
  galaxies = (Particles * )malloc( (size_t)Npart*sizeof(Particles) );
  if(galaxies == NULL){
    printf("Allocation of galaxies failed\n");
    exit(0);
  }

 aux_part = (Aux_particles * )malloc( (size_t)Npart*sizeof(Aux_particles) );
  if(aux_part == NULL){
    printf("Allocation of aux_part failed\n");
    exit(0);
  }
  galaxies[1].pos[0] = -21.95;
  galaxies[1].pos[1] = -12.58;
  galaxies[1].pos[2] = 1.517963;
  
  galaxies[1].vel[0] = -85.014319; 
  galaxies[1].vel[1] = 306.956602; 
  galaxies[1].vel[2] = -4.732407;

  galaxies[1].mass = 2.41;

  galaxies[0].pos[0] = 0.0;
  galaxies[0].pos[1] = 0.0;
  galaxies[0].pos[2] = 0.0;

  galaxies[0].vel[0] = 0.0; 
  galaxies[0].vel[1] = 0.0; 
  galaxies[0].vel[2] = 0.0;

  galaxies[0].mass =  103.656888;

  printf("Mass main galaxy = %lf\n",galaxies[0].mass);
  printf("Mass satellite galaxy = %lf\n",galaxies[1].mass);

  system_center_of_mass(Poscm, Velcm);
  printf("Pos center of mass: %lf %lf %lf\n",Poscm[0],Poscm[1],Poscm[2]);
  printf("Vel center of mass: %lf %lf %lf\n",Velcm[0],Velcm[1],Velcm[2]);
  
  translateSystem(Poscm, Velcm);
  printf("Pos galaxy 0: %lf %lf %lf\n",galaxies[0].pos[0],galaxies[0].pos[1],galaxies[0].pos[2]);
  printf("Pos galaxy 1: %lf %lf %lf\n",galaxies[1].pos[0],galaxies[1].pos[1],galaxies[1].pos[2]);
  printf("Vel galaxy 0: %lf %lf %lf\n",galaxies[0].vel[0],galaxies[0].vel[1],galaxies[0].vel[2]);
  printf("Vel galaxy 1: %lf %lf %lf\n",galaxies[1].vel[0],galaxies[1].vel[1],galaxies[1].vel[2]);
  system_center_of_mass(Poscm, Velcm);
  printf("Pos center of mass: %lf %lf %lf\n",Poscm[0],Poscm[1],Poscm[2]);
  printf("Vel center of mass: %lf %lf %lf\n",Velcm[0],Velcm[1],Velcm[2]);

 
  T_end = 1.0;
  etha = 1e-2;

  leapfrog(T_end);

  free(galaxies);
  free(aux_part);
  free(particles_set);

  return 0;
}

int counterLines(char *infile) 
{

  int NDAT,c;
  static FILE *pf;

  if((pf=fopen(infile,"r"))==NULL)
    printf("no puedo abrir archivo 1\n");

  NDAT=0;

  while((c=fgetc(pf))!= EOF)
    {
      if(c == '\n')
        {
          ++NDAT;
        }
    }

  fclose(pf);

  return NDAT;

}
