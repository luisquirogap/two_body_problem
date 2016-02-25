
double deltatmax;
double deltatmin;
int nmax;

#include"adaptative_time_step.c"

int leapfrog(double T_end)
{

  int counter;
  double t;
  double Poscm[3],Velcm[3];     
  
  FILE *forbit0,*forbit1,*fconservation,*fcmass;
  forbit0 = fopen("orbita0.dat","w");
  forbit1 = fopen("orbita1.dat","w");
  fconservation = fopen("conservation.dat","w");
  fcmass = fopen("center_of_mass.dat","w");
  
  adap_timestep = (Adaptative_timestep *)malloc( (size_t)Npart*sizeof(Adaptative_timestep) );
  if( adap_timestep == NULL){
    printf("Allocation of adap_timestep failed\n");
    exit(0);
  }
  
  aux_part[0].pos[0] = galaxies[0].pos[0];
  aux_part[0].pos[1] = galaxies[0].pos[1];
  aux_part[0].pos[2] = galaxies[0].pos[2];
  
  aux_part[0].accel[0] = 0.0;
  aux_part[0].accel[1] = 0.0;
  aux_part[0].accel[2] = 0.0;
  
  galaxies[0].pot = 0.0;
  
  aux_part[1].pos[0] = galaxies[1].pos[0];
  aux_part[1].pos[1] = galaxies[1].pos[1];
  aux_part[1].pos[2] = galaxies[1].pos[2];
  
  aux_part[1].accel[0] = 0.0;
  aux_part[1].accel[1] = 0.0;
  aux_part[1].accel[2] = 0.0;
  
  galaxies[1].pot = 0.0;
  
 
  fprintf(forbit0,"%lf %lf %lf %lf %lf %lf %lf %lf %e %lf %lf\n",0.0,
	  galaxies[0].pos[0],galaxies[0].pos[1],galaxies[0].pos[2],
	  galaxies[0].vel[0],galaxies[0].vel[1],galaxies[0].vel[2],
	  galaxies[0].pot,adap_timestep[0].timestep,adap_timestep[0].accel,
	  sqrt(galaxies[0].pos[0]*galaxies[0].pos[0] + 
	       galaxies[0].pos[1]*galaxies[0].pos[1] + 
	       galaxies[0].pos[2]*galaxies[0].pos[2]));	  
  
  
  fprintf(forbit1,"%lf %lf %lf %lf %lf %lf %lf %lf %e %lf %lf\n",0.0,
	  galaxies[1].pos[0],galaxies[1].pos[1],galaxies[1].pos[2],
	  galaxies[1].vel[0],galaxies[1].vel[1],galaxies[1].vel[2],
	  galaxies[1].pot,adap_timestep[1].timestep,adap_timestep[1].accel,
	  sqrt(galaxies[1].pos[0]*galaxies[1].pos[0] + 
	       galaxies[1].pos[1]*galaxies[1].pos[1] + 
	       galaxies[1].pos[2]*galaxies[1].pos[2]));
  
 
  // conservation(0.0, Npart, fconservation, 1);
  counter = 0;
  
  t = 0.0;
    
  while(t<=T_end)
    {

      acceleration();
      assingTimeSteps_accel();
      dt = adap_timestep[1].timestep;
      
      aux_part[0].pos[0] = aux_part[0].pos[0] + 0.5*dt*galaxies[0].vel[0];
      aux_part[0].pos[1] = aux_part[0].pos[1] + 0.5*dt*galaxies[0].vel[1];
      aux_part[0].pos[2] = aux_part[0].pos[2] + 0.5*dt*galaxies[0].vel[2];
      
      aux_part[1].pos[0] = aux_part[1].pos[0] + 0.5*dt*galaxies[1].vel[0];
      aux_part[1].pos[1] = aux_part[1].pos[1] + 0.5*dt*galaxies[1].vel[1];
      aux_part[1].pos[2] = aux_part[1].pos[2] + 0.5*dt*galaxies[1].vel[2];
      
      acceleration();
      
      galaxies[0].vel[0] = galaxies[0].vel[0] + dt*aux_part[0].accel[0];
      galaxies[0].vel[1] = galaxies[0].vel[1] + dt*aux_part[0].accel[1];
      galaxies[0].vel[2] = galaxies[0].vel[2] + dt*aux_part[0].accel[2];
      
      galaxies[1].vel[0] = galaxies[1].vel[0] + dt*aux_part[1].accel[0];
      galaxies[1].vel[1] = galaxies[1].vel[1] + dt*aux_part[1].accel[1];
      galaxies[1].vel[2] = galaxies[1].vel[2] + dt*aux_part[1].accel[2];
      
      aux_part[0].pos[0] = aux_part[0].pos[0] + 0.5*dt*galaxies[0].vel[0];
      aux_part[0].pos[1] = aux_part[0].pos[1] + 0.5*dt*galaxies[0].vel[1];
      aux_part[0].pos[2] = aux_part[0].pos[2] + 0.5*dt*galaxies[0].vel[2];
      
      aux_part[1].pos[0] = aux_part[1].pos[0] + 0.5*dt*galaxies[1].vel[0];
      aux_part[1].pos[1] = aux_part[1].pos[1] + 0.5*dt*galaxies[1].vel[1];
      aux_part[1].pos[2] = aux_part[1].pos[2] + 0.5*dt*galaxies[1].vel[2];

      galaxies[0].pos[0] = aux_part[0].pos[0];
      galaxies[0].pos[1] = aux_part[0].pos[1];
      galaxies[0].pos[2] = aux_part[0].pos[2];
      
      galaxies[1].pos[0] = aux_part[1].pos[0];
      galaxies[1].pos[1] = aux_part[1].pos[1];
      galaxies[1].pos[2] = aux_part[1].pos[2];
                  
      fprintf(forbit0,"%lf %lf %lf %lf %lf %lf %lf %lf %e %lf %lf\n",t+dt,
	      galaxies[0].pos[0],galaxies[0].pos[1],galaxies[0].pos[2],
	      galaxies[0].vel[0],galaxies[0].vel[1],galaxies[0].vel[2],
	      galaxies[0].pot,adap_timestep[0].timestep,adap_timestep[0].accel,
	      sqrt(galaxies[0].pos[0]*galaxies[0].pos[0] + 
		   galaxies[0].pos[1]*galaxies[0].pos[1] + 
		   galaxies[0].pos[2]*galaxies[0].pos[2]));	  
            
      fprintf(forbit1,"%lf %lf %lf %lf %lf %lf %lf %lf %e %lf %lf\n",t+dt,
	      galaxies[1].pos[0],galaxies[1].pos[1],galaxies[1].pos[2],
	      galaxies[1].vel[0],galaxies[1].vel[1],galaxies[1].vel[2],
	      galaxies[1].pot,adap_timestep[1].timestep,adap_timestep[1].accel,
	      sqrt(galaxies[1].pos[0]*galaxies[1].pos[0] + 
		   galaxies[1].pos[1]*galaxies[1].pos[1] + 
		   galaxies[1].pos[2]*galaxies[1].pos[2]));
                  
      system_center_of_mass(Poscm, Velcm);
      fprintf(fcmass," %lf %lf %lf\n",Poscm[0],Poscm[1],Poscm[2]);
      
      counter++;
      
      t = t +dt;     

    }
    
  printf("ENERGY DEVIATION = %e\n",(energymax-energymin)/energymax);
  
  printf("total iterations = %d\n",counter);
  
  fclose(forbit0);  
  fclose(forbit1);  
  fclose(fconservation);
  fclose(fcmass);
  

  free(adap_timestep);

  
  return 0;
  
}


