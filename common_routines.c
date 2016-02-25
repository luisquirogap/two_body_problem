
#define mypow2(x) (x*x)

int acceleration(void)
{
  
  double r, dx, dy, dz, r3;
  
  aux_part[0].accel[0] = 0.0;
  aux_part[0].accel[1] = 0.0;
  aux_part[0].accel[2] = 0.0;
  
  aux_part[1].accel[0] = 0.0;
  aux_part[1].accel[1] = 0.0;
  aux_part[1].accel[2] = 0.0;
  
  galaxies[0].pot = galaxies[1].pot = 0.0;
  
  /* COMPUTING ACCELERATIONS FOR THE HOST */
  
  dx = aux_part[0].pos[0] - aux_part[1].pos[0];
  dy = aux_part[0].pos[1] - aux_part[1].pos[1];
  dz = aux_part[0].pos[2] - aux_part[1].pos[2];
  
  r = sqrt ( mypow2(dx) + mypow2(dy) + mypow2(dz) );
  
  galaxies[0].pot = - G*galaxies[1].mass/r;
  
  r3 = r*r*r;

  aux_part[0].accel[0] = - G*galaxies[1].mass* ( dx/r3 );
  aux_part[0].accel[1] = - G*galaxies[1].mass* ( dy/r3 );
  aux_part[0].accel[2] = - G*galaxies[1].mass* ( dz/r3 );
  
  
  /* COMPUTING ACCELERATIONS FOR THE SATELLITE */

  galaxies[1].pot = - G*galaxies[0].mass/r;
  
  aux_part[1].accel[0] = - G*galaxies[0].mass* ( - dx/r3 );
  aux_part[1].accel[1] = - G*galaxies[0].mass* ( - dy/r3 );
  aux_part[1].accel[2] = - G*galaxies[0].mass* ( - dz/r3 );
  
  adap_timestep[0].accel = sqrt( aux_part[0].accel[0]*aux_part[0].accel[0] + aux_part[0].accel[1]*aux_part[0].accel[1] + aux_part[0].accel[2]*aux_part[0].accel[2] ); 
  adap_timestep[1].accel = sqrt( aux_part[1].accel[0]*aux_part[1].accel[0] + aux_part[1].accel[1]*aux_part[1].accel[1] + aux_part[1].accel[2]*aux_part[1].accel[2] ); 
  
  return 0;
}

int move_extended_body(void)
{
  int i;

  double dx, dy, dz;
   
  dx = aux_part[0].pos[0] - galaxies[0].pos[0];
  dy = aux_part[0].pos[1] - galaxies[0].pos[1];
  dz = aux_part[0].pos[2] - galaxies[0].pos[2];

  for(i=0; i<N_particles_set; i++)
    {
      particles_set[i].pos[0] = particles_set[i].pos[0] + dx;
      particles_set[i].pos[1] = particles_set[i].pos[1] + dy;
      particles_set[i].pos[2] = particles_set[i].pos[2] + dz;
    }
  
  return 0;
}

int system_center_of_mass(double poscm[3], double velcm[3])
{
  
  double Mtotal;
  
  poscm[0] = 0.0;
  poscm[1] = 0.0;
  poscm[2] = 0.0;
  
  velcm[0] = 0.0;
  velcm[1] = 0.0;
  velcm[2] = 0.0;

  poscm[0] = poscm[0] + galaxies[0].pos[0]*galaxies[0].mass + galaxies[1].pos[0]*galaxies[1].mass;
  poscm[1] = poscm[1] + galaxies[0].pos[1]*galaxies[0].mass + galaxies[1].pos[1]*galaxies[1].mass;
  poscm[2] = poscm[2] + galaxies[0].pos[2]*galaxies[0].mass + galaxies[1].pos[2]*galaxies[1].mass;
  
  velcm[0] = velcm[0] + galaxies[0].vel[0]*galaxies[0].mass + galaxies[1].vel[0]*galaxies[1].mass;
  velcm[1] = velcm[1] + galaxies[0].vel[1]*galaxies[0].mass + galaxies[1].vel[1]*galaxies[1].mass;
  velcm[2] = velcm[2] + galaxies[0].vel[2]*galaxies[0].mass + galaxies[1].vel[2]*galaxies[1].mass;

  Mtotal = galaxies[0].mass + galaxies[1].mass;
  
  poscm[0] = poscm[0]/Mtotal;
  poscm[1] = poscm[1]/Mtotal;
  poscm[2] = poscm[2]/Mtotal;
  
  velcm[0] = velcm[0]/Mtotal;
  velcm[1] = velcm[1]/Mtotal;
  velcm[2] = velcm[2]/Mtotal;
    
  return 0;
  
}

int translateSystem(double poscm[3], double velcm[3])
{
  galaxies[0].pos[0] = galaxies[0].pos[0] - poscm[0]; 
  galaxies[0].pos[1] = galaxies[0].pos[1] - poscm[1];
  galaxies[0].pos[2] = galaxies[0].pos[2] - poscm[2];
  
  galaxies[0].vel[0] = galaxies[0].vel[0] - velcm[0]; 
  galaxies[0].vel[1] = galaxies[0].vel[1] - velcm[1];
  galaxies[0].vel[2] = galaxies[0].vel[2] - velcm[2];
  
  galaxies[1].pos[0] = galaxies[1].pos[0] - poscm[0]; 
  galaxies[1].pos[1] = galaxies[1].pos[1] - poscm[1];
  galaxies[1].pos[2] = galaxies[1].pos[2] - poscm[2];
  
  galaxies[1].vel[0] = galaxies[1].vel[0] - velcm[0]; 
  galaxies[1].vel[1] = galaxies[1].vel[1] - velcm[1];
  galaxies[1].vel[2] = galaxies[1].vel[2] - velcm[2];
  
  return 0;
}

double energymin,energymax;

int conservation(double time, int Npart, FILE *file,int Tinit)
{
  double Ek,Ep,totalAngularMom[3],totalLinealMom[3];
  double P,L;
  double vel02,vel12,d;
  
  //Satellite lineal momentum
  totalLinealMom[0] = galaxies[0].mass*galaxies[0].vel[0] + galaxies[1].mass*galaxies[1].vel[0];
  totalLinealMom[1] = galaxies[0].mass*galaxies[0].vel[1] + galaxies[1].mass*galaxies[1].vel[1];
  totalLinealMom[2] = galaxies[0].mass*galaxies[0].vel[2] + galaxies[1].mass*galaxies[1].vel[2];
  
  P = sqrt( totalLinealMom[0]*totalLinealMom[0] 
		+ totalLinealMom[1]*totalLinealMom[1] 
		+ totalLinealMom[2]*totalLinealMom[2]  );

  //Satellite angular momentum
  totalAngularMom[0] = galaxies[1].mass*(galaxies[1].pos[1]*galaxies[1].vel[2]-galaxies[1].pos[2]*galaxies[1].vel[1]);
  totalAngularMom[1] = - galaxies[1].mass*(galaxies[1].pos[0]*galaxies[1].vel[2]-galaxies[1].pos[2]*galaxies[1].vel[0]);
  totalAngularMom[2] = galaxies[1].mass*(galaxies[1].pos[0]*galaxies[1].vel[1]-galaxies[1].pos[1]*galaxies[1].vel[0]);
  
  totalAngularMom[0] = totalAngularMom[0] + galaxies[0].mass*(galaxies[0].pos[1]*galaxies[0].vel[2]-galaxies[0].pos[2]*galaxies[0].vel[1]);
  totalAngularMom[1] = totalAngularMom[1] - galaxies[0].mass*(galaxies[0].pos[0]*galaxies[0].vel[2]-galaxies[0].pos[2]*galaxies[0].vel[0]);
  totalAngularMom[2] = totalAngularMom[2] + galaxies[0].mass*(galaxies[0].pos[0]*galaxies[0].vel[1]-galaxies[0].pos[1]*galaxies[0].vel[0]);
  
    
  L = sqrt( totalAngularMom[0]*totalAngularMom[0] 
	    + totalAngularMom[1]*totalAngularMom[1] 
	    +totalAngularMom[2]*totalAngularMom[2] );

  //Satellite kinetic energy
  vel02 = galaxies[1].vel[0]*galaxies[1].vel[0] 
    + galaxies[1].vel[1]*galaxies[1].vel[1] 
    + galaxies[1].vel[2]*galaxies[1].vel[2];

  vel12 = galaxies[0].vel[0]*galaxies[0].vel[0] 
    + galaxies[0].vel[1]*galaxies[0].vel[1] 
    + galaxies[0].vel[2]*galaxies[0].vel[2];      

  Ek = 0.5*galaxies[0].mass*vel02 + 0.5*galaxies[1].mass*vel12;

  d = sqrt( pow((galaxies[1].pos[0]-galaxies[0].pos[0]),2)
	    +pow((galaxies[1].pos[1]-galaxies[0].pos[1]),2)
	    +pow((galaxies[1].pos[2]-galaxies[0].pos[2]),2) );

  Ep = - G*galaxies[0].mass*galaxies[1].mass/d;  

  if(Tinit==1) //if t=0
    {
      initLinealMom = P;
      initAngularMom = L;
      initEnergy = Ek+Ep;
      printf("Lineal momentum : %e %e %e\n",totalLinealMom[0],totalLinealMom[1],totalLinealMom[2]);
      printf("Initial energy : %e\n",initEnergy);
      energymin = initEnergy; 
      energymax = initEnergy; 
    }

  if( Ek+Ep < energymin)
    energymin = Ek+Ep; 

  if( Ek+Ep > energymax)
    energymax = Ek+Ep; 

  fprintf(file,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",time,
	  totalLinealMom[0],totalLinealMom[1],totalLinealMom[2],P,
	  totalAngularMom[0],totalAngularMom[1],totalAngularMom[2],L,
	  Ek,Ep,Ek+Ep,fabs((P - initLinealMom)/initLinealMom),fabs((L - initAngularMom)/initAngularMom),fabs((Ek+Ep - initEnergy)/initEnergy) );
 
 return 0;
}


