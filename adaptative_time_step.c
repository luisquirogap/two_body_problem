

int assingTimeSteps_accel()
{

  deltatmin = etha/sqrt(adap_timestep[0].accel);
  
  adap_timestep[1].timestep = etha/sqrt(adap_timestep[1].accel);
  
  if(adap_timestep[1].timestep < deltatmin)
    deltatmin = adap_timestep[1].timestep; 
  
  adap_timestep[0].timestep = deltatmin;
  adap_timestep[1].timestep = deltatmin;
  
  return 0;
}
