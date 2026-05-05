/* File: diffusion_s.c */ 

#include "defines.h"



/************************************************************************************/
void thermal_update (struct info_param param, float *grid, float *grid_chips, int tam, int des, int pid, int npr) //falta testearlo
{
  int i,i_chips, j, a, b;
  
  // heat injection at chip positions
  for (i_chips=1; i_chips<tam+1; i_chips++){
    i=i_chips+des;
    for (j=1; j<NCOL-1; j++) 
      if (grid_chips[i_chips*NCOL+j] > grid[i*NCOL+j])
        grid[i*NCOL+j] += 0.05 * (grid_chips[i_chips*NCOL+j] - grid[i*NCOL+j]);

  }

  // air cooling at the middle of the card
  a = 0.45*(NCOL-2) + 1;
  b = 0.55*(NCOL-2) + 1;

  for (i=des+1; i<des+tam+1; i++) //cambiar los indices
  for (j=a; j<b; j++)
      grid[i*NCOL+j] -= 0.01 * (grid[i*NCOL+j] - param.t_ext);
}

/************************************************************************************/
double thermal_diffusion (struct info_param param, float *grid, float *grid_aux, int tam, int des, int pid, int npr)
{
  int    i, j;
  double  T;
  double Tfull = 0.0, Tfull_loc=0.0;
  MPI_Request req1,req2;
  //pedir/enviar las lineas necesarias
  if(pid==0){ //Seguramente cambiar a ISends + Recvs
    MPI_Isend(&grid[NCOL*(tam+des)],1,fila,pid+1,0,MPI_COMM_WORLD,&req1);
    MPI_Irecv(&grid[NCOL*(tam+des+1)],1,fila,pid+1,0,MPI_COMM_WORLD,&req2);
    //MPI_Sendrecv(&grid[NCOL*(tam+des)],1,fila,pid+1,0,&grid[NCOL*(tam+des+1)],1,fila,pid+1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE); //Es posible que de problemas por usar el mismo buffer, si eso cambiar
  }
  else if(pid==npr-1){
    MPI_Isend(&grid[NCOL*(des+1)],1,fila,pid-1,0,MPI_COMM_WORLD,&req1);
    MPI_Irecv(&grid[NCOL*(des)],1,fila,pid-1,0,MPI_COMM_WORLD,&req1);
    
    //MPI_Sendrecv(&grid[NCOL*(des+1)],1,fila,pid-1,0,&grid[NCOL*(des)],1,fila,pid-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  }
  else{
    MPI_Isend(&grid[NCOL*(tam+des)],1,fila,pid+1,0,MPI_COMM_WORLD,&req1);
    MPI_Irecv(&grid[NCOL*(tam+des+1)],1,fila,pid+1,0,MPI_COMM_WORLD,&req2);
    MPI_Isend(&grid[NCOL*(des+1)],1,fila,pid-1,0,MPI_COMM_WORLD,&req1);
    MPI_Irecv(&grid[NCOL*(des)],1,fila,pid-1,0,MPI_COMM_WORLD,&req1);

    //MPI_Sendrecv(&grid[NCOL*(tam+des)],1,fila,pid+1,0,&grid[NCOL*(tam+des+1)],1,fila,pid+1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    //MPI_Sendrecv(&grid[NCOL*(des+1)],1,fila,pid-1,0,&grid[NCOL*(des)],1,fila,pid-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  }
  for (i=des+2; i<des+tam; i++)//no procesa ni la primera ni la ultima fila 
  {
    for (j=1; j<NCOL-1; j++)
    {
      T = grid[i*NCOL+j] + 
          0.10 * (grid[(i+1)*NCOL+j]   + grid[(i-1)*NCOL+j]   + grid[i*NCOL+(j+1)]     + grid[i*NCOL+(j-1)] + 
                  grid[(i+1)*NCOL+j+1] + grid[(i-1)*NCOL+j+1] + grid[(i+1)*NCOL+(j-1)] + grid[(i-1)*NCOL+(j-1)] 
                  - 8*grid[i*NCOL+j]);
       grid_aux[i*NCOL+j] = T;
      Tfull_loc += T;
    }
  }
  if(pid==0){ 
    i=des+1;
    for (j=1; j<NCOL-1; j++)
      {
        T = grid[i*NCOL+j] + 
            0.10 * (grid[(i+1)*NCOL+j]   + grid[(i-1)*NCOL+j]   + grid[i*NCOL+(j+1)]     + grid[i*NCOL+(j-1)] + 
                    grid[(i+1)*NCOL+j+1] + grid[(i-1)*NCOL+j+1] + grid[(i+1)*NCOL+(j-1)] + grid[(i-1)*NCOL+(j-1)] 
                    - 8*grid[i*NCOL+j]);

        grid_aux[i*NCOL+j] = T;
        Tfull_loc += T;
      }
    MPI_Wait(&req2,MPI_STATUS_IGNORE);
    i=des+tam;
    for (j=1; j<NCOL-1; j++)
      {
        T = grid[i*NCOL+j] + 
            0.10 * (grid[(i+1)*NCOL+j]   + grid[(i-1)*NCOL+j]   + grid[i*NCOL+(j+1)]     + grid[i*NCOL+(j-1)] + 
                    grid[(i+1)*NCOL+j+1] + grid[(i-1)*NCOL+j+1] + grid[(i+1)*NCOL+(j-1)] + grid[(i-1)*NCOL+(j-1)] 
                    - 8*grid[i*NCOL+j]);

        grid_aux[i*NCOL+j] = T;
        Tfull_loc += T;
      }
  }
  else if(pid==npr-1){ 
    //MPI_Wait(&req2,MPI_STATUS_IGNORE);
    i=des+tam;
    for (j=1; j<NCOL-1; j++)
      {
        T = grid[i*NCOL+j] + 
            0.10 * (grid[(i+1)*NCOL+j]   + grid[(i-1)*NCOL+j]   + grid[i*NCOL+(j+1)]     + grid[i*NCOL+(j-1)] + 
                    grid[(i+1)*NCOL+j+1] + grid[(i-1)*NCOL+j+1] + grid[(i+1)*NCOL+(j-1)] + grid[(i-1)*NCOL+(j-1)] 
                    - 8*grid[i*NCOL+j]);

        grid_aux[i*NCOL+j] = T;
        Tfull_loc += T;
      }
    MPI_Wait(&req1,MPI_STATUS_IGNORE);
    i=des+1;
    for (j=1; j<NCOL-1; j++)
      {
        T = grid[i*NCOL+j] + 
            0.10 * (grid[(i+1)*NCOL+j]   + grid[(i-1)*NCOL+j]   + grid[i*NCOL+(j+1)]     + grid[i*NCOL+(j-1)] + 
                    grid[(i+1)*NCOL+j+1] + grid[(i-1)*NCOL+j+1] + grid[(i+1)*NCOL+(j-1)] + grid[(i-1)*NCOL+(j-1)] 
                    - 8*grid[i*NCOL+j]);

        grid_aux[i*NCOL+j] = T;
        Tfull_loc += T;
      }
    
  }
  else{
    MPI_Wait(&req1,MPI_STATUS_IGNORE);
    i=des+1;
    for (j=1; j<NCOL-1; j++)
      {
        T = grid[i*NCOL+j] + 
            0.10 * (grid[(i+1)*NCOL+j]   + grid[(i-1)*NCOL+j]   + grid[i*NCOL+(j+1)]     + grid[i*NCOL+(j-1)] + 
                    grid[(i+1)*NCOL+j+1] + grid[(i-1)*NCOL+j+1] + grid[(i+1)*NCOL+(j-1)] + grid[(i-1)*NCOL+(j-1)] 
                    - 8*grid[i*NCOL+j]);

        grid_aux[i*NCOL+j] = T;
        Tfull_loc += T;
      }
    MPI_Wait(&req2,MPI_STATUS_IGNORE);
    i=des+tam;
    for (j=1; j<NCOL-1; j++)
      {
        T = grid[i*NCOL+j] + 
            0.10 * (grid[(i+1)*NCOL+j]   + grid[(i-1)*NCOL+j]   + grid[i*NCOL+(j+1)]     + grid[i*NCOL+(j-1)] + 
                    grid[(i+1)*NCOL+j+1] + grid[(i-1)*NCOL+j+1] + grid[(i+1)*NCOL+(j-1)] + grid[(i-1)*NCOL+(j-1)] 
                    - 8*grid[i*NCOL+j]);

        grid_aux[i*NCOL+j] = T;
        Tfull_loc += T;
      }
  }
  //new values for the grid
  for (i=des+1; i<des+tam+1; i++)//cambiar los indices
  for (j=1; j<NCOL-1; j++)
    grid[i*NCOL+j] = grid_aux[i*NCOL+j]; 
  
  //Allreduce de Tfull (suma)
  MPI_Allreduce(&Tfull_loc,&Tfull,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  return (Tfull);
}

/************************************************************************************/
double calculate_Tmean (struct info_param param, float *grid, float *grid_chips, float *grid_aux, int tam, int des, int pid, int npr) //En principio no veo necesario tocar nada para la version paralela
{
  int    i, j, end, niter;
  double  Tfull;
  double Tmean, Tmean0 = param.t_ext;

  end = 0; niter = 0;

  while (end == 0)
  {
    niter++;
    Tmean = 0.0;

    // heat injection and air cooling 
    thermal_update (param, grid, grid_chips, tam, des, pid, npr);
    
    // thermal diffusion
    Tfull = thermal_diffusion(param, grid, grid_aux, tam, des, pid, npr);
  
    // convergence every 10 iterations
    if (niter % 10 == 0)
    {
      Tmean = Tfull / ((NCOL-2)*(NROW-2));
      if ((fabs(Tmean - Tmean0) < param.t_delta) || (niter > param.max_iter))
           end = 1;
      else Tmean0 = Tmean;
    }
  } // end while 
  if(pid==0){
    printf ("Iter (par): %d\t", niter);
  }
  return (Tmean);
}

