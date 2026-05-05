/* File: diffusion_s.c */

#include "defines.h"
#include <mpi.h>


/************************************************************************************/
void thermal_update (struct info_param param, float *grid, float *grid_chips, int tam_loc)
{
    int i, j, a, b;

      // Inyección de calor en las posiciones de los chips
      // i<=tam_loc en lugar de i<NROW-1
      for (i=1; i<=tam_loc; i++) {
          for (j=1; j<NCOL-1; j++) {
              if (grid_chips[i*NCOL+j] > grid[i*NCOL+j]) {
                  grid[i*NCOL+j] += 0.05 * (grid_chips[i*NCOL+j] - grid[i*NCOL+j]);
              }
          }
      }

      // Enfriamiento en el centro de la tarjeta
      a = 0.45*(NCOL-2) + 1;
      b = 0.55*(NCOL-2) + 1;

      for (i=1; i<=tam_loc; i++) {
          for (j=a; j<b; j++) {
              grid[i*NCOL+j] -= 0.01 * (grid[i*NCOL+j] - param.t_ext);
          }
      }
    }

/************************************************************************************/
double thermal_diffusion (struct info_param param, float *grid, float *grid_aux, int tam_loc)
{
  int    i, j;
  double  T;
  double Tfull_local = 0.0;

  for (i=1; i<=tam_loc; i++)
    for (j=1; j<NCOL-1; j++)
    {
      T = grid[i*NCOL+j] +
          0.10 * (grid[(i+1)*NCOL+j]   + grid[(i-1)*NCOL+j]   + grid[i*NCOL+(j+1)]     + grid[i*NCOL+(j-1)] +
                  grid[(i+1)*NCOL+j+1] + grid[(i-1)*NCOL+j+1] + grid[(i+1)*NCOL+(j-1)] + grid[(i-1)*NCOL+(j-1)]
                  - 8*grid[i*NCOL+j]);

      grid_aux[i*NCOL+j] = T;
      Tfull_local += T;;
    }

    //new values for the grid
    for (i=1; i<=tam_loc; i++)
    for (j=1; j<NCOL-1; j++)
      grid[i*NCOL+j] = grid_aux[i*NCOL+j];

    return (Tfull_local);
}

/************************************************************************************/
double calculate_Tmean (struct info_param param, float *grid, float *grid_chips, float *grid_aux, int tam_loc, int pid, int npr, MPI_Status info)
{
  int    i, j, end, niter;
  double Tfull_global;
  double  Tfull_local;
  double Tmean, Tmean0 = param.t_ext;
  end = 0; niter = 0;

  while (end == 0)
  {
    niter++;
    Tmean = 0.0;

    // heat injection and air cooling
    thermal_update (param, grid, grid_chips,tam_loc);

    // 1. envia las filas pares envian su ultima fila a las filas impares las filas impares sus
    //    primeras filas a las pares
    // 2. las pares envian sus primeras filas a las impares
    //    las impares sus ultimas filas a la pares.
    //pares a sus siguientes
    if (pid % 2 == 0) {
        if (pid != npr-1) { // Si el par NO es el último proceso
            // Envío mi última fila (tam_loc) y recibo en mi fantasma inferior (tam_loc + 1)
            MPI_Ssend(&grid_local[tam_loc * NCOL], NCOL, MPI_FLOAT, pid + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(&grid_local[(tam_loc + 1) * NCOL], NCOL, MPI_FLOAT, pid + 1, 1, MPI_COMM_WORLD, &info);
        }
    }
    //impares a sus anteriores
    else {
        // Los impares reciben del par de arriba y luego le contestan
        MPI_Recv(&grid_local[0 * NCOL], NCOL, MPI_FLOAT, pid - 1, 0, MPI_COMM_WORLD, &info);
        MPI_Ssend(&grid_local[1 * NCOL], NCOL, MPI_FLOAT, pid - 1, 1, MPI_COMM_WORLD);
    }
    //impares a sus siguientes
    if (pid % 2 != 0) {
        if (pid != npr-1) { // Si el impar NO es el último proceso
            // Envío mi última fila (tam_loc) y recibo en mi inferior (tam_loc + 1)
            MPI_Ssend(&grid_local[tam_loc * NCOL], NCOL, MPI_FLOAT, pid + 1, 2, MPI_COMM_WORLD);
            MPI_Recv(&grid_local[(tam_loc + 1) * NCOL], NCOL, MPI_FLOAT, pid + 1, 3, MPI_COMM_WORLD, &info);
        }
    }
    //pares a sus anteriores
    else {
        if (pid > 0) { // Si el par NO es el proceso 0 entonces envia a su anterior
            // Los pares ahora miran hacia arriba, reciben del impar y le contestan
            MPI_Recv(&grid_local[0 * NCOL], NCOL, MPI_FLOAT, pid - 1, 2, MPI_COMM_WORLD, &info);
            MPI_Ssend(&grid_local[1 * NCOL], NCOL, MPI_FLOAT, pid - 1, 3, MPI_COMM_WORLD);
        }
    }

    // thermal diffusion
    Tfull_local = thermal_diffusion(param, grid, grid_aux,tam_loc);


    // convergence every 10 iterations
    if (niter % 10 == 0)
    {
        //juntamos la suma de calor de todos los procesos para saber si hay que parar
        MPI_Allreduce(&Tfull_local, &Tfull_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        Tmean = Tfull_global / ((NCOL-2)*(NROW-2));
        if ((fabs(Tmean - Tmean0) < param.t_delta) || (niter > param.max_iter))
            end = 1;
        else Tmean0 = Tmean;
    }
} // end while
  if (pid == 0) printf ("Iter (par): %d\t", niter);
  return (Tmean);
}
