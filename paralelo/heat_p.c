/* heat_s.c

     Difusion del calor en 2 dimensiones      Version en serie

     Se analizan las posiciones de los chips en una tarjeta, para conseguir la temperatura minima
     de la tarjeta. Se utiliza el metodo de tipo Poisson, y la tarjeta se discretiza en una rejilla
     de puntos 2D.

     Entrada: card > la definicion de la tarjeta y las configuraciones a simular
     Salida: la mejor configuracion y la temperatura media
	      card_s.chips: situacion termica inicial
        card_s.res: la situacion termica final

     defines.h: definiciones de ciertas variables y estructuras de datos

     Compilar con estos dos ficheros:
       diffusion.c: insertar calor, difundir y calcular la temperatura media hasta que se estabilice
       faux.c: ciertas funciones auxiliares

************************************************************************************************/

#include <stdio.h>
#include <values.h>
#include <time.h>
#include <mpi.h>

#include "defines.h"
#include "faux_p.h"
#include "diffusion_p.h"

/************************************************************************************/
void init_grid_chips (int conf, struct info_param param, struct info_chips *chips, int **chip_coord, float *grid_chips) //no paralelizar, uno inicializa y reparte
{
  int i, j, n;

  for (i=0; i<NROW; i++)
  for (j=0; j<NCOL; j++)
    grid_chips[i*NCOL+j] = param.t_ext;

  for (n=0; n<param.nchip; n++)
  for (i = chip_coord[conf][2*n]   * param.scale; i- < (chip_coord[conf][2*n] + chips[n].h) * param.scale; i++)
  for (j = chip_coord[conf][2*n+1] * param.scale; j < (chip_coord[conf][2*n+1]+chips[n].w) * param.scale; j++)
    grid_chips[(i+1)*NCOL+(j+1)] = chips[n].tchip;
}

/************************************************************************************/
void init_grids_parallel (struct info_param param, float *grid, float *grid_aux, int tam) //paralelizar, cada uno inicializa su parte
{
  int i,j;
  for (i = 0; i<tam; i++)
  for (j=0; j<NCOL; j++)
    grid[i*NCOL+j] = grid_aux[i*NCOL+j] = param.t_ext;
}
/************************************************************************************/












/************************************************************************************/
int main (int argc, char *argv[])
{
  struct info_param param;
  struct info_chips *chips;
  int	 **chip_coord;

  float *grid, *grid_chips, *grid_aux;
  struct info_results BT;

  int    conf, i;
  double t0, t1;
  double tej, Tmean;
  int pid, npr;







  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &npr);
  MPI_Comm_rank (MPI_COMM_WORLD, &pid);

  if(pid==0){
      // reading initial data file
    if (argc != 2) {
        printf ("\n\nERROR: needs a card description file \n\n");
        MPI_Abort(MPI_COMM_WORLD, -1); //para abortar la ejecucion en caso de que haya un error en los parametros
    }

  read_data (argv[1], &param, &chips, &chip_coord);

  if(pid==0){
    printf ("\n  ===================================================================");
    printf ("\n    Thermal diffusion - parallel version ");
    printf ("\n    %d x %d points, %d chips", RSIZE*param.scale, CSIZE*param.scale, param.nchip);
    printf ("\n    T_ext = %1.1f, Tmax_chip = %1.1f, T_delta: %1.3f, Max_iter: %d", param.t_ext, param.tmax_chip, param.t_delta, param.max_iter);
    printf ("\n  ===================================================================\n\n");
}

  //enviamos los datos cargados, cargar los datos tarda tiempo, asi nos lo ahorramos para el resto de procesos.
  MPI_Bcast(&param, sizeof(param), MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(&chips, sizeof(chips), MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(&chip_coord, sizeof(chip_coord), MPI_CHAR, 0, MPI_COMM_WORLD);

//clock_gettime (CLOCK_REALTIME, &t0);
t0 = MPI_Wtime ();

//calculo del numero de filas que tiene que procesar cada procesador
    int tam_loc,tmean_loc;
    int * tam = malloc(npr * sizeof(int));
    int * dist = malloc(npr * sizeof(int));
    int base =  (NROW-2)/npr;
    int resto = (NROW-2) % npr;

    tam_loc = base;
    if (pid<resto) tam_loc++;

    const int N_VECINOS = 2;

    float * grid_local = malloc((base + N_VECINOS) * NCOL * sizeof(float));
    float * grid_chips_local = malloc((base + N_VECINOS) * NCOL * sizeof(float));
    float * grid_aux_local = malloc((base + N_VECINOS) * NCOL * sizeof(float));

    if (pid == 0) {
        for (i=0; i<npr;  i++){
        tam[i] = base;
        if (i < resto ) tam[i]++;
        if (i==0) dis[i] = 0;
        else dis[i] = dis[i-1]+tam[i-1];
        }

        //TODO pid 0 tiene que mandar los datos de grid_chips a los demas procesos, tiene que separar todo en filas y mandar las respectivas filas a sus procesos.

        grid = malloc(NROW*NCOL * sizeof(float));
        grid_chips = malloc(NROW*NCOL * sizeof(float));
        grid_aux = malloc(NROW*NCOL * sizeof(float));
        for (conf=0; conf<param.nconf; conf++)
        {
            init_grid_chips(conf, param, grid_chips, chip_coord, grid_chips);
            //mando los datos de cada fila de chips a su respectivo proceso
            MPI_Scatterv(&grid_chips[dis[pid]*NCOL], tam[pid]*NCOL, MPI_FLOAT, &grid_chips_local[0][0], tam[pid]*NCOL, MPI_FLOAT, 0, MPI_COMM_WORLD, &info);

        }
    }
    else {
        for (conf=0; conf<param.nconf; conf++)
        {
            //inicializa sus matrices locales
            init_grids (param, grid_local, grid_aux_local);
            //recibe los datos de los chips desde el proceso 0
            MPI_Recv (&grid_chips_local[0][0], tam[pid]*NCOL, MPI_FLOAT, 0, 0,MPI_COMM_WORLD, &info);
            //procesa la simulacion localmente

        }
    }




    BT.bgrid = malloc(NROW*NCOL * sizeof(float));
    BT.cgrid = malloc(NROW*NCOL * sizeof(float));
    BT.Tmean = MAXDOUBLE;

    // loop to process chip configurations
    for (conf=0; conf<param.nconf; conf++)
    {
        // inintial values for grids
        init_grid_chips (conf, param, chips, chip_coord, grid_chips);
        init_grids (param, grid, grid_aux);

        // main loop: thermal injection/disipation until convergence (t_delta or max_iter)
        Tmean = calculate_Tmean (param, grid, grid_chips, grid_aux);
        printf ("  Config: %2d    Tmean: %1.2f\n", conf + 1, Tmean);

        // processing configuration results
        results_conf (conf, Tmean, param, grid, grid_chips, &BT);
    }



    t1 = MPI_Wtime ();
  //clock_gettime (CLOCK_REALTIME, &t1);

  tej = (t1 - t0) + (t1 - t0)/(double)1e9;
  printf ("\n\n >>> Best configuration: %2d    Tmean: %1.2f\n", BT.conf + 1, BT.Tmean);
  //printf ("   > Time (parallel): %1.3f s \n\n", tej);
  // writing best configuration results
  results (param, &BT, argv[1]);


  free (grid);free (grid_chips);free (grid_aux);
  free (BT.bgrid);free (BT.cgrid);
  free (chips);
  for (i=0; i<param.nconf; i++) free (chip_coord[i]);
  free (chip_coord);

  return (0);
}
