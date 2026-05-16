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
#include "faux.h"
#include "diffusion.h"

/************************************************************************************/
void init_grid_chips (int conf, struct info_param param, struct info_chips *chips, int **chip_coord, float *grid_chips) //no paralelizar, uno inicializa y reparte
{
  int i, j, n;

  for (i=0; i<NROW; i++)
  for (j=0; j<NCOL; j++)
    grid_chips[i*NCOL+j] = param.t_ext;

  for (n=0; n<param.nchip; n++)
  for (i = chip_coord[conf][2*n]   * param.scale; i < (chip_coord[conf][2*n] + chips[n].h) * param.scale; i++)
  for (j = chip_coord[conf][2*n+1] * param.scale; j < (chip_coord[conf][2*n+1]+chips[n].w) * param.scale; j++)
    grid_chips[(i+1)*NCOL+(j+1)] = chips[n].tchip;
}

/************************************************************************************/
void init_grids_parallel (struct info_param param, float *grid, float *grid_aux, int tam_loc) //paralelizar, cada uno inicializa su parte
{
  int i,j;
  for (i = 0; i<tam_loc; i++)
  for (j=0; j<NCOL; j++)
    grid[i*NCOL+j] = grid_aux[i*NCOL+j] = param.t_ext;
}
/************************************************************************************/












/************************************************************************************/
int main (int argc, char *argv[])
{
    struct info_param param;
    struct info_chips *chips;
    int     **chip_coord;

    float *grid, *grid_chips, *grid_aux;
    struct info_results BT;

    int    conf, i;
    double t0, t1;
    double tej, Tmean;

    const int N_VECINOS = 2;
    int pid, npr;
    MPI_Status info;
    int tam_loc;
    int base, resto;
    int * tam, * dist;
    float * grid_local, * grid_chips_local, * grid_aux_local;

    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &npr);
    MPI_Comm_rank (MPI_COMM_WORLD, &pid);

    if(pid==0){
        // reading initial data file
        if (argc != 2) {
            printf ("\n\nERROR: needs a card description file \n\n");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }

        read_data (argv[1], &param, &chips, &chip_coord);

        printf ("\n  ===================================================================");
        printf ("\n    Thermal diffusion - parallel version ");
        printf ("\n    %d x %d points, %d chips", RSIZE*param.scale, CSIZE*param.scale, param.nchip);
        printf ("\n    T_ext = %1.1f, Tmax_chip = %1.1f, T_delta: %1.3f, Max_iter: %d", param.t_ext, param.tmax_chip, param.t_delta, param.max_iter);
        printf ("\n  ===================================================================\n\n");
    }

    // enviamos los datos cargados a todos los procesos
    MPI_Bcast(&param, sizeof(param), MPI_CHAR, 0, MPI_COMM_WORLD);

    t0 = MPI_Wtime ();

    // calculo del numero de filas que tiene que procesar cada procesador
    base =  (NROW-2)/npr;
    resto = (NROW-2) % npr;

    tam_loc = base;
    if (pid<resto) tam_loc++;

    grid_local = malloc((tam_loc + N_VECINOS) * NCOL * sizeof(float));
    grid_chips_local = malloc((tam_loc + N_VECINOS) * NCOL * sizeof(float));
    grid_aux_local = malloc((tam_loc + N_VECINOS) * NCOL * sizeof(float));

    if (pid == 0) {
        tam = malloc(npr * sizeof(int));
        dist = malloc(npr * sizeof(int));
        for (i=0; i<npr;  i++){
            tam[i] = base * NCOL;
            if (i < resto ) tam[i] += NCOL;

            if (i==0){
                dist[i] = 1 * NCOL;
            } else {
                dist[i] = dist[i-1] + tam[i-1];
            }
        }

        grid = malloc(NROW*NCOL * sizeof(float));
        grid_chips = malloc(NROW*NCOL * sizeof(float));
        grid_aux = malloc(NROW*NCOL * sizeof(float));

        BT.bgrid = malloc(NROW*NCOL * sizeof(float));
        BT.cgrid = malloc(NROW*NCOL * sizeof(float));
        BT.Tmean = MAXDOUBLE;
    }

    for (conf=0; conf<param.nconf; conf++)
    {
        if (pid==0){
            init_grid_chips(conf, param, chips, chip_coord, grid_chips);
        }

        init_grids_parallel (param, grid_local, grid_aux_local, tam_loc + N_VECINOS);

        for (i = 0; i < NCOL; i++) {
            grid_chips_local[0 * NCOL + i] = param.t_ext;
            grid_chips_local[(tam_loc + 1) * NCOL + i] = param.t_ext;
        }

        MPI_Scatterv(grid_chips, tam, dist, MPI_FLOAT,
                     &grid_chips_local[1 * NCOL], tam_loc * NCOL, MPI_FLOAT, 0, MPI_COMM_WORLD);

        Tmean = calculate_Tmean (param, grid_local, grid_chips_local, grid_aux_local, tam_loc, pid, npr, info);

        MPI_Gatherv(&grid_local[1 * NCOL], tam_loc * NCOL, MPI_FLOAT, grid, tam,
                    dist, MPI_FLOAT, 0, MPI_COMM_WORLD);

        // Guardar mejor resultado
        if (pid == 0) {
            printf("  Config: %2d    Tmean: %1.2f\n", conf + 1, Tmean);
            results_conf(conf, Tmean, param, grid, grid_chips, &BT);
        }
    }
    t1 = MPI_Wtime ();

    if(pid==0){
        tej = t1 - t0;
        printf ("\n\n >>> Best configuration: %2d    Tmean: %1.2f\n", BT.conf + 1, BT.Tmean);
        printf ("   > Time (parallel): %1.3f s \n\n", tej);

        // writing best configuration results
        results (param, &BT, argv[1]);

        free (grid); free (grid_chips); free (grid_aux);
        free (BT.bgrid); free (BT.cgrid);
        free (tam);
        free (dist);
        for (i=0; i<param.nconf; i++) free (chip_coord[i]);
        free (chips);
        free (chip_coord);
    }

    free (grid_local); free (grid_chips_local); free (grid_aux_local);
    MPI_Finalize ();
    return (0);
}
