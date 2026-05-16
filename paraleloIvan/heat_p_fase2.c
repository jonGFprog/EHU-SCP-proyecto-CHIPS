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

#define TAG_READY 0 // Proceso libre para recibir trabajo
#define TAG_RESULT 1 // Resultados procesados + proceso libre
#define TAG_WORK 2 // Nueva configuracion
#define TAG_STOP 3 // Señal de stop para todos los procesos
#define TAG_BEST 4 //señal de mejor resultado


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
    struct info_results BT;
    const int N_VECINOS = 2;

    int     **chip_coord;
    int * tam, * dist;
    float *grid, *grid_chips, *grid_aux;
    float * grid_local, * grid_chips_local, * grid_aux_local;

    int    conf, i;
    double t0, t1;
    int pid, npr;
    int tam_loc;
    int base, resto;
    double tej, Tmean;
    MPI_Status info;

    //grupos de procesos FASE2
    int local_pid, local_npr,P, n_groups, n_workers;
    int color;
    MPI_Comm   comm_group;

    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &npr);
    MPI_Comm_rank (MPI_COMM_WORLD, &pid);
    if (argc != 3) {
            if (pid == 0) printf ("\n\nERROR: Uso: %s <archivo_tarjeta> <tam_grupo_P>\n\n", argv[0]);
            MPI_Finalize();
            return -1;
    }
    P = atoi(argv[2]);
    n_workers = npr-1;
    n_groups = npr/P;

    if(pid==0){
        if (n_workers <= 0 || n_workers % P != 0) {
            printf("\nERROR: El numero de procesos workers (%d) no es divisible por P (%d)\n", n_workers, P);
            printf("Calculo de los procesos : 1 + N*P procesos (ej: 1 + 3*4 = 13)\n\n");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }


        read_data (argv[1], &param, &chips, &chip_coord);

        printf ("\n  ===================================================================");
        printf ("\n    Thermal diffusion - parallel version ");
        printf ("\n    npr = %d, groups = %d, processors by group = %d", npr, n_groups, P);
        printf ("\n    %d x %d points, %d chips", RSIZE*param.scale, CSIZE*param.scale, param.nchip);
        printf ("\n    T_ext = %1.1f, Tmax_chip = %1.1f, T_delta: %1.3f, Max_iter: %d", param.t_ext, param.tmax_chip, param.t_delta, param.max_iter);
        printf ("\n  ===================================================================\n\n");

        color = MPI_UNDEFINED;
    }
    else {
        color = (pid-1) / P+1; //calculo del color del grupo al que pertenece el proceso
    }
    MPI_Comm_split(MPI_COMM_WORLD, color, pid, &comm_group);
    if (comm_group != MPI_COMM_NULL) {
        // Solo los workers entran aquí
        MPI_Comm_rank(comm_group, &local_pid);
        MPI_Comm_size(comm_group, &local_npr);
    }
    // enviamos los datos cargados a todos los procesos
    MPI_Bcast(&param, sizeof(param), MPI_CHAR, 0, MPI_COMM_WORLD);

    if (pid != 0) {
        chips = malloc(param.nchip * sizeof(struct info_chips));

        chip_coord = malloc(param.nconf * sizeof(int*));
        for (i = 0; i < param.nconf; i++) {
            chip_coord[i] = malloc(2 * param.nchip * sizeof(int));
        }
    }

    MPI_Bcast(chips, param.nchip * sizeof(struct info_chips), MPI_BYTE, 0, MPI_COMM_WORLD);

    for (i = 0; i < param.nconf; i++) {
        // Cada fila tiene 2 coordenadas (x,y) por cada chip
        MPI_Bcast(chip_coord[i], 2 * param.nchip, MPI_INT, 0, MPI_COMM_WORLD);
    }

    t0 = MPI_Wtime ();


    if (pid == 0) {
        grid = malloc(NROW*NCOL * sizeof(float));
        grid_chips = malloc(NROW*NCOL * sizeof(float));

        BT.bgrid = malloc(NROW*NCOL * sizeof(float));
        BT.cgrid = malloc(NROW*NCOL * sizeof(float));
        BT.Tmean = MAXDOUBLE;
    }
    if (pid != 0) {
        // calculo del numero de filas que tiene que procesar cada procesador
        base =  (NROW-2)/local_npr;
        resto = (NROW-2) % local_npr;

        tam_loc = base;
        if (local_pid<resto) tam_loc++;

        grid_local = malloc((tam_loc + N_VECINOS) * NCOL * sizeof(float));
        grid_chips_local = malloc((tam_loc + N_VECINOS) * NCOL * sizeof(float));
        grid_aux_local = malloc((tam_loc + N_VECINOS) * NCOL * sizeof(float));

        if (local_pid == 0) {
            tam = malloc(local_npr * sizeof(int));
            dist = malloc(local_npr * sizeof(int));
            for (i=0; i<local_npr;  i++){
                tam[i] = base * NCOL;
                if (i < resto ) tam[i] += NCOL;

                if (i==0){
                    dist[i] = 1 * NCOL;
                } else {
                    dist[i] = dist[i-1] + tam[i-1];
                }
            }
        }
    }

    if (pid == 0) {
        //TODO mover declaraciones arriba
        conf = 0; //conf actual
        int jefe_grupo, rec_tag, conf_finalizadas = 0;
        int grupos_activos = n_groups;
        int mejor_grupo = -1;
        double tmean_recibido;
        MPI_Status status;
        int tmp_val = -1;
        float mejor_tmean = 100000; //numero alto (temeperatura imposible)
        int conf_terminada;
        int conf_asignada[npr];

        while (grupos_activos > 0) {
            double tmean_recibido;

            MPI_Recv(&tmean_recibido, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            jefe_grupo = status.MPI_SOURCE;
            rec_tag = status.MPI_TAG;
            //init_grid_chips(conf, param, chips, chip_coord, grid_chips);

            //si lo recibido es un dato de resultado, lo procesamos, guardamos el emisor, el dato y la conf y actualizamos el contador
            if (rec_tag == TAG_RESULT){
                conf_finalizadas++;
                conf_terminada = conf_asignada[jefe_grupo];
                printf("  Config: %2d    Tmean: %1.2f\n", conf_terminada + 1, tmean_recibido);

                if (tmean_recibido < BT.Tmean) {
                    mejor_grupo = jefe_grupo;
                    BT.Tmean = tmean_recibido;
                    BT.conf = conf_terminada;
                }
            }
            //si necesitan nueva configuracion y quedan configuraciones las enviamos
            if ((rec_tag == TAG_RESULT || rec_tag == TAG_READY) && conf < param.nconf) {
                conf_asignada[jefe_grupo] = conf;
                MPI_Send(&conf, 1, MPI_INT, jefe_grupo, TAG_WORK, MPI_COMM_WORLD);
                conf++;
            }
            //si ya no quedan configuraciones las detenemos
            else if (conf >= param.nconf){
                MPI_Send(&tmp_val, 1, MPI_INT, jefe_grupo, TAG_STOP, MPI_COMM_WORLD);
                grupos_activos--;
            }
        }
        int orden_ganador = 1;
        int orden_perdedor = 0;

        // Avisamos a TODOS los grupos de si deben mandar datos o detenerse
        for (i = 1; i <= n_groups; i++) {
            int pid_jefe_actual = (i - 1) * P + 1;
            if (pid_jefe_actual == mejor_grupo) {
                MPI_Send(&orden_ganador, 1, MPI_INT, pid_jefe_actual, TAG_BEST, MPI_COMM_WORLD);
            } else {
                MPI_Send(&orden_perdedor, 1, MPI_INT, pid_jefe_actual, TAG_BEST, MPI_COMM_WORLD);
            }
        }
        //Recibimos las matrices bgrid y cgrid para results
        MPI_Recv(BT.bgrid, NROW * NCOL, MPI_FLOAT, mejor_grupo, TAG_BEST, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(BT.cgrid, NROW * NCOL, MPI_FLOAT, mejor_grupo, TAG_BEST, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else{
            int mi_conf;
            double mi_tmean = 0.0;
            int tag_envio = TAG_READY;
            MPI_Status status_worker;

            struct info_results BT_local;
            BT_local.Tmean = 100000;

            if (local_pid == 0) {
                BT_local.bgrid = malloc(NROW * NCOL * sizeof(float));
                BT_local.cgrid = malloc(NROW * NCOL * sizeof(float));
                grid_chips = malloc(NROW * NCOL * sizeof(float));
                grid = malloc(NROW * NCOL * sizeof(float));
            }

            while (1) {
                if (local_pid == 0) {
                    MPI_Send(&mi_tmean, 1, MPI_DOUBLE, 0, tag_envio, MPI_COMM_WORLD);
                    MPI_Recv(&mi_conf, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status_worker);

                    if (status_worker.MPI_TAG == TAG_STOP) {
                        mi_conf = -1;
                    }
                }
                MPI_Bcast(&mi_conf, 1, MPI_INT, 0, comm_group);

                if (mi_conf == -1) {
                    break;
                }

                if (local_pid == 0) {
                    init_grid_chips(mi_conf, param, chips, chip_coord, grid_chips);
                }

                init_grids_parallel (param, grid_local, grid_aux_local, tam_loc + N_VECINOS);

                for (i = 0; i < NCOL; i++) {
                    grid_chips_local[0 * NCOL + i] = param.t_ext;
                    grid_chips_local[(tam_loc + 1) * NCOL + i] = param.t_ext;
                }

                MPI_Scatterv(grid_chips, tam, dist, MPI_FLOAT,
                             &grid_chips_local[1 * NCOL], tam_loc * NCOL, MPI_FLOAT, 0, comm_group);

                mi_tmean = calculate_Tmean (param, grid_local, grid_chips_local, grid_aux_local, tam_loc, local_pid, local_npr, info, comm_group);

                MPI_Gatherv(&grid_local[1 * NCOL], tam_loc * NCOL, MPI_FLOAT,
                            grid, tam, dist, MPI_FLOAT, 0, comm_group);

                if (local_pid == 0) {
                    if (mi_tmean < BT_local.Tmean) {
                        BT_local.Tmean = mi_tmean;
                        BT_local.conf = mi_conf;

                        memcpy(BT_local.bgrid, grid, NROW * NCOL * sizeof(float));
                        memcpy(BT_local.cgrid, grid_chips, NROW * NCOL * sizeof(float));
                    }
                    tag_envio = TAG_RESULT;
                }
            } // Fin del while

            if (local_pid == 0) {
                int orden;
                MPI_Recv(&orden, 1, MPI_INT, 0, TAG_BEST, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                if (orden == 1) {

                    MPI_Send(BT_local.bgrid, NROW * NCOL, MPI_FLOAT, 0, TAG_BEST, MPI_COMM_WORLD);
                    MPI_Send(BT_local.cgrid, NROW * NCOL, MPI_FLOAT, 0, TAG_BEST, MPI_COMM_WORLD);
                }
                free(BT_local.bgrid);
                free(BT_local.cgrid);
                free(grid_chips);
                free(grid);
            }
        }




    t1 = MPI_Wtime ();
















    if (pid == 0) {
            tej = t1 - t0;
            // Escribimos los resultados en el archivo
            results (param, &BT, argv[1]);
            printf ("   > Time (parallel): %1.3f s \n\n", tej);

            free (grid);
            free (grid_chips);
            free (BT.bgrid);
            free (BT.cgrid);

            for (i=0; i<param.nconf; i++) free (chip_coord[i]);
            free (chips);
            free (chip_coord);

        } else {
            free (grid_local);
            free (grid_chips_local);
            free (grid_aux_local);

            for (i=0; i<param.nconf; i++) free (chip_coord[i]);
            free (chips);
            free (chip_coord);

            if (local_pid == 0) {
                free (tam);
                free (dist);
            }
        }

        MPI_Finalize ();
        return (0);
    }
