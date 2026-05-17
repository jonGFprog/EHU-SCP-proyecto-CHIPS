/* heat_fase2_p.c

     Difusión del calor en 2 dimensiones - Versión Paralela (Fase 2)
     Arquitectura: Master-Worker Asíncrona Multinivel (Afín al Hardware)

     Entrada: archivo de tarjeta y Tamaño del Grupo (P)
     Salida: la mejor configuración y la temperatura media

********************************************************************************/

#include <stdio.h>
#include <values.h>
#include <time.h>
#include <stdlib.h>
#include <mpi.h>

#include "defines.h"
#include "faux.h"
#include "diffusion.h"

/* ===================================================================== */
/* FUNCIONES AUXILIARES DE INICIALIZACIÓN                                */
/* ===================================================================== */

void init_grid_chips (struct info_param param, struct info_chips *chips, int *chip_coord, float *grid_chips)
{
    int i, j, n;

    for (i=0; i<NROW; i++)
        for (j=0; j<NCOL; j++)
            grid_chips[i*NCOL+j] = param.t_ext;

    for (n=0; n<param.nchip; n++)
        for (i = chip_coord[2*n] * param.scale; i < (chip_coord[2*n] + chips[n].h) * param.scale; i++)
            for (j = chip_coord[2*n+1] * param.scale; j < (chip_coord[2*n+1]+chips[n].w) * param.scale; j++)
                grid_chips[(i+1)*NCOL+(j+1)] = chips[n].tchip;
}

void init_grids_parallel (struct info_param param, float *grid, float *grid_aux, int tam_loc)
{
    int i, j;
    for (i = 0; i<tam_loc; i++)
        for (j=0; j<NCOL; j++)
            grid[i*NCOL+j] = grid_aux[i*NCOL+j] = param.t_ext;
}

/* ===================================================================== */
/* PROGRAMA PRINCIPAL                                                    */
/* ===================================================================== */

int main (int argc, char *argv[])
{
    /* Variables Globales */
    int pid, npr;
    double t0, t1, tej;
    struct info_param param;
    struct info_chips *chips;
    int **chip_coord;
    struct info_results BT;
    const int N_VECINOS = 2;
    int P, ngroups, n_workers;

    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &npr);
    MPI_Comm_rank (MPI_COMM_WORLD, &pid);


    /* ===================================================================== */
    /* 1. LECTURA Y DISTRIBUCIÓN DE DATOS INICIALES                          */
    /* ===================================================================== */


    if (argc != 3) {
        if (pid == 0) printf ("\n\nERROR: Uso: %s <archivo_tarjeta> <tam_grupo_P>\n\n", argv[0]);
        MPI_Finalize();
        return -1;
    }

    P = atoi(argv[2]);
    n_workers = npr - 1;

    if (pid == 0) {
        if (n_workers <= 0 || n_workers % P != 0) {
            printf ("\n\nERROR: Los %d trabajadores no se pueden dividir en grupos exactos de %d procesos.\n", n_workers, P);
            printf ("Fórmula correcta: Procesos totales = 1 Manager + (N_Grupos * %d)\n\n", P);
            MPI_Abort(MPI_COMM_WORLD, -1);
        }

        read_data (argv[1], &param, &chips, &chip_coord);
        ngroups = n_workers / P;

        printf ("\n  ===================================================================");
        printf ("\n    Thermal diffusion - parallel version ");
        printf ("\n    %d x %d points, %d chips", RSIZE*param.scale, CSIZE*param.scale, param.nchip);
        printf ("\n    T_ext = %1.1f, Tmax_chip = %1.1f, T_delta: %1.3f, Max_iter: %d", param.t_ext, param.tmax_chip, param.t_delta, param.max_iter);
        printf ("\n    %d groups of %d workers", ngroups, P);
        printf ("\n  ===================================================================\n\n");
    }

    MPI_Bcast(&param, sizeof(param), MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ngroups, 1, MPI_INT, 0, MPI_COMM_WORLD);
    t0 = MPI_Wtime ();


    /* ===================================================================== */
    /* 2. CREACIÓN DE COMUNICADORES (TOPOLOGÍA AFÍN AL HARDWARE)             */
    /* ===================================================================== */

    MPI_Comm grupo, lideres;
    int color, pid_lideres, pid_grupo, npr_grupo;

    /* Comunicador de Grupos de Trabajo (Reparto por Bloques Contiguos) */
    if (pid == 0) color = MPI_UNDEFINED;
    else color = (pid - 1) / P;

    MPI_Comm_split(MPI_COMM_WORLD, color, 0, &grupo);

    /* Comunicador exclusivo de Líderes (Manager + Jefes de Grupo) */
    if (pid == 0 || (pid - 1) % P == 0) color = 0;
    else color = MPI_UNDEFINED;

    MPI_Comm_split(MPI_COMM_WORLD, color, pid, &lideres);

    /* Asignación de rangos locales */
    if (pid == 0) {
        pid_grupo = -1;
        npr_grupo = -1;
    } else {
        MPI_Comm_rank(grupo, &pid_grupo);
        MPI_Comm_size(grupo, &npr_grupo);
    }

    if (color == 0) {
        MPI_Comm_rank(lideres, &pid_lideres);
    } else {
        pid_lideres = -1;
    }


    /* ===================================================================== */
    /* 3. RESERVA DE MEMORIA Y CÁLCULO DE REPARTO                            */
    /* ===================================================================== */

    int tam_loc, base, resto;
    int *tam, *dist;
    float *grid, *grid_chips;
    float *grid_local, *grid_chips_local, *grid_aux_local;
    int i, k;

    /* Memoria para Trabajadores */
    if (pid != 0) {
        base = (NROW - 2) / npr_grupo;
        resto = (NROW - 2) % npr_grupo;
        tam_loc = base + (pid_grupo < resto ? 1 : 0);

        grid_local = malloc((tam_loc + N_VECINOS) * NCOL * sizeof(float));
        grid_chips_local = malloc((tam_loc + N_VECINOS) * NCOL * sizeof(float));
        grid_aux_local = malloc((tam_loc + N_VECINOS) * NCOL * sizeof(float));

        chips = (struct info_chips *) malloc(param.nchip * sizeof(struct info_chips));
    }

    if (pid_lideres >= 0) {
        MPI_Bcast(chips, param.nchip * sizeof(struct info_chips), MPI_BYTE,
                          0, lideres);
    }

    /* Memoria para Jefes de Grupo */
    if (pid_grupo == 0) {
        tam = malloc(npr_grupo * sizeof(int));
        dist = malloc(npr_grupo * sizeof(int));
        for (i = 0; i < npr_grupo; i++) {
            tam[i] = base * NCOL;
            if (i < resto) tam[i] += NCOL;
            dist[i] = (i == 0) ? 1 * NCOL : dist[i-1] + tam[i-1];
        }

        grid = malloc(NROW * NCOL * sizeof(float));
        grid_chips = malloc(NROW * NCOL * sizeof(float));

        /* Limpieza de bordes físicos para no recolectar "basura" de la RAM */
        for(k = 0; k < NROW * NCOL; k++){
            grid[k] = param.t_ext;
            grid_chips[k] = param.t_ext;
        }

        BT.bgrid = malloc(NROW * NCOL * sizeof(float));
        BT.cgrid = malloc(NROW * NCOL * sizeof(float));
        BT.Tmean = MAXDOUBLE;
    }

    /* Memoria para el Manager */
    if (pid == 0) {
        BT.bgrid = malloc(NROW * NCOL * sizeof(float));
        BT.cgrid = malloc(NROW * NCOL * sizeof(float));
    }


    /* ===================================================================== */
    /* 4. BUCLE PRINCIPAL (ASIGNACIÓN ASÍNCRONA)                             */
    /* ===================================================================== */

    char sol = 's';
    int fin = 1;
    int conf;
    MPI_Status status;

    /* Buffer grande para evitar desbordamiento por overhead de MPI_Pack */
    int tam_buf = sizeof(int) * 2 * param.nchip + sizeof(int);
    char buf[tam_buf];
    int pos = 0;

    if (pid == 0) {
        /* ----------------------- MANAGER ----------------------- */
        for (conf = 0; conf < param.nconf; conf++) {
            MPI_Recv(&sol, 1, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, lideres, &status);

            pos = 0;
            MPI_Pack(&conf, 1, MPI_INT, buf, tam_buf, &pos, lideres);
            MPI_Pack(chip_coord[conf], 2 * param.nchip, MPI_INT, buf, tam_buf, &pos, lideres);
            MPI_Send(buf, tam_buf, MPI_PACKED, status.MPI_SOURCE, 19, lideres);
        }

        /* Señal de fin bloqueante para que no se destruya la variable antes de llegar */
        for (i = 0; i < ngroups; i++) {
            MPI_Send(&fin, 1, MPI_INT, i+1, 12, lideres);
        }

    } else {
        /* ----------------------- WORKERS ----------------------- */
        fin = 0;
        int tag = -1;
        int chip_coord_loc[2 * param.nchip];
        double Tmean;
        MPI_Status info;

        while (!fin) {
            pos = 0;
            if (pid_grupo == 0) {
                MPI_Send(&sol, 1, MPI_CHAR, 0, 0, lideres);
                MPI_Probe(0, MPI_ANY_TAG, lideres, &status);

                if (status.MPI_TAG == 12) {
                    MPI_Recv(&fin, 1, MPI_INT, 0, 12, lideres, &status);
                    tag = 12;
                } else {
                    MPI_Recv(buf, tam_buf, MPI_PACKED, 0, 19, lideres, &status);
                    MPI_Unpack(buf, tam_buf, &pos, &conf, 1, MPI_INT, lideres);
                    MPI_Unpack(buf, tam_buf, &pos, chip_coord_loc, 2 * param.nchip, MPI_INT, lideres);
                    tag = 0;
                }
            }
            MPI_Bcast(&tag, 1, MPI_INT, 0, grupo);

            if (tag == 12) {
                fin = 1;
            } else {
                if (pid_grupo == 0) {
                    init_grid_chips(param, chips, chip_coord_loc, grid_chips);
                }

                init_grids_parallel (param, grid_local, grid_aux_local, tam_loc + N_VECINOS);

                for (i = 0; i < NCOL; i++) {
                    grid_chips_local[0 * NCOL + i] = param.t_ext;
                    grid_chips_local[(tam_loc + 1) * NCOL + i] = param.t_ext;
                }

                MPI_Scatterv(grid_chips, tam, dist, MPI_FLOAT,
                            &grid_chips_local[1 * NCOL], tam_loc * NCOL, MPI_FLOAT, 0, grupo);

                /* Cálculo llamando al subgrupo local (pid_grupo, npr_grupo, grupo) */
                Tmean = calculate_Tmean (param, grid_local, grid_chips_local, grid_aux_local, tam_loc, pid_grupo, npr_grupo, info, grupo);

                MPI_Gatherv(&grid_local[1 * NCOL], tam_loc * NCOL, MPI_FLOAT, grid, tam, dist, MPI_FLOAT, 0, grupo);

                if (pid_grupo == 0) {
                    printf("  Config: %2d    Tmean: %1.2f\n", conf + 1, Tmean);
                    results_conf(conf, Tmean, param, grid, grid_chips, &BT);
                }
            }
        }
    }

    /* ===================================================================== */
    /* 5. FASE DE VOTACIÓN Y RECOLECCIÓN DEL GANADOR FINAL                   */
    /* ===================================================================== */

    int tam_bufBT =  NROW * NCOL * sizeof(float) * 2 + sizeof(double) + sizeof(int);
    char *bufBT = malloc(tam_bufBT); // Uso de memoria dinámica para el megabuffer
    pos = 0;

    if (pid == 0) {
        double Tmeans[ngroups];
        int min = 0;

        for (i = 0; i < ngroups; i++) {
            MPI_Probe(MPI_ANY_SOURCE, 16, lideres, &status);
            MPI_Recv(&Tmeans[status.MPI_SOURCE - 1], 1, MPI_DOUBLE, status.MPI_SOURCE, 16, lideres, &status);
        }

        for (i = 1; i < ngroups; i++) {
            if (Tmeans[min] > Tmeans[i]) {
                min = i;
            }
        }

        for (i = 0; i < ngroups; i++) {
            MPI_Send(&min, 1, MPI_INT, i + 1, 1, lideres);
        }

        MPI_Recv(bufBT, tam_bufBT, MPI_PACKED, min + 1, 2, lideres, &status);
        MPI_Unpack(bufBT, tam_bufBT, &pos, &BT.conf, 1, MPI_INT, lideres);
        MPI_Unpack(bufBT, tam_bufBT, &pos, &BT.Tmean, 1, MPI_DOUBLE, lideres);
        MPI_Unpack(bufBT, tam_bufBT, &pos, BT.bgrid, NROW*NCOL, MPI_FLOAT, lideres);
        MPI_Unpack(bufBT, tam_bufBT, &pos, BT.cgrid, NROW*NCOL, MPI_FLOAT, lideres);

    } else if (pid_lideres > 0) {
        int min;

        MPI_Send(&BT.Tmean, 1, MPI_DOUBLE, 0, 16, lideres);
        MPI_Recv(&min, 1, MPI_INT, 0, 1, lideres, &status);

        if (min + 1 == pid_lideres) {
            MPI_Pack(&BT.conf, 1, MPI_INT, bufBT, tam_bufBT, &pos, lideres);
            MPI_Pack(&BT.Tmean, 1, MPI_DOUBLE, bufBT, tam_bufBT, &pos, lideres);
            MPI_Pack(BT.bgrid, NROW*NCOL, MPI_FLOAT, bufBT, tam_bufBT, &pos, lideres);
            MPI_Pack(BT.cgrid, NROW*NCOL, MPI_FLOAT, bufBT, tam_bufBT, &pos, lideres);
            MPI_Send(bufBT, tam_bufBT, MPI_PACKED, 0, 2, lideres);
        }
    }


    /* ===================================================================== */
    /* 6. FINALIZACIÓN Y LIMPIEZA                                            */
    /* ===================================================================== */

    t1 = MPI_Wtime ();

    if (pid == 0) {
        tej = t1 - t0;
        printf ("\n\n >>> Best configuration: %2d    Tmean: %1.2f\n", BT.conf + 1, BT.Tmean);
        printf ("   > Time (parallel): %1.3f s \n\n", tej);

        results (param, &BT, argv[1]);

        free(BT.bgrid);
        free(BT.cgrid);
        for (i = 0; i < param.nconf; i++) free(chip_coord[i]);
        free(chip_coord);
        free(chips);

    } else {
        free(grid_local);
        free(grid_chips_local);
        free(grid_aux_local);
        free(chips);

        if (pid_grupo == 0) {
            free(grid);
            free(grid_chips);
            free(BT.bgrid);
            free(BT.cgrid);
            free(tam);
            free(dist);
        }
    }

    free(bufBT);
    MPI_Finalize();
    return 0;
}
