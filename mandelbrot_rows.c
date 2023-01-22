#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <math.h>
#include <memory.h>
#include <time.h>

void compute_mandelbrot(unsigned char* color, int x, int y, int xMax, int yMax);

int main(int argc, char** argv){
    /* screen ( integer) coordinate */
	int iXmax = 8000; // default
	int iYmax = 8000; // default

	/* color component ( R or G or B) is coded from 0 to 255 */
	/* it is 24 bit color RGB file */
	const int MaxColorComponentValue = 255; 
	FILE * fp;
	char *filename = "Mandelbrot_rows.ppm";
    /* comment should start with # */
	char *comment = "# ";	

	/* Clock information */
	double start_all, end_all, start_comm, end_comm;
	double cpu_time_used, comm_time;
    
    MPI_Init(&argc, &argv);

    /* world rand and world size */ 
    int world_rank;
    int world_size;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // RGB color array
    unsigned char* color;
    color = (unsigned char*)malloc(3 * (int)sizeof(unsigned char));
    // local colors
    unsigned char* local_colors;
    // global colors array
    unsigned char* global_colors;

    int rows_per_processor;
    int rows_remained;
    rows_per_processor = iYmax / world_size;
    rows_remained = iYmax % world_size;

    // root node create a binary file
    if (world_rank == 0){
		/*create new file,give it a name and open it in binary mode  */
		fp = fopen(filename, "wb"); /* b -  binary mode */

		/*write ASCII header to the file (PPM file format)*/
		fprintf(fp,"P6\n %s\n %d\n %d\n %d\n", comment, iXmax, iYmax, MaxColorComponentValue);

		printf("File: %s successful ly opened for writing.\n", filename);
		printf("Computing Mandelbrot Set. Please wait...\n");

        // Start the computation time (which covers the communication time)
		start_all = MPI_Wtime();
    }

    // allocate memory to local_colors
    if (rows_remained != 0){
        local_colors = (unsigned char*)malloc((rows_per_processor+1)*iXmax*(int)sizeof(unsigned char)*3);
    } else {
        local_colors = (unsigned char*)malloc(rows_per_processor*iXmax*(int)sizeof(unsigned char)*3);
    }

    int offset = 0;
    // each processor computes mandebrot set
    for (int j = world_rank; j < iYmax; j += world_size){
            for (int i = 0; i < iXmax; i++){
                compute_mandelbrot(color, i, j, iXmax, iYmax);
                memcpy(local_colors + offset, color, 3);
                offset += 3;
            }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    // root node create an array to store all color values
    if (world_rank == 0) {
        // allocate memory to global_colors
        if (rows_remained != 0){
            global_colors = (unsigned char*)malloc((rows_per_processor+1) * world_size * iXmax * 3 *(int)sizeof(unsigned char));
        } else {
            global_colors = (unsigned char*)malloc(iXmax * iYmax * 3 * (int)sizeof(unsigned char));
        }  
        start_comm = MPI_Wtime();       
    }

    int localOffset = 0;
    int pointOffset = 0;
    // root node gather color values from all processors
    if (rows_remained != 0){
        for (int i = 0; i < (rows_per_processor+1); i++) {
            MPI_Gather(local_colors+localOffset, iXmax*3, MPI_UNSIGNED_CHAR, global_colors+pointOffset, iXmax*3, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
            localOffset += iXmax*3;
            pointOffset += iXmax*3*world_size;
        }
    } else {
        for (int i = 0; i < rows_per_processor; i++) {
            MPI_Gather(local_colors+localOffset, iXmax*3, MPI_UNSIGNED_CHAR, global_colors+pointOffset, iXmax*3, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
            localOffset += iXmax*3;
            pointOffset += iXmax*3*world_size;
        }
    }

    //root node write file
    if (world_rank == 0){
        end_comm = MPI_Wtime();
        for (int k = 0; k < iXmax * iYmax * 3; k +=3){
            color[0]=global_colors[k];
            color[1]=global_colors[k+1];
            color[2]=global_colors[k+2];
            /* write color to the file */
			fwrite(color, 1, 3, fp);
        }
		end_all = MPI_Wtime();
        cpu_time_used = end_all-start_all;
        comm_time = end_comm-start_comm; 
        printf("Mandelbrot communication time: %lf\n", comm_time);
		printf("Mandelbrot computational process time: %lf\n", cpu_time_used);
        free(global_colors);

	}
    free(color);
    free(local_colors);
    MPI_Finalize();
    return (0);
}

/* This function compute the mandelbrot set based on given values and store the corresponding corlor value*/
void compute_mandelbrot(unsigned char* color, int x, int y, int xMax, int yMax) {
    /* screen ( integer) coordinate */
	int iXmax = xMax;
	int iYmax = yMax;

    /* world ( double) coordinate = parameter plane*/
	double Cx, Cy;
	const double CxMin = -2.5;
	const double CxMax = 1.5;
	const double CyMin = -2.0;
	const double CyMax = 2.0;

    /* each pixel width and height*/
	double PixelWidth = (CxMax - CxMin)/iXmax;
	double PixelHeight = (CyMax - CyMin)/iYmax;

	/* color component ( R or G or B) is coded from 0 to 255 */
	/* it is 24 bit color RGB file */
	const int MaxColorComponentValue = 255; 

	/* Z = Zx + Zy*i;	Z0 = 0 */
    /* Zx2 = Zx*Zx;  Zy2 = Zy*Zy  */
	double Zx, Zy;
	double Zx2, Zy2; 

	int Iteration;

    // default 
	const int IterationMax = 2000; 

    /* bail-out value , radius of circle ;  */
	const double EscapeRadius = 400;
	double ER2 = EscapeRadius * EscapeRadius;
	Cy = CyMin + (y * PixelHeight);
    
    if (fabs(Cy) < (PixelHeight / 2)){
        Cy = 0.0; /* Main antenna */
	}

    Cx = CxMin + (x * PixelWidth);
    /* initial value of orbit = critical point Z= 0 */
    Zx = 0.0;
    Zy = 0.0;
    Zx2 = Zx * Zx;
    Zy2 = Zy * Zy;

    for(Iteration = 0; Iteration < IterationMax && ((Zx2 + Zy2) < ER2); Iteration++){
        Zy = (2 * Zx * Zy) + Cy;
        Zx = Zx2 - Zy2 + Cx;
        Zx2 = Zx * Zx;
        Zy2 = Zy * Zy;
	};

    /* compute  pixel color (24 bit = 3 bytes) */
    if (Iteration == IterationMax){
        // Point within the set. Mark it as black
        color[0] = 0;
        color[1] = 0;
        color[2] = 0;
    } else{
        // Point outside the set. Mark it as white
        double c = 3*log((double)Iteration)/log((double)(IterationMax) - 1.0);
        if (c < 1){
            color[0] = 0;
            color[1] = 0;
            color[2] = 255*c;
        } else if (c < 2){
            color[0] = 0;
            color[1] = 255*(c-1);
            color[2] = 255;
        } else{
            color[0] = 255*(c-2);
            color[1] = 255;
            color[2] = 255;
        }
    }
}
