#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <memory.h>

// Main program
int main()
 {
	/* screen ( integer) coordinate */
	int iX,iY;
	const int iXmax = 8000; // default
	const int iYmax = 8000; // default

	/* world ( double) coordinate = parameter plane*/
	double Cx, Cy;
	const double CxMin = -2.5;
	const double CxMax = 1.5;
	const double CyMin = -2.0;
	const double CyMax = 2.0;

	/* */
	double PixelWidth = (CxMax - CxMin)/iXmax;
	double PixelHeight = (CyMax - CyMin)/iYmax;

	/* color component ( R or G or B) is coded from 0 to 255 */
	/* it is 24 bit color RGB file */
	const int MaxColorComponentValue = 255; 
	FILE * fp;
	char *filename = "Mandelbrot.ppm";
	char *comment = "# ";	/* comment should start with # */

	// RGB color array
	static unsigned char color[3];

	/* Z = Zx + Zy*i;	Z0 = 0 */
	double Zx, Zy;
	double Zx2, Zy2; /* Zx2 = Zx*Zx;  Zy2 = Zy*Zy  */
	/*  */
	int Iteration;
	const int IterationMax = 2000; // default

	/* bail-out value , radius of circle ;  */
	const double EscapeRadius = 400;
	double ER2 = EscapeRadius * EscapeRadius;
	
	/* Clock information */
	clock_t start_all, end_all, start_compute, end_compute;
	double cpu_time_used, compute_time;

	/* a offset pointer which assist to store color in points */ 
	int offset = 0;

	// Get current clock time.
	start_all = clock();

	/*create new file,give it a name and open it in binary mode  */
	fp = fopen(filename, "wb"); /* b -  binary mode */

	/* declare a dynamic array to store all points color */
	unsigned char* points;
	points = (unsigned char*)malloc(iXmax * iYmax * 3 * (int)sizeof(unsigned char));

	/*write ASCII header to the file (PPM file format)*/
	fprintf(fp,"P6\n %s\n %d\n %d\n %d\n", comment, iXmax, iYmax, MaxColorComponentValue);

	printf("File: %s successfully opened for writing.\n", filename);
	printf("Computing Mandelbrot Set. Please wait...\n");

	// Get current clock time.
	start_compute = clock();

	/* compute and write image data bytes to the file */
	for(iY = 0; iY < iYmax; iY++)
	{
		Cy = CyMin + (iY * PixelHeight);
		if (fabs(Cy) < (PixelHeight / 2))
		{
			Cy = 0.0; /* Main antenna */
		}

		for(iX = 0; iX < iXmax; iX++)
		{
			Cx = CxMin + (iX * PixelWidth);
			/* initial value of orbit = critical point Z= 0 */
			Zx = 0.0;
			Zy = 0.0;
			Zx2 = Zx * Zx;
			Zy2 = Zy * Zy;

			/* */
			for(Iteration = 0; Iteration < IterationMax && ((Zx2 + Zy2) < ER2); Iteration++)
			{
				Zy = (2 * Zx * Zy) + Cy;
				Zx = Zx2 - Zy2 + Cx;
				Zx2 = Zx * Zx;
				Zy2 = Zy * Zy;
			};

			/* compute  pixel color (24 bit = 3 bytes) */
			if (Iteration == IterationMax)
			{
				// Point within the set. Mark it as black
				color[0] = 0;
				color[1] = 0;
				color[2] = 0;
			}
			else 
			{
				// Point outside the set. Mark it as white
				double c = 3*log((double)Iteration)/log((double)(IterationMax) - 1.0);
				if (c < 1)
				{
					color[0] = 0;
					color[1] = 0;
					color[2] = 255*c;
				}
				else if (c < 2)
				{
					color[0] = 0;
					color[1] = 255*(c-1);
					color[2] = 255;
				}
				else
				{
					color[0] = 255*(c-2);
					color[1] = 255;
					color[2] = 255;
				}
			}

			/* write color to the file */
			//fwrite(color, 1, 3, fp);
			memcpy(points + offset, color, 3);
			offset += 3;
		}
	}
	// Get the clock current time
	// Subtract endCompute from startCompute to get the computation time
	printf("Completed Computing Mandelbrot Set.\n");
	end_compute = clock();
	compute_time = ((double)(end_compute - start_compute)) / CLOCKS_PER_SEC;
	printf("Mandelbrot set computation time: %lf\n", compute_time);

	// Write color to the file
	for (int i = 0; i < iXmax * iYmax * 3; i +=3) {
		color[0]=points[i];
		color[1]=points[i+1];
		color[2]=points[i+2];
		/* write color to the file */
		fwrite(color, 1, 3, fp);	
	}

	// Get the clock current time again
	// Subtract end from start to get the CPU time used.
	start_all = clock();
	cpu_time_used = ((double)(end_all - start_all)) / CLOCKS_PER_SEC;

	fclose(fp);

	printf("Completed Computing Mandelbrot Set.\n");
	printf("File: %s successfully closed.\n", filename);
	printf("Mandelbrot computational process time: %lf\n", cpu_time_used);

	return 0;
 }
