/*BUCUR SORIN 333CB TEMA3 APD*/
/*main.c*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <complex.h>


int *Mandelbrot(double x_min, double y_min, double resolution, int width, int height, int MAX_STEPS, int rank) {

	int step;

	int *image = (int *)malloc((width * height + 2) * sizeof(int));

	int i = 0;
	double cy = y_min;

	while (i < height) {
		
		double cx = x_min;
		int j = 0;

		while(j < width) {

			double complex z = 0 + 0 * I;
			double complex c = cx + cy * I;
			step = 0;

			while ((cabs(z) < 2) && (step < MAX_STEPS)) {
				z = z*z + c;
				step ++;
			}

			image[i * width + j] = step % 256;
			j++;
			cx += resolution;
		}

		i++;
		cy += resolution;

	}

	int start = rank * height;
	int end = start + height;

	image[width * height] = start;
	image[width * height + 1] = end;

	return image;
}

int *Julia(double x_min, double y_min, double resolution, int width, int height, int MAX_STEPS, double complex c, int rank) {
	
	int step;

	int *image = (int *)malloc((width * height + 2) * sizeof(int));

	double zy = y_min;
	int i = 0;

	while(i < height) {

		double zx = x_min;
		int j = 0;

		while(j < width) {

			double complex z = zx + zy * I;
			step = 0;

			while ((cabs(z) < 2) && (step < MAX_STEPS)) {
				z = z * z + c;
				step ++;
			}

			image[i * width + j] = step % 256;
			j++;
			zx += resolution;
		}

		i++;
		zy += resolution;
	}

	int start = rank * height;
	int end = start + height;

	image[width * height] = start;
	image[width * height + 1] = end;

	return image;
}

void write_file(char *filename, int **image, int width, int height) {

	FILE *f = fopen(filename, "w");

	fprintf(f, "P2\n");
	fprintf(f, "%d %d\n", width, height);
	fprintf(f, "%d\n", 255);

	for (int i = height - 1; i >= 0; i--) {
		for (int j = 0; j < width; j++) {
			fprintf(f, "%d ", image[i][j]);
		}
		fprintf(f,"\n");
	}

	fclose(f);
}

int main(int argc, char *argv[]) {

	if (argc < 3) {
		printf("Programul nu a fost rulat bine");
		return -1;
	}

	int type;
	double x_min, x_max, y_min, y_max;
	double resolution;
	int MAX_STEPS;
	double cx, cy;
	int width, height;
	int rank;
	int nr_tasks;

	int **image;
	int *tmp_image;

	double buffer[20];

	MPI_Init(&argc, &argv);
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nr_tasks);

	if (rank == 0) {

		FILE *f = fopen(argv[1], "r");

		fscanf(f, "%i", &type);
		fscanf(f, "%lf %lf %lf %lf", &x_min, &x_max, &y_min, &y_max);
		fscanf(f, "%lf", &resolution);
		fscanf(f, "%i", &MAX_STEPS);

		if (type == 1) {
			fscanf(f, "%lf %lf", &cx, &cy);
		}

		fclose(f);

		width = floor((x_max-x_min)/resolution);
		height = floor((y_max-y_min)/resolution);

		//height-ul matricei trimisa la fiecare worker
		int compute_height = ceil(height / nr_tasks);

		for (int i = 1; i < nr_tasks; i ++) {
			MPI_Send(&type, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(&x_min, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			MPI_Send(&x_max, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			MPI_Send(&y_min, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			MPI_Send(&y_max, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			MPI_Send(&resolution, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			MPI_Send(&compute_height, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(&width, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(&MAX_STEPS, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			if (type == 1) {
				MPI_Send(&cx, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
				MPI_Send(&cy, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			}
		}
		
		image = (int **)malloc(height * sizeof(int *));
		for (int i = 0; i < height; i ++) {
			image[i] = (int *)malloc(width * sizeof(int *));
		}

		if (type == 0) 
			tmp_image = Mandelbrot(x_min, y_min, resolution, width, compute_height, MAX_STEPS, rank);
		else {
			double complex c = cx + cy * I;
			tmp_image = Julia(x_min, y_min, resolution, width, compute_height, MAX_STEPS, c, rank);
		}

		int start = tmp_image[width * compute_height];
		int end = tmp_image[width * compute_height + 1];
		int color = 0;

		for (int i = start; i < end; i++) {
			for (int j = 0; j < width; j++) {
				image[i][j] = tmp_image[color];
				color++;
			}
		}

		for (int i = 1; i < nr_tasks; i++) {
			//informatia primita de la workeri si scrierea in matrice
			MPI_Recv(tmp_image, compute_height * width + 2, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
			
			start = tmp_image[compute_height * width];
			end = tmp_image[compute_height * width + 1];

			int color = 0;
			
			for (int j = start; j < end; j++) {
				for (int k = 0; k < width; k++) {
					image[j][k] = tmp_image[color];
					color++;
				}
			}
		}
		
		write_file(argv[2], image, width, height);

	}
	//fiecare wroker va procesa partea lui de matrice
	else {
		MPI_Recv(&type, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&x_min, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&x_max, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&y_min, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&y_max, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&resolution, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&height, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&width, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&MAX_STEPS, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		if (type == 1) {
			MPI_Recv(&cx, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(&cy, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
		}
		//de unde incep valorile pt y_min
		y_min += rank * height * resolution;

		if (type == 0) 
			tmp_image = Mandelbrot(x_min, y_min, resolution, width, height, MAX_STEPS, rank);
		else {
			double complex c = cx + cy * I;
			tmp_image = Julia(x_min, y_min, resolution, width, height, MAX_STEPS, c, rank);
		}
		//trimite masterului matricea procesata pt a o afisa
		MPI_Send(tmp_image, width * height + 2, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}

	MPI_Finalize();

	return 0;
}