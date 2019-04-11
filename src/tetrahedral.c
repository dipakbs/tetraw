// Written by Dipak Sanap, email: dipakbs@udel.edu
// This code computes oriental tetrahedral order paramter for water from an input csv file.
// CSV file should have the following format
// atom_index,x_coordinate,y_coordinate,z_coordinate

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//struct for input csv data
struct oxygen_coordinates
{
	unsigned int index; //index of an atom
	//x,y and z coordinates of atom	
	float x;
	float y;
	float z;
};

//Given the data in a line in a an inputfile, process it to put it in oxygen_coordinates struct
struct oxygen_coordinates * line_splitter(struct oxygen_coordinates *data, const char *input)
{
	return (sscanf(input, "%u,%f,%f,%f", &data->index, &data->x, &data->y, &data->z) != 7) 
		? NULL : data;
}
//Distance function for two pints in a struct
float getDistance(struct oxygen_coordinates a, struct oxygen_coordinates b)
{
    float distance;
    distance = sqrt((a.x - b.x) * (a.x - b.x) + (a.y-b.y) *(a.y-b.y) + (a.z - b.z) * (a.z - b.z));
    return distance;
}

float getCosine(struct oxygen_coordinates a, struct oxygen_coordinates c, struct oxygen_coordinates b)
{
	// Compute cosine of the angle between CA andCB
	//                   C
	//                  / \
	//                 A---B                  
	//
	float cosangle,q;
	float ab,ca,cb,a2,b2,c2;
	// Distanc
	ab=getDistance(a,b);
	ca=getDistance(c,a);
	cb=getDistance(c,b);
	a2=ca*ca;
	b2=cb*cb;
	c2=ab*ab;
	// Calculate cosine of angle C

	cosangle=(a2+b2-c2)/(2*ca*cb);
	q=pow((cosangle+(1.0f/3.0f)),2);
	return q;
}
//struct for neighbour distance and their indices
struct nbrs_ind {
	float value;
	int index;
};

// comparision function for sorting the neighbours -> qsort library
int cmp(const void *pa, const void *pb)
{
	struct nbrs_ind *pa1 = (struct nbrs_ind *)pa;
	struct nbrs_ind *pa2 = (struct nbrs_ind *)pb;
	if ((*pa1).value < (*pa2).value)
		return -1;
	else if ((*pa1).value > (*pa2).value)
		return 1;
	else
		return 0;
}

//Retrun to the six angle combination from four nearest neighbours.
void angle_comb( int input_array[],int combinations[] ,int central_atom )
{
	int input_natom=4;
        int i,j,k=0; //loop variablesi
        for (i = 0; i < input_natom-1; i++) {
                for (j = i+1; j < input_natom; j++) {
//                      printf("%d\t%d\t%d\n",input_array[i],central_atom,input_array[j]);
                        combinations[k] = input_array[i] ;
                        combinations[k+1] = central_atom ;
                        combinations[k+2] = input_array[j];
                        if(j >= 0) {k = k+3; }
                }
        }     
}

//int calculateDistanceMatrix(float **distMatrix, oxygen_coordinates *coords, int Natoms) {
//}



//main program
int  main(int argc, char *argv[])
{
	FILE *stream; // file pointer
	char *line = NULL; //line pointer
	size_t len = 0;
	ssize_t nread;
	struct oxygen_coordinates * atom_data = NULL; //pointer to oxygen_coordinate struct
	int numatoms = 0; // counter variable for number of atoms

	int i,j,k,p ; // loop initilizers
	//Check for correct number of arguments
	if (argc !=2 ) {
		fprintf(stderr, "Usage: %s <inputfile> <outputfile>\n", argv[0]);
		exit(EXIT_FAILURE);
	}
	//Open the input csv file given in the first argument
	stream = fopen(argv[1], "r");
	if (stream == NULL) {
		perror("fopen");
		exit(EXIT_FAILURE);
	}
	while ((nread = getline(&line, &len, stream)) != -1) {
		if ((atom_data = realloc(atom_data, (size_t) (numatoms + 1) * sizeof(struct oxygen_coordinates))) == NULL) {
			fprintf(stderr, "error not enough memory");
			exit(EXIT_FAILURE);
		}
		line_splitter(&atom_data[numatoms], line);
		numatoms = numatoms + 1;
	}
	free(line);
	fclose(stream);

	// All the data is read in memory in atom_data
	printf("-------------------------------------------\n");
	printf("There are %d atoms in the input file. \n", numatoms);
	printf("-------------------------------------------\n");
	//	printf("The data for atom 10 is %u\t %f\t %f\t %f\n",atom_data[10].index, atom_data[10].x,atom_data[10].y,atom_data[10].z);

	float dist_mat[numatoms][numatoms] ; // create n by n matrix for  distances
//	int index[numatoms]; //index of rows in distance matrix
	int nbrs_index[numatoms][5];	//create n by 5 matrix for atom and its neighbours
	int nbrs_values[numatoms][5];//create n by 5 matrix for atom and its neighbour distances

	//// calculate 2-D distance matrix        
	for(j=0; j < numatoms; j++){ 	//loop over atoms
		for(k=0; k < numatoms; k++) { //loop over atoms
			dist_mat[j][k] = getDistance(atom_data[j], atom_data[k]); //Get the paiwise distances in dist_mat		
#ifdef DEBUG
			printf("%f\t",dist_mat[j][k]);
#endif
		}
#ifdef DEBUG
		printf("\n");
#endif
	}

	struct nbrs_ind atom_nbrs[numatoms] ; //Initialize neighbours and their indices struct

	/// Find nearest neighbours from a distance matrix. 
	for (k=0; k < numatoms; k++) { // Loops over all atoms (k) 
		for (j=0; j < numatoms; j++) { // Loop over all atoms
			atom_nbrs[j].index = j;
			atom_nbrs[j].value = dist_mat[k][j];
		}         
			qsort(atom_nbrs, numatoms, sizeof(atom_nbrs[0]),cmp); //sort the rows which are distances from an item in distace matrix 	

			for(i=0; i <=4 ;i++) {
				nbrs_index[k][i] = atom_nbrs[i].index;
				nbrs_values[k][i] = atom_nbrs[i].value;					
#ifdef DEBUG
				printf("Closest four neighbours are %f\t ", nbrs_index[k][i] );
#endif
			}
#ifdef DEBUG
		printf("\n");
//		printf("Closest four neighbours are %f\t %f\t %f\t \n",nbrs_index[k] );
#endif
}
	//    printf("The atom is %d", numatoms);
	//#ifdef DEBUG
	//    for (i=0; i<=4; i++) {
	//      printf("Closest four neighbours are %f\t %f\t %f\t",nbrs_values[k][i] );	
	//    }
	//#endif
	//    printf("\n");
	//    float cosines_atom[numatoms][6];
	//    float Q_atom[numatoms];
	//    float tmp_cosine[6];
	//    int tmp_nbrs[4];
	//    int tmp_nbrs_comb[18];
	//    for(k=0; k<numatoms ; k++) {
	//      for(j=0;j<4;j++) {
	//	tmp_nbrs[j] = nbrs_index[k][j+1];
	//      }
	//      angle_comb(tmp_nbrs,tmp_nbrs_comb,k);
	//      for(i=0;i<6;i++) {
	//	cosines_atom[k][i]=getCosine(atom_data[tmp_nbrs_comb[3*i]],atom_data[tmp_nbrs_comb[3*i+1]],atom_data[tmp_nbrs_comb[3*i+2]]);
	//      }
	//      
	//      //			printf("The six cosine angles for atom %d are: \n",k);	
	//      float sum_cosine = 0;
	//      for(i=0;i<6;i++) {
	//	sum_cosine += cosines_atom[k][i];
	//      }
	//      Q_atom[k] = (3.0f/32.0f)*sum_cosine;
	//    }
	//    float sum_q = 0;
	//    FILE 	*outfile;
	//    outfile  = fopen(argv[2],"w");
	//    for (i=0; i<numatoms;i++) {
	//      fprintf(outfile,"%d\t%.4f \n ",atom_data[i].index,Q_atom[i]);
	//    }
	//    
	//    for (i=0; i<numatoms;i++) {
	//      sum_q += Q_atom[i];
	//    }
	//    fprintf(outfile, "Tetrahedral Order Parameter is %.4f. \n ",sum_q/numatoms);
	//    printf("Tetrahedral Order Parameter is %.4f. \n ",sum_q/numatoms);
	//    fclose(outfile)
	free(atom_data);
	exit(EXIT_SUCCESS);
}
 
