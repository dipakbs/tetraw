// run as follows:
// gcc tetra.c -lm 
// ./a.out input_coord.dat

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//struct for input csv data
struct O_data
{
   unsigned int index;
   float x;
   float y;
   float z;
};

//Given the data in a line in a an inputfile, process it to put it in O_data struct
struct O_data * deserialize_data(struct O_data *data, const char *input)
{
   return (sscanf(input, "%u,%f,%f,%f", &data->index, &data->x, &data->y, &data->z) != 7) 
     ? NULL : data;
}
//Distance function for two pints in a struct
double getDistance(struct O_data a, struct O_data b)
{
    double distance;
    distance = sqrt((a.x - b.x) * (a.x - b.x) + (a.y-b.y) *(a.y-b.y) + (a.z - b.z) * (a.z - b.z));
    return distance;
}

double getCosine(struct O_data a, struct O_data c, struct O_data b)
{
      float ang, angle;
      ang = (getDistance(c,a)*getDistance(c,a) + getDistance(c,b)*getDistance(c,b)- getDistance(a,b)*getDistance(a,b))/(2*getDistance(c,a)*getDistance(c,b));
      angle = (ang+((1.0)/3))*(ang+((1.0)/3));	
//    SIN = sqrt(1-CS*CS);
//    CS_1_3 = CS*cos(1/3) - SIN*sin(1/3);		
//    angle = CS_1_3*CS_1_3;
    return angle;
}
//struct for neighbour distance and their indices
struct nbrs_ind {
        float value;
	int index;
};

// comparision function for sorting the neighbours
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
void angle_comb( int var1[], int var2[],int var3)
{
    int x;	
    for(x=0; x<4; x++)
    {
	var2[x] = var1[x];
    }
	var2[0] = var1[0] ; var2[1] = var3; var2[2] = var1[1];var2[3] = var1[0];var2[4] = var3;
	var2[5] = var1[2];var2[6] = var1[0];var2[7] = var3;var2[8] = var1[3];var2[9] = var1[1];
	var2[10] = var3;var2[11] = var1[2];var2[12] = var1[1];var2[13] = var3;var2[14] = var1[3] ;
        var2[15] = var1[2];var2[16] = var3;var2[17] = var1[3] ;	
}


//main program
int  main(int argc, char *argv[])
{
   FILE *stream; // file pointer
   char *line = NULL; //line pointer
   size_t len = 0;
   ssize_t nread;
   struct O_data * somedata = NULL; //pointer to inout struct
   size_t ctr = 0; // count variable for the input data

   int i,j,k,p ; // loop initilizers
		 	
   if (argc != 3) {
       fprintf(stderr, "Usage: %s <inputfile> %s <outputfile>\n", argv[0]);
       exit(EXIT_FAILURE);
   }
   stream = fopen(argv[1], "r");
   if (stream == NULL) {
       perror("fopen");
       exit(EXIT_FAILURE);
   }
	while ((nread = getline(&line, &len, stream)) != -1) {
	if ((somedata = realloc(somedata, (ctr + 1) * sizeof(struct O_data))) == NULL) {
		fprintf(stderr, "error not enough memory");
		exit(EXIT_FAILURE);
	}
	ctr = ctr + 1;
	deserialize_data(&somedata[ctr - 1], line);
 	}
// All the data is read in memory in somedata
	p = ctr -1 ;
//	printf("-------------------------------------------\n");
//	printf("There are %d atoms in the input file. \n", ctr);
//	printf("-------------------------------------------\n");
        free(line);
        fclose(stream);

        float dist_mat[ctr][ctr] ; // create n by n matrix for  distances
        int index[ctr]; //index of rows in distance matrix
        int nbrs_index[ctr][5];	//create n by 5 matrix for atom and its neighbours
// fill 2-D distnace matrix        
 	for(j=0; j < ctr; j++){
	for(k=0; k < ctr; k++) {
			dist_mat[j][k] = getDistance(somedata[j], somedata[k]);			
		}
	}
//	printf("The angle between 3 points is: %lf \n", getCosine(somedata[0],somedata[1],somedata[2]));

////	Print the atom indices from somedata 
//        for (i=0; i<ctr; i++) {
//	index[i] = i;	
//		printf("The somedata is %u \n", somedata[i]);
//	}

//	printf("-------------------------------------------\n");
//	printf("The distance matrix is below: \n");
//	printf("-------------------------------------------\n");

// 	for(j=0; j <ctr; j++){
//        	for(k=0;k<ctr;k++) {
//	        	printf("%f \t",dist_mat[j][k]);
//		}
//	                	printf("\n");
//	}

	struct nbrs_ind objects[ctr] ; //Initialize neighbours and their indices struct
	
 //       objects[ctr] = (struct nbrs_ind*) malloc (ctr * sizeof(struct nbrs_ind));
	for (k=0; k < ctr; k++) {
               for (j=0; j < ctr; j++) {
        		objects[j].index = j;
          		objects[j].value = dist_mat[k][j];
          	
        }      
  
//       printf("\n");
//       printf("The neighbour distances for atom number %d are: \n", k);
//       for (i=0; i < ctr ; i++) {
//        	printf("%f \t",objects[i].value);
//       }  
         qsort(objects, ctr, sizeof(objects[0]),cmp); //sort the rows which are distances from an item in distace matrix 
//       printf("\n");            
//       printf("THe sorted values are: \n");
//       for (i = 0; i < ctr; i++)	
//       printf("%f \t", objects[i].value);
//       printf("\n");

//	 for(i=0; i <=4 ;i++) {
//		nbrs_index[k][i] = objects[i].index;			
//	 }

//       printf("THe sorted indices are: \n");
//       for (i = 0; i < ctr; i++)
//       printf("%d \t", objects[i].index);
//	 printf("\n");

	 for(i=0; i <=4 ;i++) {
		nbrs_index[k][i] = objects[i].index;			
	 }

//	for(i=0; i < ctr ;i++)	
//        nbrs_index[k][i] = objects[i].index; 

}
//	for(i=0; i < ctr ;i++) {
//		for(j=0; j <=4 ;j++) {
//		nbrs_index[k][i] = objects[i].index;			
//		}
//	}	
//      nbrs_index[k][i] = objects[i].index; 
//	printf("-------------------------------------------\n");
//	printf("The neighbours index matrix is :\n");
//	printf("-------------------------------------------\n");
//	for (k=0; k < ctr; k++) {
//		for (i=0; i <=4 ;i++) {
//			printf("%d \t ", nbrs_index[k][i]);
//		}
//			printf("\n");
//	}
//	printf("The angle combinations are:\n");	

	float cosines_atom[ctr][6];
	float Q_atom[ctr];
	for(k=0; k<ctr ; k++) {
		int *tmp_nbrs = malloc(4*sizeof(int)) ; 	
		int *tmp_nbrs_comb = malloc(18*sizeof(int)) ;
//		printf("Tmp neighbours for atom %d are: \t ",k);
		for(j=0;j<4;j++) {
			tmp_nbrs[j] = nbrs_index[k][j+1];
//			printf("%d \t", tmp_nbrs[j]);
		}
		
		angle_comb(tmp_nbrs,tmp_nbrs_comb,k);
//			printf("%d\t  %d\t %d\t \n",tmp_nbrs_comb[0],tmp_nbrs_comb[1], tmp_nbrs_comb[2]);	
//			printf("%d\t  %d\t %d\t \n",tmp_nbrs_comb[3],tmp_nbrs_comb[4], tmp_nbrs_comb[5]);
//     			printf("%d\t  %d\t %d\t \n",tmp_nbrs_comb[6],tmp_nbrs_comb[7], tmp_nbrs_comb[8]);	
//			printf("%d\t  %d\t %d\t \n",tmp_nbrs_comb[9],tmp_nbrs_comb[10], tmp_nbrs_comb[11]);
//			printf("%d\t  %d\t %d\t \n",tmp_nbrs_comb[12],tmp_nbrs_comb[13],tmp_nbrs_comb[14]);
//			printf("%d\t  %d\t %d\t \n",tmp_nbrs_comb[15],tmp_nbrs_comb[16],tmp_nbrs_comb[17]);
			for(i=0;i<6;i++) {
			if(i==0) {
	 			cosines_atom[ctr][i] = getCosine(somedata[tmp_nbrs_comb[0]],somedata[tmp_nbrs_comb[1]],somedata[tmp_nbrs_comb[2]]);
			} else if(i==1){
				cosines_atom[ctr][i] = getCosine(somedata[tmp_nbrs_comb[3]],somedata[tmp_nbrs_comb[4]],somedata[tmp_nbrs_comb[5]]);	
			} else if(i==2){
				cosines_atom[ctr][i] = getCosine(somedata[tmp_nbrs_comb[6]],somedata[tmp_nbrs_comb[7]],somedata[tmp_nbrs_comb[8]]);	
			} else if(i==3){
				cosines_atom[ctr][i] = getCosine(somedata[tmp_nbrs_comb[9]],somedata[tmp_nbrs_comb[10]],somedata[tmp_nbrs_comb[11]]);	
			} else if(i==4){
				cosines_atom[ctr][i] = getCosine(somedata[tmp_nbrs_comb[12]],somedata[tmp_nbrs_comb[13]],somedata[tmp_nbrs_comb[14]]);	
			} else if(i==5){
				cosines_atom[ctr][i] = getCosine(somedata[tmp_nbrs_comb[15]],somedata[tmp_nbrs_comb[16]],somedata[tmp_nbrs_comb[17]]);
			} else {
			exit(EXIT_FAILURE);
			}
			
			}
//			printf("The six cosine angles for atom %d are: \n",k);	
			float *tmp_cosine = malloc(6*sizeof(float)) ; 
			float sum_cosine = 0;
//			for(i=0;i<6;i++) {
//				printf("%f\t",cosines_atom[ctr][i]);
//			}
//			printf("\n");
			for(i=0;i<6;i++) {
				sum_cosine += cosines_atom[ctr][i];
			}
//			printf("iteration %d \n",k);
//			printf("The sum of cosines is %f \n",sum_cosine);
			Q_atom[k] = (3/32.0)*sum_cosine;
//		printf("\n");
		free(tmp_nbrs);
		free(tmp_nbrs_comb);
		free(tmp_cosine);
	}
	float sum_q = 0;
	FILE 	*outfile;
	outfile  = fopen(argv[2],"w");
	for (i=0; i<ctr;i++) {
//		printf("The Q for atom %d is %f \n",somedata[i].index, Q_atom[i]);
		fprintf(outfile,"%d\t%.4f \n ",somedata[i].index,Q_atom[i]);
	}
//	fclose(outfile);

	for (i=0; i<ctr;i++) {
		sum_q += Q_atom[i];
	}
	fprintf(outfile, "Tetrahedral Order Parameter is %.4f. \n ",sum_q/ctr);
	printf("Tetrahedral Order Parameter is %.4f \n ",sum_q/ctr);
	fclose(outfile);
	free(somedata);
        exit(EXIT_SUCCESS);
}

