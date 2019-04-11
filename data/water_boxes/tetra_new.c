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
float getDistance(struct O_data a, struct O_data b)
{
    float distance;
    distance = sqrt((a.x - b.x) * (a.x - b.x) + (a.y-b.y) *(a.y-b.y) + (a.z - b.z) * (a.z - b.z));
    return distance;
}

float getCosine(struct O_data a, struct O_data c, struct O_data b)
{
      float ang, angle;
      float aa,bb,ab,b2,a2;
      aa=getDistance(c,a);
      bb=getDistance(c,b);
      ab=aa*bb;
      b2=bb*bb;
      a2=aa*aa;
      ang = (a2 + b2 - a2/(2*ab));
      angle = (ang+((1.0)/3))*(ang+((1.0)/3));	
      return angle;
}

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
  FILE *inputfile; // file pointer
  char *line = NULL; //line pointer
  size_t len = 0;
  ssize_t nread;
  struct O_data * somedata = NULL; //pointer to inout struct
  int ctr = 0; // count variable for the input data
  
  int i,j,k,p ; // loop initilizers
  
  if (argc != 2) {
    fprintf(stderr, "Usage: %s <inputfile> <outputfile>\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  inputfile = fopen(argv[1], "r");
  if (inputfile == NULL) {
    perror("fopen");
    exit(EXIT_FAILURE);
  }
  while ((nread = getline(&line, &len, inputfile)) != -1) {
    if ((somedata = realloc(somedata, (size_t) (ctr + 1) * sizeof(struct O_data))) == NULL) {
      fprintf(stderr, "error not enough memory");
      exit(EXIT_FAILURE);
    }
    deserialize_data(&somedata[ctr], line);
    ctr = ctr + 1;
  }
  // All the data is read in memory in somedata
  p = ctr -1 ;
  	printf("-------------------------------------------\n");
  	printf("There are %d atoms in the input file. \n", ctr);
  	printf("-------------------------------------------\n");
  free(line);
  fclose(inputfile);
  free(somedata);	
}
