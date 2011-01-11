#include <stdio.h>
#include <math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_permutation.h>
#define PRINT 1
#define CNT_SQUARE 2
#define CNT_HEX_ARM 4
#define HORN_SQUARE 8
#define HORN_HEX 16
#define NUMPRINT 32
#define FINISHED 64
#define DOS 128
#define SORTED 256
#define CNT_HEX_ZIG 512
#define READFROMFILE 1024

#define INTERVAL 1 

void PrintMatrix(double * p1,int n);
void PopulateMatrix(double * p1,int n);
void Diag(double * p1, int n, int offset, int size, double value);
void ReflectMatrix(double * p1, int n);
void Calculate(double * vec1, int n, double * ret_vec,int * states);
int comparedouble(const void * a, const void * b);

int flags = 0;
double Alpha,Gamma,Epsilon,division;

int main(int argc, char * argv[])
{
	int num,i,*states,num_e,j,k;
	double *entry, *dos_vec, *eval_vec;
	double coef;
	char tmpc;
	gsl_vector *eval;
	gsl_matrix_view m;
	gsl_eigen_symm_workspace *w;
	gsl_matrix_complex * matrix2;
	gsl_complex z;
	gsl_permutation * p;
	char name[25];
	FILE *out;

	
	if ((out = fopen(argv[1],"r")) != NULL)
	{
		num = 0;
		printf("reading from file: %s....\n",argv[1]);
		while(fscanf(out,"%lf",&coef) != EOF) num++;
		if (modf(sqrt(num),&coef) != 0)
		{
			printf("not a square matrix!\n\n");
			return 1;
		}
		//printf("%d\n",num);
		fseek(out,0,SEEK_SET);
		num = sqrt(num);
		entry = (double *)malloc(sizeof(double) * num * num);
		i = 0;
		while (i < (num * num))
		{
			*(entry + i) = 0;
			i++;
		}
		i = 0;
		printf("loading matrix into memory...\n");
		while(fscanf(out,"%lf",&coef) != EOF)
		{
			*(entry + i) = coef;
			printf("entry %d = %g\n",i + 1,coef);
			i++;
		}
		flags += NUMPRINT;
		flags += READFROMFILE;
	}
	
	if (argc < 6 && (flags & READFROMFILE) == 0)
	{
		printf("usage: no. of layers, slope of cone/no. of atoms per layer, on-site energy, hopping parameter, model, (option)\n");
		printf("\n Models are:\n1: Nanotube, square lattice\n2: Nanotube, hex lattice (armchair)\n");
		printf("3: Nanohorn, square lattice\n4: Nanohorn hex lattice\n5: Nanotube, hex lattice (zig-zag)\n");
		printf("\nOptions are: p - print symbolically\n             n - print numerically\n\n");
		return 1;
	}
	else if ((flags & READFROMFILE) != 0)
	{
		
	}
	else
	{
		num = atoi(argv[1]); 
		Epsilon = atol(argv[3]);
		Gamma = atol(argv[4]);
		switch(atoi(argv[5]))
		{
			case 1: flags = flags + CNT_SQUARE;
				Alpha = atoi(argv[2]);
				num *= Alpha;
				entry = (double *)malloc(sizeof(double) * num * num);
				break;
			case 2: flags = flags + CNT_HEX_ARM;
				Alpha = atoi(argv[2]);
				i = Alpha;
				if ((i % 2) != 0)
				{
					printf("no. of atoms per layer must be even, changing to %d\n",i + 1);
					Alpha += 1;
				}
				num *= (Alpha);
				entry = (double *)malloc(sizeof(double) * num * num);
				break;
			case 3: flags = flags + HORN_SQUARE;
				i = 1;
				num_e = 0;
				while (i <= num)
				{
					num_e += i;
					i++;
			}
				num = num_e;
				entry = (double *)malloc(sizeof(double) * num * num);
				break;
			case 4: flags = flags + HORN_HEX;
				num = 0;
				entry = NULL;
				break;
			case 5: flags = flags + CNT_HEX_ZIG;
				Alpha = atoi(argv[2]);
				i = Alpha;
				while ((i % 4) != 0)
				{
					i++;
				}
				Alpha = i;
				num *= Alpha;
				entry = (double *)malloc(sizeof(double) * num * num);	
				break;
			default: printf("No valid model choice!\n");
				return 1;
				break;
				  
		}
		i = 0;
		while (i < (num * num))
		{
			*(entry + i) = 0;
			i++;
		}
	printf("Alpha = %g\n",Alpha);
	}
	printf("%d atoms\n",num);
	if (argc > 6)
	{
		if (strcmp("p\0",argv[6]) == 0) 
		{
			//printf("print");
			flags = flags + PRINT;
		}
		if (strcmp("n\0", argv[6]) == 0)
		{
			flags = flags + NUMPRINT;
		}	
	}
	eval = gsl_vector_alloc(num);


	//printf("flags %d\n",flags); 
	//k = (flags & PRINT);
	//printf("print ? %d\n",k); 
	//printf("%lf",Alpha);
	//printf("%d",num);
	
	if (entry == NULL)
	{
		printf("could not allocate memory\n");
		return 1;

	}
	i = 0;
	if ((flags & READFROMFILE) == 0) 
	{
		printf("Populating Matrix...\n");
		PopulateMatrix(entry,num);
		out = fopen("matrix","w");
		i = 0;
		while (i < (num * num))
		{
			fprintf(out,"%3.5lf ",*(entry + i));
			if (((i + 1) % num) == 0) fprintf(out,"\n");
			i++;
		}
		fclose(out);
	}
	else
	{
		
	}
	printf("Printing Matrix\n");
	
	PrintMatrix(entry,num);

	//CalcEigenvalues(entry,num,eval);
	
	
	
	//Eigenvalue Calculation
	m = gsl_matrix_view_array (entry,num,num);
	//printf("87\n");
        w = gsl_eigen_symm_alloc (num);
	//printf("89\n");
	printf("Solving Eigenvalue equation....\n");
        if(gsl_eigen_symm (&m.matrix, eval, w)) 
	{
		printf("an error occurred solving the eigenvalue equation\n");
		return 1;
	}	
	//printf("\nworkspace location = %#x",w);
	gsl_eigen_symm_free(w);


	
	printf("\n\n");
	//if (flags & NUMPRINT != 0 || flags & PRINT != 0)
	printf("Eigenvalues found. ");
	eval_vec = malloc(sizeof(double) * num);
	i = 0;
	while (i < num)
	{
		*(eval_vec + i) = gsl_vector_get(eval,i);
		i++;
	}
	gsl_vector_free(eval);
	tmpc = 0;
	printf("Ordering Eigenvalues...\n");
	qsort(eval_vec,num,sizeof(double),comparedouble);
	division = ((*(eval_vec + num - 1) - *eval_vec) / num) * 10 * INTERVAL;
	num_e = ((*(eval_vec + num - 1) - *eval_vec))/ division;
	num_e++; //just in case
	while ((flags & FINISHED) == 0)
	{
		if (tmpc != -1) printf("(s)ave in file, (q)uit and discard, (p)rint (!), (d)ensity of states: ");
		scanf("%c",&tmpc);
		if (tmpc == 'p' || tmpc == 'P')
		{
			i = 0;
			while (i < num)
			{
				printf("Eigenvalue %d\t: %g\n",(i + 1),*(eval_vec + i));
				i++;
			}
			printf("\n");
			tmpc = -1;
		}
		else if (tmpc == 'q' || tmpc == 'Q')
		{
			printf("\nexiting...\n\n");
			flags += FINISHED;
		}
		else if (tmpc == 'd' || tmpc == 'D')
		{
			tmpc = 0;
			dos_vec = malloc(sizeof(double) * num_e);
			states = malloc(sizeof(int) * num_e);
			printf("\nCalculating density of states...\n");
			Calculate(eval_vec,num,dos_vec,states);
			i = 0;
			printf("\n");
			out = fopen("dostube","w");
			while (i < num_e && *(states + i) >= 0)
			{
				printf("Energy: % .5lf\tNo. of states: %d\n",*(dos_vec + i),*(states + i));
				fprintf(out,"% .5lf\t%d\n",*(dos_vec + i),*(states + i));
				i++;
			}
			num_e = i;
			fclose(out);

			flags += FINISHED;
			flags += DOS;
		}
		else if (tmpc == -1) //WHY?!!
		{
			tmpc = 0;
		}
		else if (tmpc == 's' || tmpc == 'S') //broken for some reason
		{
			printf("\nchoose filename: ");

			scanf("%s",name);
			//printf("%s",name);
			out = fopen(name,"r");
			if (out != NULL)
			{
				printf("File already exists!\n");
				fclose(out);
			}
			else 
			{
				fclose(out);
				out = fopen(name,"w");
				printf("Writing to file: %s",name);
				i = 0;
				while (i < num)
				{
					fprintf(out,"Eigenvalue %d\t: %g\n",(i + 1),*(eval_vec + i));
					i++;
				}
				fclose(out);
				tmpc = -1;
			}
		}
		else
		{
			printf("unrecognised option!\n");
		}
	}
	if ((flags & DOS) != 0)
	{
		while (j < num_e)
		{
			matrix2 = gsl_matrix_complex_alloc(sizeof(gsl_complex) * num,sizeof(gsl_complex) * num);
			i = 0;
			gsl_matrix_complex_set_identity(matrix2);
			z.dat[0] = *(dos_vec + j);
			z.dat[1] = 0.0001;
			gsl_matrix_complex_scale(matrix2,z);
			while (i < num)
			{
				k = 0;
				while(k < num)
				{
					z = gsl_matrix_complex_get(matrix2,i,k);
					z.dat[0] -= *(entry + i + (k * num));
					gsl_matrix_complex_set(matrix2,i,k,z);
					k++;		
				}
				i++;
			}
			out = fopen("matrix","w");
			gsl_matrix_complex_fprintf(out,matrix2,"% 2.2lf");
			fclose(out);
			free(entry);
		}
		
	}

	//printf("\neval location = %#x",eval);
	return 0;
}

void PrintMatrix(double * p1,int n)
{
	int i = 0, j = 0;
	if ((flags & NUMPRINT) != 0)
	{	
		while (i < n)
		{
			j = 0;
			while (j < n)
			{
				printf("% 2.2lf ",*(p1 + j + (i * n)));
				j++;
			}
			printf("\n");
			i++;
		}
		return;
	}
	if ((flags & PRINT) != 0)
	{
		
		while (i < n)
		{
			j = 0;
			while (j < n)
			{
				if(*(p1 + j + (i * n)) == Gamma && Gamma != 0) 
				{
					printf("%c ",116);
				}
				else if(*(p1 + j + (i * n)) == Epsilon && Epsilon != 0)
				{
					printf("%c ",101);
				}
				else if(*(p1 + j + (i * n)) == 0) 
				{
					printf("0 ");
				}				
				j++;
			}
			printf("\n");
			i++;
		}
		return;
	}
	
	//printf("did not print");
	return;	

}
void PopulateMatrix(double *p1,int n)
{
	int i = 0,j,k;
	//double * block = NULL; do i need this? i don't know
	if ((flags & CNT_SQUARE) != 0)
	{
		j = Alpha;
		while (i < n)
		{
			Diag(p1,n,i + (i * n),j,Epsilon);	//Diagonals = epsilon (on-site energy)
			*(p1 + i + j - 1 + (i * n)) = Gamma;
			Diag(p1,n,i + 1 + (i * n),j - 1,Gamma);
			Diag(p1,n,i + j + (i * n),j,Gamma);
			i += j;
		}
		ReflectMatrix(p1,n);
	}
	else if ((flags & CNT_HEX_ARM) != 0)
	{
		j = Alpha;
		while ((i * j) < n)
		{
			Diag(p1,n,(i * j) + (i * j * n),j,Epsilon);
			if (i % 2 == 0)
			{
				*(p1 + ((2 + i) * j) - 1 + (i * j * n)) = Gamma;
				Diag(p1,n,((i + 1) * j) + (((i * j) + 1) * n),j - 1,Gamma);
				//printf("even\n");
			}	
			else
			{
				*(p1 + ((i + 1) * j) + (((i + 1) * j) - 1) * n) = Gamma;
				Diag(p1,n,((i + 1) * j) + 1 + (i * j * n),j - 1,Gamma); 				
				//printf("odd\n");
			}
			i += 1;
		}
		i = 0;
		while (i < n)
		{
			*(p1 + i + 1 + (i * n)) = Gamma;
			i += 2;
		}
		ReflectMatrix(p1,n);
	}
	else if ((flags & CNT_HEX_ZIG) != 0)
	{
		j = Alpha;
		i = 0;
		while (i < n)
		{
			Diag(p1,n,i + (i * n),j,Epsilon);
			Diag(p1,n,i + 1 + (i * n),j - 1,Gamma);
			*(p1 + i + j - 1 + (i * n)) = Gamma;
			if ((i / j) % 2 == 0)
			{
				//printf("even");
				k = 0;
				while (k < j)
				{
					*(p1 + i + k + j +((i + k) * n)) = Gamma;
					k += 2;
				}
			}	
			else
			{
				//printf("odd");
				k = 0;
				while (k < j)
				{
					*(p1 + i + 1 + k + j + ((i + k + 1) * n)) = Gamma;
					k += 2;
				}
			}
			i += j;
		}
		ReflectMatrix(p1,n);
	}
	else if ((flags & HORN_SQUARE) != 0)
	{
		j = 1;
		i = 0;
		while(i < n)
		{
			//printf("%d",i);
			Diag(p1,n,i + (i * n),j,Epsilon);
			if (i >= 3) Diag(p1,n,i + 1 + (i * n),(j - 1),Gamma);
			if (i > 0) *(p1 + i + j - 1 + (i * n)) = Gamma;
			k = 0;
			if ((i + j) < n) 
			{
				while (k < j)
				{
					*(p1 + i + k + j + ((i + k) * n)) = Gamma;
					*(p1 + i + k + j + 1 + ((i + k) * n)) = Gamma;
					k++;
				}
			}	
			i += j;
			j++;
		}
		ReflectMatrix(p1,n);
	}

	return;
}	

void Diag(double * p1, int n, int offset, int size, double value)
{
	int i = 0;
	while (i < size)
	{
		*(p1 + i + (i * n) + offset) = value;
		i++;
	}
	return;
}


void ReflectMatrix(double * p1, int n)
{
	int i = 1,j;
	while (i < n)
	{
		j = 0;
		while (j < n)
		{
			if(i - j > 0) *(p1 + j + (i * n)) = *(p1 + i + (j * n));
			j++;
		}
		i++;
	}
}
/* void Calculate(double * vec1,int n, double * ret_vec) //Calculation of DOS, not sure about this yet
{
	int i = 0;
	double tmpd;
	tmpd = *vec1 - *(vec1 + 1);
	if (tmpd == 0) tmpd = 1/MAXDOS;
	while (i + 1 < n)
	{
		tmpd = (*(vec1 + i) - *(vec1 + i + 1));
		if (tmpd == 0) tmpd = 1 / MAXDOS;
		*(ret_vec + i) = 1 / tmpd;
		*(ret_vec + i + n) = *(vec1 + i);
		i++;
	}
	tmpd = *(vec1 + n - 2) - *(vec1 + n - 1);
	if (tmpd == 0) tmpd = 1/MAXDOS;
	*(ret_vec + n - 1) = 1 / tmpd;
	return;
}
*/

void Calculate(double * vec1, int n, double * ret_vec,int * states)
{
	double E = *vec1;
	int i,j = 0;
	while (E < (*(vec1 + n - 1))) 
	{
		*(states + j) = 0;
		i = 0;
		while (i < n)
		{
			if (*(vec1 + i) < (E + division) && *(vec1 + i) > E) *(states + j) += 1;
			i++;
		}
		*(ret_vec + j) = E;
		j++;
		E += division;
	}
	return;
}

int comparedouble(const void * a, const void * b)
{
	double * arga = (double *) a;
	double * argb = (double *) b;
	if (*arga < *argb) return -1;
	else if (*arga == *argb) return 0;
	else return 1;
}
