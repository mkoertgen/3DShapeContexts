#ifndef GLOBALMATCHING_H
#define GLOBALMATCHING_H

#include <iostream>
//#include <stdio.h>
//#include <stdlib.h>
//#include <fstream>

#include "MathClasses.h"
#include "SampleArr.h"
#include "ShapeMatrix.h"

/**
* \class GlobalMatching
* \brief template GlobalMatching class
*
* This class is a template GLOBALMATCHING class used the Hungarin Method
*/
template<class Type> class GlobalMatching
{
	public:
		template<class Type> class Inhalt
		{
			public:
				Type value;  //  values of the matrix elements
				bool assign;
				bool zeroassignment;
				bool mark;
				bool linie;
				bool markzeile;
				bool markspalte;	
		};

		Inhalt<Type> **matrix;    /// matrix elements
		int rows;    /// number of rows for the matrix
		int columns; /// number of rows for the matrix
        int dim;

		int* liste;  // list for permutation
		int* minassign; // minimum assignment list
		Type Smin;    // solution for Global Matching with permutation

		typedef Type cost;		
		
		typedef int row;
		typedef int col;

		cost **assigncost, *u, *v, lapcost;
		row i, *colsol;
		col j, *rowsol;


		int* k;  //  list for multi assignment
        int**  minassign2;  //minimum assignment list for Multi Assignment
		Type** w;           // weights for costs     
                       


		
	/// constructor, sets all entries 0
	GlobalMatching() 
	{
		rows=columns=dim=0;
		liste=minassign=NULL;
	}

	/// destructor
	~GlobalMatching()
	{
		if (matrix)
		{
			deleteEntries();
			delete [] matrix;
		}
	}


//**************************************************************************************************************
//------------- GlobalMatching with shortest augmenting path algorithm for linear assignment problem -----------
//**************************************************************************************************************

	Type LinearAssignmentMethod(const ShapeMatrix<Type>* A)
	{	
        if (A->rows > A->columns) 
			dim = A->rows;
		else
			dim = A->columns;

		assigncost = new cost*[dim];
		for (i = 0; i < dim; i++)
			assigncost[i] = new cost[dim];

		rowsol = new col[dim];
		colsol = new row[dim];
		u = new cost[dim];
		v = new cost[dim];

		for (i = 0; i < A->rows; i++)
			for (j = 0; j < A->columns; j++)
					assigncost[i][j] = A->m[i][j]+1;  // add 1 to prevent too little assigncost

        if (A->rows > A->columns)			
			for (i = 0; i < A->rows; i++)
				for (j = A->columns; j < A->rows; j++)
					assigncost[i][j] = (Type)1;       // not zero to prevent too little assigncost

		if (A->columns > A->rows)
			for (i = A->rows; i< A->columns; i++)
				for (j = 0; j < A->columns; j++)
					assigncost[i][j] = (Type)1;

/*		for (i=0; i<dim; i++)
			if (assigncost[i])
			{
				printf("[ %.2f",assigncost[i][0]);
				for (j=1; j<dim; j++)
					printf(", %.2f", assigncost[i][j]);
				printf(" ]\n");
			} else printf("NULL\n");*/

		lapcost = lap(dim, assigncost, rowsol, colsol, u, v);

		delete[] assigncost;
		delete[] rowsol;
		delete[] colsol;
		delete[] u;
		delete[] v;

		return lapcost;
	}

	Type lap(int dim, cost **assigncost,col *rowsol, row *colsol, cost *u, cost *v)
	{
		bool unassignedfound;
		row  i, imin, numfree = 0, prvnumfree, f, i0, k, freerow, *pred, *free;
		col  j, j1, j2, endofpath, last, low, up, *collist, *matches;
		cost min, h, umin, usubmin, v2, *d, BIG;

		free = new row[dim];       // list of unassigned rows.
		collist = new col[dim];    // list of columns to be scanned in various ways.
		matches = new col[dim];    // counts how many times a row could be assigned.
		d = new cost[dim];         // 'cost-distance' in augmenting path calculation.
		pred = new row[dim];       // row-predecessor of column in augmenting/alternating path.

		// init how many times a row will be assigned in the column reduction.
		for (i = 0; i < dim; i++)  
			matches[i] = 0;

		// COLUMN REDUCTION 
		for (j = dim-1; j >= 0; j--)    // reverse order gives better results.
		{
			// find minimum cost over rows.
			min = assigncost[0][j]; 
			imin = 0;
			for (i = 1; i < dim; i++)  
				if (assigncost[i][j] < min) 
				{ 
					min = assigncost[i][j]; 
					imin = i;
				}
			v[j] = min; 

			if (++matches[imin] == 1) 
			{ 
				// init assignment if minimum row assigned for first time.
				rowsol[imin] = j; 
				colsol[j] = imin; 
			}
			else
				colsol[j] = -1;        // row already assigned, column not assigned.
		}

		BIG = 0;
		for (i=0; i < dim; i++)
			for (j=0; j < dim; j++)
				BIG += assigncost[i][j];  // to find a upper limit 

		// REDUCTION TRANSFER
		for (i = 0; i < dim; i++) 
			if (matches[i] == 0)     // fill list of unassigned 'free' rows.
				free[numfree++] = i;
			else
			if (matches[i] == 1)   // transfer reduction from rows that are assigned once.
			{
				j1 = rowsol[i]; 
				min = BIG;
				for (j = 0; j < dim; j++)  
					if (j != j1)
						if (assigncost[i][j] - v[j] < min) 
						min = assigncost[i][j] - v[j];
				v[j1] = v[j1] - min;
		}

		// AUGMENTING ROW REDUCTION 
		int loopcnt = 0;           // do-loop to be done twice.
		do
		{
			loopcnt++;

		    // scan all free rows.
			// in some cases, a free row may be replaced with another one to be scanned next.
			k = 0; 
			prvnumfree = numfree; 
			numfree = 0;             // start list of rows still free after augmenting row reduction.
			while (k < prvnumfree)
			{
				i = free[k]; 
				k++;

				// find minimum and second minimum reduced cost over columns.
				umin = assigncost[i][0] - v[0]; 
				j1 = 0; 
				usubmin = BIG;
				for (j = 1; j < dim; j++) 
				{
					 h = assigncost[i][j] - v[j];
					 if (h < usubmin)
						if (h >= umin) 
						{ 
							usubmin = h; 
							j2 = j;
						}
					else 
					{ 
						usubmin = umin; 
						umin = h; 
						j2 = j1; 
						j1 = j;
					}
				}

				i0 = colsol[j1];
				if (umin < usubmin) 
				// change the reduction of the minimum column to increase the minimum
				// reduced cost in the row to the subminimum.
						v[j1] = v[j1] - (usubmin - umin);
				else                   // minimum and subminimum equal.
					if (i0 >= 0)         // minimum column j1 is assigned.
					{ 
					// swap columns j1 and j2, as j2 may be unassigned.
						j1 = j2; 
						i0 = colsol[j2];
					}

					 // (re-)assign i to j1, possibly de-assigning an i0.
			        rowsol[i] = j1; 
                    colsol[j1] = i;

                    if (i0 >= 0)           // minimum column j1 assigned earlier.
						if (umin < usubmin) 
						// put in current k, and go back to that k.
						// continue augmenting path i - j1 with i0.
							free[--k] = i0; 
						else 
			           // no further augmenting reduction possible.
					   // store i0 in list of free rows for next phase.
							free[numfree++] = i0; 
			}
		}
		while (loopcnt < 2);       // repeat once.

		// AUGMENT SOLUTION for each free row.
		for (f = 0; f < numfree; f++) 
		{
			freerow = free[f];       // start row of augmenting path.

	       // Dijkstra shortest path algorithm.
		   // runs until unassigned column added to shortest path tree.
		   for (j = 0; j < dim; j++)  
		   { 
				d[j] = assigncost[freerow][j] - v[j]; 
				pred[j] = freerow;
				collist[j] = j;        // init column list.
		   }

		   low=up=0; // columns in 0..low-1 are ready, now none.
		             // columns in low..up-1 are to be scanned for current minimum, now none.
				     // columns in up..dim-1 are to be considered later to find new minimum, 
					 // at this stage the list simply contains all columns 
		   unassignedfound = false;
		   do
		   {
				if (up == low)         // no more columns to be scanned for current minimum.
				{
					last = low - 1; 

					// scan columns for up..dim-1 to find all indices for which new minimum occurs.
					// store these indices between low..up-1 (increasing up). 
					min = d[collist[up++]]; 
					for (k = up; k < dim; k++) 
					{
						j = collist[k]; 
						h = d[j];
						if (h <= min)
						{
							if (h < min)     // new minimum.
							{ 
								up = low;      // restart list at index low.
								min = h;
							}
							// new index with same minimum, put on undex up, and extend list.
							collist[k] = collist[up]; 
							collist[up++] = j; 
						}
					}

					// check if any of the minimum columns happens to be unassigned.
					// if so, we have an augmenting path right away.
					for (k = low; k < up; k++) 
					if (colsol[collist[k]] < 0) 
					{
						endofpath = collist[k];
						unassignedfound = true;
						break;
					}
				}

				if (!unassignedfound) 
				{
					// update 'distances' between freerow and all unscanned columns, via next scanned column.
					j1 = collist[low]; 
					low++; 
					i = colsol[j1]; 
					h = assigncost[i][j1] - v[j1] - min;

					for (k = up; k < dim; k++) 
					{
						j = collist[k]; 
						v2 = assigncost[i][j] - v[j] - h;
						if (v2 < d[j])
						{
							pred[j] = i;
							if (v2 == min)   // new column found at same minimum value
								if (colsol[j] < 0) 
								{
									// if unassigned, shortest augmenting path is complete.
									endofpath = j;
									unassignedfound = true;
									break;
								}
								// else add to list to be scanned right away.
								else 
								{ 
									collist[k] = collist[up]; 
									collist[up++] = j; 
								}
							d[j] = v2;
						}
					}
				} 
			}
			while (!unassignedfound);

			// update column prices.
			for (k = 0; k <= last; k++)  
			{ 
				j1 = collist[k]; 
				v[j1] = v[j1] + d[j1] - min;
			}

			// reset row and column assignments along the alternating path.
			do
			{
				i = pred[endofpath]; 
				colsol[endofpath] = i; 
				j1 = endofpath; 
				endofpath = rowsol[i]; 
				rowsol[i] = j1;
			}
			while (i != freerow);
		}

		// calculate optimal cost.
		cost lapcost = 0;
		cost max = 0.0;

		minassign = new int[dim];

		for (i = 0; i < dim; i++)  
		{
			j = rowsol[i];
//            printf("%1d%1d ",i,j) ;
			minassign[i] = j;
			u[i] = assigncost[i][j] - v[j];
			lapcost += assigncost[i][j]-1; 
		}

		// free reserved memory.
		delete[] pred;
		delete[] free;
		delete[] collist;
		delete[] matches;
		delete[] d;

		lapcost /= (Type)dim;
		return lapcost;
	}

//**********************************************************************************************
// ----------------- GlobalMatching with Hungarian Method --------------------------------------
//**********************************************************************************************

	Type HungarianMethod(const ShapeMatrix<Type>* A)
	{
		rows = A->rows;   // initialize 
		columns = A->columns;

        if (A->rows > A->columns) 
			dim = A->rows;
		else
			dim = A->columns;

		matrix = new Inhalt<Type>*[dim];

		for (i=0; i<dim; i++)
			matrix[i] = new Inhalt<Type>[dim];
    		
		deleteall(matrix);

		for (i=0; i<A->rows; i++) 
			for (j=0; j<A->columns; j++) 
				matrix[i][j].value=A->m[i][j]+1;

        if (A->rows > A->columns)			
			for (i = 0; i < A->rows; i++)
				for (j = A->columns; j < A->rows; j++)
					matrix[i][j].value = (Type)1;       // not zero to prevent too little assigncost

		if (A->columns > A->rows)
			for (i = A->rows; i< A->columns; i++)
				for (j = 0; j < A->columns; j++)
					matrix[i][j].value = (Type)1;

		reducedmatrix(matrix);  // reduce matrix

		Type result = (Type)0;

		minassign = new int[dim];

		for (i=0; i< dim; i++)
			for (j=0;j<dim; j++)
				if (matrix[i][j].zeroassignment == true)
				{
					minassign[i] = j;
					result+= A->m[i][j];
				}

		result /= (Type) rows;
		return result;
	}

private:
    void deleteEntries()
	{
		for (i=0; i<rows; i++)
			if (matrix[i])
			{
				delete [] matrix[i];
				matrix[i] = NULL;
			}
	}

	void allocEntries(int rows, int columns)
	{
		if (rows>0 && columns>0)
		{
			this->rows = rows;
			this->columns = columns;
			matrix = new Inhalt<Type> *[rows];
			for (i=0; i < rows; i++) 
				matrix[i] = new Inhalt<Type>[columns];
			this->null();
		}
	}

// -------------- to find the minimum in a row ----------------------------------	
	Type findminrow (Inhalt<Type>** matrix, int row)
	{
		Type min = matrix[row][0].value;
		for (j = 1; j< dim; j++)
		if (min > matrix[row][j].value) 
			min = matrix[row][j].value;
//		cout << "min: " << min << " ";
		return min;
	}

//-------------- to find the minimum in a column --------------------------------
	Type findmincolumn(Inhalt<Type>** matrix, int column)
	{
		Type min = matrix[0][column].value;
		for (i=1; i< dim; i++)
		 if (min > matrix[i][column].value) min = matrix[i][column].value;
		return min;
	}

//-------------- to delete all assigns and marks ------------------------------------------
	void deleteall(Inhalt<Type>** matrix)
	{
		for (i=0;i<dim;i++)
			for (j=0;j<dim;j++)
			{
				matrix[i][j].mark=false;
				matrix[i][j].markspalte=false;
				matrix[i][j].markzeile=false;
				matrix[i][j].assign = false;
				matrix[i][j].zeroassignment = false;
			}
	}

//-------------- to assign the row and the column --------------------------
	void assign(Inhalt<Type>** matrix,int row, int column)
	{
		for (i=0;i<dim;i++)
		{
			matrix[i][column].assign = true;
		}
		for (j=0;j<dim;j++)
		{
			matrix[row][j].assign = true;
		}
	}

//-------------- to find a zero,which stands alone in a row -------------------
	bool SearchZeroAloneZeile(Inhalt<Type>** matrix) 
	{
		int i=0;
		int j;
		int S;
		int Zero=0;  // number of zeroes

		bool found=false;
		while (i<dim && found == false)
		{
			for (j=0;j<dim;j++)
			{
				if (matrix[i][j].value == 0 && matrix[i][j].assign != true)
				{
					Zero++;
					S=j;
				}
			}
			if (Zero == 1)
			{
				found = true;
				matrix[i][S].zeroassignment = true;

				assign(matrix,i,S);
			}
			else 
				Zero = 0;
			i++;
		}
		return found;	
	}

// ------------- to find a zero, which stands alone in a column -----------------
	bool SearchZeroAloneSpalte(Inhalt<Type>** matrix)
	{
		int j=0;
		int i=0;
		int k=0;
		int AnzahlZero=0;
		bool found = false;

		while (j<dim && found == false)
		{
			for (i=0; i<dim; i++)
			{
				if (matrix[i][j].value == 0 && matrix[i][j].assign != true)
				{
					AnzahlZero++;
					k=i;
				}
			}

			if (AnzahlZero == 1)
			{
				found = true;
				matrix[k][j].zeroassignment = true;
				assign(matrix,k,j);
			}
			else 
				AnzahlZero = 0;
		
			j++;
		}
		return found;
	}

// --------------------- to find a unassigned zero -----------------------------
	bool SearchAnyZero(Inhalt<Type>** matrix)
	{
		int number=0;

		for (i=0;i<dim;i++)
			for (j=0;j<dim;j++)
				if (matrix[i][j].value == 0 &&matrix[i][j].assign != true)
					number++;

		int z = rand()%(number-1);
//		cout << " Zufallszahl: " << z << endl;
        
		number =0;

		for (i=0;i<dim;i++)
			for (j=0;j<dim;j++)
				if (matrix[i][j].value == 0 &&matrix[i][j].assign != true)
				{
					if (number != z)
						number++;
					else
					{
						matrix[i][j].zeroassignment = true;
						assign(matrix,i,j);
						return true;   // zero found
					}
				}
		return false;  // no zero found
	}

// --------------------  to find Zeros in a reduced matrix ------------------------
	int SearchZeroes(Inhalt<Type>** matrix)
	{
		int i=0,j=0;
		int assignment =0;
		bool found = true;
	
		while (found)
		{
			while (SearchZeroAloneZeile(matrix))
				assignment++;

			while (SearchZeroAloneSpalte(matrix))
				assignment++;

			if( SearchAnyZero(matrix))
			{
				assignment ++;
			}
			else
				found = false;  // no zero anymore
		}
		return assignment;
	}	
	
// ------------------- to mark a row -------------------------------------
	void markiereZeile(Inhalt<Type>** matrix,int i)
	{
		for (j=0;j<dim;j++)
		{
			matrix[i][j].markzeile = true;
		}
	}

// ------------------- to mark a column ----------------------------------
	void markiereSpalte(Inhalt<Type>** matrix,int j)
	{

		for (i=0;i<dim;i++)
		{
			matrix[i][j].markspalte = true;
		}

		for (i=0;i<dim;i++)
		{
			if (matrix[i][j].zeroassignment == true)
			{
				if (matrix[i][j].markzeile != true)
				{
					markiereZeile(matrix,i);
					findSpalte(matrix);
				}
			}
		}	
	}	

// -------------------- to find zeros in marked columns --------------------------  
	void findSpalte(Inhalt<Type>** matrix)
	{
		for (i=0;i<dim;i++)
		{
			if (matrix[i][0].markzeile == true)
			{
				for(j=0;j<dim;j++)
				{
					if(matrix[i][j].value == 0)
					{
						if (matrix[i][j].markspalte != true)
						{
							markiereSpalte(matrix,j);
						}
					}
				}
			}
		}
	}	

// ------------------ to find unassigned rows and mark them -----------------------------
	void findrows(Inhalt<Type>** matrix)
	{
		for (i=0; i<dim; i++)
		{
			for (j=0; j<dim; j++)
				if(!matrix[i][j].assign)
				{
					markiereZeile(matrix,i);
				}
		}
	}

// -------------------- to find the minimum of the uncovered elements -----------------------
	Type findmin(Inhalt<Type>** matrix)
	{
		Type min=0;
		for (i=0; i<dim; i++)
		{
			for (j=0;j<dim;j++)
			{
				if (matrix[i][j].markspalte == false && matrix[i][j].markzeile == true)
				{
					if (min == 0)
						min = matrix[i][j].value;
					else 
						if (min > matrix[i][j].value)
							min = matrix[i][j].value;
				}
			
			}
		}

		return min;
	}


// -------------- add the minimum to all elements that are double covered (horizontal and vertical covered)
// -------------- subtract from all elements that are not covered  -------------------------------------
	void subadd(Inhalt<Type>** matrix, Type min)
	{
		for (i=0; i<dim; i++)
		{
			if (matrix[i][0].markzeile == true)
				for (j=0;j<dim;j++)
					matrix[i][j].value=matrix[i][j].value-min;
		}	
		for (j=0;j<dim; j++)
		{
			if (matrix[0][j].markspalte == true)
				for (i=0;i<dim;i++)
					matrix[i][j].value=matrix[i][j].value+min;
		}
	}


// ---------------- reduce the matrix by finding minimums in a row or column and substract them from all matrix elements
// ---------------- next find a cover (means find the minimum of rows and columns that cover the zeroes)
	void reducedmatrix(Inhalt<Type>** matrix)
	{
		int i,j;
		Type min;
		int assignment=0;
    
		deleteall(matrix);

		for (j=0; j< dim; j++)
		{
			min = findminrow(matrix,j);
			for (i=0;i< dim; i++) 
				matrix[j][i].value= matrix[j][i].value - min;
		}									

	    for (i=0; i< dim; i++) 
		{
		     min = findmincolumn(matrix,i);
			 for (j=0; j<dim; j++) 
				 matrix[j][i].value = matrix[j][i].value - min;
		}

		assignment=SearchZeroes(matrix);  // if number of assignments is equal to the number of dim stop
	  
		if (assignment < dim)  // not enough assignments, so subadd the matrix (look above) and go on
		{
			findrows(matrix);	  
			findSpalte(matrix);	          
			min=findmin(matrix);
			subadd(matrix,min);
			reducedmatrix(matrix);
		}
	}	



//****************************************************************************************
//  ------------------- GlobalMatching with Permutation ----------------------------------
//****************************************************************************************

	
	public:

	Type EasyMethod(const ShapeMatrix<Type>* A)
	{
		rows = A->rows;   // initialize 
		columns = A->columns;

        if (A->rows > A->columns) 
			dim = A->rows;
		else
			dim = A->columns;



		liste = new int[dim];
		for (i=0; i<dim; i++)
			liste[i]=i;


		assigncost = new cost*[dim];

		for (i=0; i<dim; i++)
			assigncost[i] = new cost[dim];


		
		for (i = 0; i < A->rows; i++)
			for (j = 0; j < A->columns; j++)
					assigncost[i][j] = A->m[i][j];  

        if (A->rows > A->columns)			
			for (i = 0; i < A->rows; i++)
				for (j = A->columns; j < A->rows; j++)
					assigncost[i][j] = (Type)0;       

		if (A->columns > A->rows)
			for (i = A->rows; i< A->columns; i++)
				for (j = 0; j < A->columns; j++)
					assigncost[i][j] = (Type)0;

        Smin = (Type)0;

		for (i=0; i<dim; i++)
			Smin += assigncost[i][i];  // to set a minimum value
        
		perm(liste,0,dim-1);


		delete [] liste;

		return Smin;
	}


	void minimum (int* was, int n)
	{ 
		int s ;
		Type min=0;
		minassign = new int[dim];
		
		for ( s = 0 ; s <= n ; s ++ )
		{
			min += assigncost[s][was[s]];
		}
		min= min/(Type)(n+1);



		if (min < Smin)
		{
			Smin = min;

			for (i = 0; i <= n; i++)  
			{
				minassign[i] = was[i];

			}

		}
		return;
	}


	void perm (int* wer, int k, int n)
	{ 
		int a, b ;

        b = wer [k] ;
        for ( a = k ; a <= n ; a ++ )
        { 
			wer [k] = wer [a] ; 
			wer [a] = b ;
            if ( k != n ) 
				perm (wer, k + 1, n) ;
			else 
				minimum(wer, n) ;
            
			wer [a] = wer [k] ;
         }
		wer [k] = b ;
		return ;
   }



//****************************************************************************************
//  ------------------- GlobalMatching with multiple assignments ----------------------------------
//****************************************************************************************











	Type MultiAssignmentMethod(const ShapeMatrix<Type>* A)
	{
		Type* min;
		Type* max;
		Type* eps;

		int i,j;
		rows = A->rows;   // initialize 
		columns = A->columns;

		assigncost = new cost*[rows];
        minassign2 = new int*[rows];

		min = new cost[A->rows];
		max = new cost[A->rows];
        eps = new cost[A->rows];
		w   = new cost*[A->rows];

		k = new int[A->rows];

  		for (i = 0; i < A->rows; i++)
		{
			// erwartungswert und streuung
     		Type mid,sigma = 0;

			// find row min/max and compute variance in row[i]
			max[i]=min[i] = A->m[i][0];

			for (j = 0; j < A->columns; j++)
			{
				mid += A->m[i][j];
				if (A->m[i][j] < min[i]) min[i] = A->m[i][j];
				if (A->m[i][j] > max[i]) max[i] = A->m[i][j];
			}
			mid /= (Type)A-> columns; // normalize

			// compute sigma
			for (j=0; j<A->columns; j++)
				sigma += (SQR(min[i]-A->m[i][j]));
			sigma = sqrt(sigma);
//			sigma /= (Type)A->columns;

//			printf("\n sigma(%i) = %f",i, sigma);

			// compute threshold for row i
//			eps[i] = A->rows/(Type)A->columns * SQR(max[i]-min[i]);
			eps[i] = sigma * (max[i]-min[i]);
		}

		// compute candidates for each row
		Type kmid=0;
        for (i=0; i < A->rows; i++)
		{
			int counter = 0;
			for (j=0; j < A->columns; j++)
				if (A->m[i][j] <= min[i]+eps[i])
					counter++;
			k[i]= counter;
			kmid += counter;

			// alloc candidates
			assigncost[i] = new cost[k[i]];  // counter-1] ?
			minassign2[i] = new int[k[i]];   // counter-1] ?
			w[i] = new cost[k[i]];           // counter-1] ?

            // save candidate matching values
			counter =0;
			for (j=0; j < A->columns; j++)
				if (A->m[i][j] <= min[i]+eps[i])
				{
					assigncost[i][counter]=A->m[i][j];
					minassign2[i][counter]=j;
					counter++;
				}
		}
		kmid /= A->rows;
//		printf("\n kmid = %f",kmid);

		// compute weights for each candidate in each row
		for (i=0; i< A->rows; i++)
			for (j=0; j<k[i]; j++)
				// TODO: weights should sum up to 1
				w[i][j] = (min[i]+eps[i] - assigncost[i][j]) / eps[i];


        					
		cost sum = lapcost = (Type)0;

		for (i=0; i<A->rows; i++)
		{	
			sum=0;
			for (j=0; j<k[i]; j++)
			{
				sum+= w[i][j]*assigncost[i][j];
			}
			sum = sum / k[i];

			lapcost+=sum;



		}



        lapcost /= (Type) A->rows;


		return lapcost;
	}





};

/// \typedef GlobalMatchingf - A float GlobalMatching
typedef GlobalMatching<float> GlobalMatchingf;
/// \typedef GlobalMatchingd - A double GlobalMatching
typedef GlobalMatching<double> GlobalMatchingd;

#endif //GLOBALMATCHING_H





