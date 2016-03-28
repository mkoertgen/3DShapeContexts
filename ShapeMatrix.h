/** 
* \author Marcel Koertgen (marcel.koertgen@gmail.com)
* \date 02-02-2003
* \version 1.0
*/
#ifndef SHAPEMATRIX_H
#define SHAPEMATRIX_H

/**
* \class ShapeMatrix
* \brief template ShapeMatrix class
*
* This class is a template ShapeMatrix class used as a basis for any global
* matching method.
* The class ShapeMatrix is used as a (quadratic) matrix saving the values
* from the former local matching.
* Once these values are calculated a highly discrimanitive global matching
* would try to minimize the assignment-costs for all <i,j> entries.
* This could be achieved by the hungarian method: O(n³) or
* clustering: O(n²logn).
* One should think over a randomized or Monte-Carlo method that would save
* lots of time in an early comparison where many candidates are available.
* After reducing the number of candidates to an adaptive/iterative threshold
* the more expensive matching methods can be applied.
* Another idea is to use soft assignments instead of hard 1-1 correspondences.
* Therefore the ShapeMatrix is not constrained to be square.
*/
template<class Type> class ShapeMatrix
{
	public:
		Type **m;    /// matrix elements
		int rows;    /// number of rows for the matrix
		int columns; /// number of rows for the matrix

	/// constructor, sets all entries 0
	ShapeMatrix(int rows, int columns) 
	{
		allocEntries(rows,columns);
	}

	/// destructor
	~ShapeMatrix()
	{
		deleteEntries();
	}

	/// Copy-constructor
	ShapeMatrix(const ShapeMatrix& in) 
	{ 
		if (rows!= in.rows || columns!=in.columns)
		{
			deleteEntries();
			allocEntries(in.rows, in.columns);
		}
		int i,j;
		for (i=0; i<rows; i++)
		{
			for (j=0; j<columns; j++)
				m[i][j] = in.m[i][j];
		}
	}

	/// Atribuition operator
	void operator=(const ShapeMatrix& in) 
	{ 
		if (rows!= in.rows || columns!=in.columns)
		{
			deleteEntries();
			allocEntries(in.rows, in.columns);
		}
		int i,j;
		for (i=0; i<rows; i++)
		{
			for (j=0; j<columns; j++)
				m[i][j] = in.m[i][j];
		}
	}

	/// Nullify all elements
	inline void null()
	{
		int i,j;
		for (i=0; i<rows; i++)
		{
			if (m[i]) {
				for (j=0; j<columns; j++) m[i][j] = (Type)0.0;
			}
		}
	}

	/// print matrix elements (uses printf)
	inline void print()
	{
		int i,j;
		for (i=0; i<rows; i++)
			if (m[i])
			{
				printf("[ %.2f",m[i][0]);
				for (j=1; j<columns; j++)
					printf(", %.2f", m[i][j]);
				printf(" ]\n");
			} else printf("NULL\n");
	}

	/// print to ostream
	inline void Print(std::ostream &output) const 
	{  
		int i,j;
		for (i=0; i<rows; i++)
			if (m[i])
			{
				output << "[ " << m[i][0];
				for (j=1; j<columns; j++)
					output << ", " << m[i][j];
				output << " ]" << std::endl;
			}
			else output << "NULL" << std::endl;
	}

	/// ostream "<<" operator
	friend inline std::ostream & operator << (std::ostream& output, const ShapeMatrix<Type> &s)
	{
	  s.Print(output);
	  return output; 
	}

private:
	void deleteEntries()
	{
		int i;
		for (i=0; i<rows; i++)
			if (m[i]) {
				delete [] m[i];
			}
	}

	void allocEntries(int rows, int columns)
	{
		if (rows>0 &&
			columns>0)
		{
			this->rows = rows;
			this->columns = columns;
			m = new Type *[rows];
			int i;
			for (i=0; i < rows; i++) m[i] = new Type[columns];
			this->null();
		}
	}

};

/// \typedef ShapeMatrixf - A float ShapeMatrix
typedef ShapeMatrix<float> ShapeMatrixf;

/// \typedef ShapeMatrixd - A double ShapeMatrix
typedef ShapeMatrix<double> ShapeMatrixd;

#endif //SHAPEMATRIX_H
