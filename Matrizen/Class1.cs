using System;
using System.Collections.Generic;

namespace Matrices
{
    enum RorC { Row, Col};

    /// <summary>
    /// complex or real element of Matrix.
    /// </summary>
    public class Element
    {
        public double re;
        public double im;

        public static Element zero = new Element(0, 0);
        public static Element real1 = new Element(1);
        public static Element comp1 = new Element(1, 1);

        /// <summary>
        /// create complex element of Matrix.
        /// </summary>
        /// <param name="real"></param>
        /// <param name="imaginary"></param>
        public Element(double real, double imaginary)
        {
            re = real;
            im = imaginary;
        }

        /// <summary>
        /// create real element of Matrix.
        /// </summary>
        /// <param name="real"></param>
        public Element(double real)
        {
            re = real;
            im = 0;
        }

        public override string ToString()
        {
            if (im == 0)
            {
                return string.Format("{0}", re);
            }
            else
            {
                return string.Format("{0} + i * {1}", re, im);
            }
        }

        public static Element operator +(Element a)
            => a;
        public static Element operator -(Element a)
            => new Element(-a.re, -a.im);
        public static Element operator +(Element a, Element b)
            => new Element(a.re + b.re, a.im + b.im);
        public static Element operator +(double a, Element b)
            => new Element(a + b.re, b.im);
        public static Element operator +(Element a, double b)
            => new Element(a.re + b, a.im);
        public static Element operator -(Element a, Element b)
            => new Element(a.re - b.re, a.im - b.im);
        public static Element operator -(double a, Element b)
            => new Element(a - b.re, b.im);
        public static Element operator -(Element a, double b)
            => new Element(a.re - b, a.im);
        public static Element operator *(Element a, Element b)
            => new Element(a.re * b.re - a.im * b.im, a.re * b.im + b.re * a.im);
        public static Element operator *(double a, Element b)
            => new Element(a * b.re, a * b.im);
        public static Element operator *(Element a, double b)
            => new Element(a.re * b, b * a.im);
        public static Element operator /(Element a, Element b)
        {
            if(b.re == 0 && b.im == 0)
            {
                throw new DivideByZeroException();
            }
            return new Element((a.re * b.re + a.im * b.im) / (b.re * b.re + b.im * b.im), (b.re * a.im - a.re * b.im) / (b.re * b.re + b.im * b.im));
        }
        public static Element operator /(double a, Element b)
        {
            if (b.re == 0 && b.im == 0)
            {
                throw new DivideByZeroException();
            }
            return new Element((a * b.re) / (b.re * b.re + b.im * b.im), (a * b.im) / (b.re * b.re + b.im * b.im));
        }
        public static Element operator /(Element a, double b)
        {
            if (b == 0)
            {
                throw new DivideByZeroException();
            }
            return new Element(a.re / b, a.im / b);
        }
        public static bool operator ==(Element a, Element b)
        {
            if(a.re == b.re && a.im == b.im)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        public static bool operator !=(Element a, Element b)
        {
            if (a.re == b.re && a.im == b.im)
            {
                return false;
            }
            else
            {
                return true;
            }
        }
    }
    
    /// <summary>
    /// Matrix m with complex or real elements.
    /// </summary>
    public class Matrix
    {
        protected int rows;
        protected int cols;
        protected bool isQuadratic;
        protected Element[,] elements;
        public Matrix()
        {
            rows = 0;
            cols = 0;
            isQuadratic = false;
        }

        /// <summary>
        /// create Matrix m with complex array.
        /// </summary>
        /// <param name="ele"></param>
        public Matrix(Element[,] ele)
        {
            rows = ele.GetLength(0);
            cols = ele.GetLength(1);

            elements = ele;

            if (rows == cols)
            {
                isQuadratic = true;
            }
            else
            {
                isQuadratic = false;
            }
        }

        /// <summary>
        /// create special quadratic Matrix m.
        /// </summary>
        /// <param name="specialquadratic">0 => m with only 0 || 1 => m is idendity matrix
        /// </param>
        /// <param name="n"></param>
        public Matrix(string specialquadratic, int n, bool isComplex)
        {
            rows = n;
            cols = n;
            isQuadratic = true;

            switch (specialquadratic)
            {
                case "0":
                    elements = new Element[rows, cols];
                    for (int row = 0; row < rows; row++)
                    {
                        for (int col = 0; col < rows; col++)
                        {
                            elements[row, col] = new Element(0);
                        }
                    }
                    break;
                case "1":
                    elements = new Element[rows, cols];
                    for (int row = 0; row < rows; row++)
                    {
                        for (int col = 0; col < rows; col++)
                        {
                            if (row != col)
                            {
                                elements[row, col] = new Element(0);
                            }
                            else
                            {
                                if (isComplex)
                                {
                                    elements[row, col] = Element.comp1;
                                }
                                else
                                {
                                    elements[row, col] = Element.real1;
                                }
                            }
                        }
                    }
                    break;
            }
        }

        /// <summary>
        /// create special quadratic Matrix m with diagonals from element array. diagonals.Length = n.
        /// </summary>
        /// <param name="diagonals"></param>
        /// <param name="n"></param>
        public Matrix(Element[] diagonals, int n)
        {
            if (diagonals.Length != n)
            {
                throw new ArgumentException("ArrayLength has to be equal to n.");
            }

            rows = n;
            cols = n;
            isQuadratic = true;

            elements = new Element[rows, cols];
            for (int row = 0; row < rows; row++)
            {
                for (int col = 0; col < rows; col++)
                {
                    if (row == col)
                    {
                        elements[row, col] = diagonals[row];
                    }
                    else
                    {
                        elements[row, col] = new Element(0, 0);
                    }
                }
            }
        }

        /// <summary>
        /// create 0 Matrix m with i rows and k columns.
        /// </summary>
        /// <param name="i"></param>
        /// <param name="k"></param>
        public Matrix(int i, int k)
        {
            rows = i;
            cols = k;

            elements = new Element[rows, cols];

            for (int row = 0; row < rows; row++)
            {
                for (int col = 0; col < rows; col++)
                {
                    elements[row, col] = new Element(0, 0);
                }
            }

            if (rows == cols)
            {
                isQuadratic = true;
            }
            else
            {
                isQuadratic = false;
            }
        }

        /// <summary>
        /// create SubMatrix deleting i. row and k. column. (indexing starts with 1)
        /// </summary>
        /// <param name="i"></param>
        /// <param name="k"></param>
        /// <returns></returns>
        public Matrix SubMatrix(int i, int k)
        {
            Element[,] sub = new Element[rows - 1, cols - 1];
            int newrow = 0;
            int newcol;
            i--; k--;
            for (int row = 0; row < rows; row++)
            {
                if (row != i)
                {
                    newcol = 0;
                    for (int col = 0; col < cols; col++)
                    {
                        if (col != k)
                        {
                            sub[newrow, newcol] = elements[row, col];
                            newcol++;
                        }
                    }
                    newrow++;
                }
            }
            return new Matrix(sub);
        }

        /// <summary>
        /// calculate transpose of Matrix.
        /// </summary>
        /// <returns></returns>
        public Matrix Transpose()
        {
            Element[,] transpose = new Element[cols, rows];
            for(int row = 0; row<rows; row++)
            {
                for(int col = 0; col<cols; col++)
                {
                    transpose[col, row] = elements[row, col];
                }
            }
            return new Matrix(transpose);
        }

        /// <summary>
        /// calculate inverse of Matrix.
        /// </summary>
        /// <returns></returns>
        public Matrix Inverse()
        {
            if(!IsRegular())
            {
                throw new ArgumentException("cannot calculate inverse of non regular matrix.");
            }
            if(!isQuadratic)
            {
                throw new ArgumentException("cannot calculate inverse of non quadratic matrix.");
            }
            Determinant determinant = new Determinant(this);
            Element[,] inverse = new Element[rows, cols];
            for(int row = 1; row <= rows; row++)
            {
                for(int col = 1; col <= cols; col++)
                {
                    Determinant subdeterminant = new Determinant(this.SubMatrix(row, col));
                    Element adjunct;
                    if((row + col) % 2 == 0)
                    {
                        adjunct = subdeterminant.det;
                    }
                    else
                    {
                        adjunct = -(subdeterminant.det);
                    }
                    inverse[row - 1, col - 1] = adjunct / determinant.det;
                }
            }
            return new Matrix(inverse);
        }

        /// <summary>
        /// calculate conjugate of Matrix.
        /// </summary>
        /// <returns></returns>
        public Matrix Conjugate()
        {
            Element[,] conjugate = new Element[rows, cols];
            for (int row = 0; row < rows; row++)
            {
                for (int col = 0; col < cols; col++)
                {
                    conjugate[row, col] = new Element(elements[row, col].re, - elements[row, col].im);
                }
            }
            return new Matrix(conjugate);
        }

        /// <summary>
        /// calculate conjugate transpose of Matrix.
        /// </summary>
        /// <returns></returns>
        public Matrix ConjugateTranspose()
        {
            return this.Conjugate().Transpose();
        }

        /// <summary>
        /// calculate complex or real trace of Matrix.
        /// </summary>
        /// <returns></returns>
        public Element Trace()
        {
            if(!isQuadratic)
            {
                throw new ArgumentException("cannot calculate trace of non quadratic matrix.");
            }
            Element trace = new Element(0);
            for(int i = 0; i<rows; i++)
            {
                trace += elements[i, i];
            }
            return trace;
        }

        /// <summary>
        /// specify wether Matrix m is symmetric. 
        /// returns true if m = m.transpose.
        /// </summary>
        /// <returns></returns>
        public bool IsSymmetric()
        {
            if (isQuadratic)
            {
                if (this == this.Transpose())
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }

        /// <summary>
        /// specify wether Matrix m is skewsymmetric. 
        /// returns true if m = - m.transpose.
        /// </summary>
        /// <returns></returns>
        public bool IsSkewSymmetric()
        {
            if (isQuadratic)
            {
                if (this == -(this.Transpose()))
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }

        /// <summary>
        /// specify wether Matrix m is complex. 
        /// returns true if there is at least one complex element.
        /// </summary>
        /// <returns></returns>
        public bool IsComplexMatrix()
        {
            bool ret = false;
            foreach(var element in elements)
            {
                if(element.im != 0)
                {
                    ret = true;
                }
            }
            return ret;
        }

        /// <summary>
        /// specify wether Matrix m is orthogonal. 
        /// returns true if m * m.transpose = idendity matrix.
        /// </summary>
        /// <returns></returns>
        public bool IsOrthogonal()
        {
            Matrix E = new Matrix("1", this.rows, IsComplexMatrix());

            if (isQuadratic)
            {
                if (this * this.Transpose() == E)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }

        /// <summary>
        /// specify wether Matrix m is quadratic.
        /// </summary>
        /// <returns></returns>
        public bool IsQuadratic()
        {
            return isQuadratic;
        }

        /// <summary>
        /// specify wether Matrix m is regular. 
        /// returns true if det(m) != 0.
        /// </summary>
        /// <returns></returns>
        public bool IsRegular()
        {
            if(isQuadratic)
            {
                Determinant determinant = new Determinant(this);
                if(determinant.det != new Element(0))
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }

        /// <summary>
        /// specify wether Matrix m is hermitic. 
        /// returns true if m = m.conjugatetranspose.
        /// </summary>
        /// <returns></returns>
        public bool IsHermitic()
        {
            if (isQuadratic)
            {
                if (this == this.ConjugateTranspose())
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }

        /// <summary>
        /// specify wether Matrix m is skewhermitic. 
        /// returns true if m = - m.conjugatetranspose.
        /// </summary>
        /// <returns></returns>
        public bool IsSkewHermitic()
        {
            if (isQuadratic)
            {
                if (this == (-this.ConjugateTranspose()))
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }

        /// <summary>
        /// specify wether Matrix m is unitary. 
        /// returns true if m * m.conjugatetranspose = idendity matrix.
        /// </summary>
        /// <returns></returns>
        public bool IsUnitary()
        {
            if (isQuadratic)
            {
                if (this * this.ConjugateTranspose() == new Matrix("1", rows, IsComplexMatrix()))
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return false;
            }
        }

        /// <summary>
        /// returns row count of Matrix.
        /// </summary>
        /// <returns></returns>
        public int GetRows()
        {
            return rows;
        }
        /// <summary>
        /// returns column count of Matrix.
        /// </summary>
        /// <returns></returns>
        public int GetCols()
        {
            return cols;
        }
        /// <summary>
        /// returns element in (row, col). (Indexing starts at 1) 
        /// </summary>
        /// <param name="row"></param>
        /// <param name="col"></param>
        /// <returns></returns>
        public Element GetElement(int row, int col)
        {
            if(row > rows || col > cols || row < 1 || col < 1)
            {
                throw new IndexOutOfRangeException();
            }
            return elements[row-1, col-1];
        }

        //***
        //RankCalculation
        //
        /// <summary>
        /// calculate rank of Matrix m.
        /// </summary>
        /// <returns></returns>
        public int Rank()
        {
            if(IsRegular())
            {
                return rows;
            }
            else
            {
                return RankRecursion(this);
            }
        }

        private Matrix SubMatrix(int todelete, RorC roworcolumn)
        {
            Element[,] sub;
            int newrow;
            int newcol;
            switch(roworcolumn)
            {
                case RorC.Row:
                    sub = new Element[rows - 1, cols];
                    todelete--;
                    newrow = 0;
                    for (int row = 0; row < rows; row++)
                    {
                        if (row != todelete)
                        {
                            newcol = 0;
                            for (int col = 0; col < cols; col++)
                            {
                                if (true)
                                {
                                    sub[newrow, newcol] = elements[row, col];
                                    newcol++;
                                }
                            }
                            newrow++;
                        }
                    }
                    break;
                case RorC.Col:
                    sub = new Element[rows, cols - 1];
                    todelete--;
                    newrow = 0;
                    for (int row = 0; row < rows; row++)
                    {
                        if (true)
                        {
                            newcol = 0;
                            for (int col = 0; col < cols; col++)
                            {
                                if (col != todelete)
                                {
                                    sub[newrow, newcol] = elements[row, col];
                                    newcol++;
                                }
                            }
                            newrow++;
                        }
                    }
                    break;
                default:
                    throw new ArgumentException("row or column to delete not specified!");
            }
            return new Matrix(sub);
        }

        private int RankRecursion(Matrix m)
        {
            int r = m.GetRows();
            int c = m.GetCols();
            if(r == c)
            {
                if(m.IsRegular())
                {
                    return r;
                }
                else if(r == 2)
                {
                    if(m.GetElement(1,1) != Element.zero || m.GetElement(1,2) != Element.zero || m.GetElement(2,1) != Element.zero || m.GetElement(2,2) != Element.zero)
                    {
                        return 1;
                    }
                    else
                    {
                        return 0;
                    }
                }
                else
                {
                    List<Matrix> mlist = new List<Matrix>();
                    for(int i = 1; i<=r; i++)
                    {
                        for(int k = 1; k<=c; k++)
                        {
                            mlist.Add(m.SubMatrix(i, k));
                        }
                    }
                    List<int> ilist = new List<int>();
                    foreach(var element in mlist)
                    {
                        ilist.Add(RankRecursion(element));
                    }
                    return MaxInList(ilist);
                }
            }
            else if(r > c)
            {
                List<Matrix> mlist = new List<Matrix>();
                for(int i = 1; i<=r; i++)
                {
                    mlist.Add(m.SubMatrix(i, RorC.Row));
                }
                List<int> ilist = new List<int>();
                foreach(var element in mlist)
                {
                    ilist.Add(RankRecursion(element));
                }
                return MaxInList(ilist);
            }
            else
            {
                List<Matrix> mlist = new List<Matrix>();
                for (int i = 1; i <= c; i++)
                {
                    mlist.Add(m.SubMatrix(i, RorC.Col));
                }
                List<int> ilist = new List<int>();
                foreach (var element in mlist)
                {
                    ilist.Add(RankRecursion(element));
                }
                return MaxInList(ilist);
            }
        }

        private int MaxInList(List<int> ilist)
        {
            int maximum = int.MinValue;
            foreach(var element in ilist)
            {
                if(element > maximum)
                {
                    maximum = element;
                }
            }
            return maximum;
        }
        //
        //***

        public static Matrix operator +(Matrix A)
        {
            return A;
        }
        public static Matrix operator -(Matrix A)
        {
            Element[,] negativ = new Element[A.rows, A.cols];
            for(int row = 0; row<A.rows; row++)
            {
                for(int col = 0; col<A.cols; col++)
                {
                    negativ[row, col] = -A.elements[row, col];
                }
            }
            return new Matrix(negativ);
        }
        public static Matrix operator +(Matrix A, Matrix B)
        {
            if(A.rows != B.rows || A.cols != B.cols)
            {
                throw new ArgumentException("cannot add two Matrices of different type.");
            }
            Element[,] sum = new Element[A.rows, A.cols];
            for(int row = 0; row<A.rows; row++)
            {
                for(int col = 0; col<A.cols; col++)
                {
                    sum[row, col] = A.elements[row, col] + B.elements[row, col];
                }
            }
            return new Matrix(sum);
        }
        public static Matrix operator -(Matrix A, Matrix B)
        {
            if (A.rows != B.rows || A.cols != B.cols)
            {
                throw new ArgumentException("cannot add two Matrices of different type.");
            }
            Element[,] newele = new Element[A.rows, A.cols];
            for (int row = 0; row < A.rows; row++)
            {
                for (int col = 0; col < A.cols; col++)
                {
                    newele[row, col] = A.elements[row, col] - B.elements[row, col];
                }
            }
            return new Matrix(newele);
        }
        public static Matrix operator *(double A, Matrix B)
        {
            Element[,] newele = new Element[B.rows, B.cols];
            for (int row = 0; row < B.rows; row++)
            {
                for (int col = 0; col < B.cols; col++)
                {
                    newele[row, col] = A * B.elements[row, col];
                }
            }
            return new Matrix(newele);
        }
        public static Matrix operator *(Matrix A, Matrix B)
        {
            if(A.cols != B.rows)
            {
                throw new ArgumentException("Cannot multiply Matrices of that type.");
            }
            int n = A.cols;
            Element[,] matrixproduct = new Element[A.rows, B.cols];
            for(int row = 0; row<A.rows; row++)
            {
                for(int col = 0; col<B.cols; col++)
                {
                    Element rowsum = new Element(0);
                    for(int i = 0; i<n; i++)
                    {
                        rowsum += (A.elements[row, i] * B.elements[i, col]);
                    }
                    matrixproduct[row, col] = rowsum;
                }
            }
            return new Matrix(matrixproduct);
        }
        public static bool operator ==(Matrix A, Matrix B)
        {
            if(A.rows != B.rows || A.cols != B.cols)
            {
                return false;
            }
            bool ret = true;
            for(int row = 0; row<A.rows; row++)
            {
                for(int col = 0; col<A.cols; col++)
                {
                    if(A.elements[row, col] != B.elements[row, col])
                    {
                        ret = false;
                    }
                }
            }
            return ret;
        }
        public static bool operator !=(Matrix A, Matrix B)
        {
            if (A.rows != B.rows || A.cols != B.cols)
            {
                return true;
            }
            bool ret = false;
            for (int row = 0; row < A.rows; row++)
            {
                for (int col = 0; col < A.cols; col++)
                {
                    if (A.elements[row, col] != B.elements[row, col])
                    {
                        ret = true;
                    }
                }
            }
            return ret;
        }
    }

    /// <summary>
    /// Determinant of Matrix m.
    /// </summary>
    public class Determinant
    {
        /// <summary>
        /// real or complex value of determinant of Matrix m.
        /// </summary>
        public Element det;

        private Determinant(Matrix matrix, int row, int col)
        {
            det = Calculate(matrix.SubMatrix(row, col));
        }

        /// <summary>
        /// calculate determinant of Matrix m.
        /// </summary>
        /// <param name="matrix"></param>
        public Determinant(Matrix matrix)
        {
            det = Calculate(matrix);
        }

        private Element Calculate(Matrix m)
        {
            if(m.IsQuadratic())
            {
                if(m.GetRows() == 2)
                {
                    return (m.GetElement(1, 1) * m.GetElement(2, 2) - m.GetElement(1, 2) * m.GetElement(2, 1));
                }
                else
                {
                    Element sum = new Element(0, 0);
                    for(int col = 1; col <= m.GetCols(); col++)
                    {
                        Determinant sub = new Determinant(m, 1, col);
                        if(col + 1 % 2 == 0)
                        {
                            sum += m.GetElement(1, col) * sub.det;
                        }
                        else
                        {
                            sum += m.GetElement(1, col) * (-sub.det);
                        }
                    }
                    return sum;
                }
            }
            else
            {
                throw new ArgumentException("cannot calculate determinant of non quadratic matrix.");
            }
        }
    }
}
