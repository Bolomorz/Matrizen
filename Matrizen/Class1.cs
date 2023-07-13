using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO.Compression;
using System.Threading;
using System.Threading.Tasks;

namespace Matrices
{
    public enum RorC { Row, Col};

    /// <summary>
    /// complex or real element of Matrix.
    /// </summary>
    public class Element
    {
        public double re;
        public double im;

        public bool isNull;

        public static Element zero = new Element(0, 0);
        public static Element one = new Element(1);
        public static Element NULL = new Element(true);
        public static Element tol = new Element(double.Epsilon);

        /// <summary>
        /// create complex element of Matrix.
        /// </summary>
        /// <param name="real"></param>
        /// <param name="imaginary"></param>
        public Element(double real, double imaginary)
        {
            re = real;
            im = imaginary;
            isNull = false;
        }

        /// <summary>
        /// create real element of Matrix.
        /// </summary>
        /// <param name="real"></param>
        public Element(double real)
        {
            re = real;
            im = 0;
            isNull = false;
        }

        protected Element(bool isnull)
        {
            re = 0;
            im = 0;
            isNull = isnull;
        }

        public Element SquareRoot()
        {
            if (im == 0)
            {
                return new Element(Math.Sqrt(re));
            }
            else
            {
                double real = Math.Sqrt((re + Math.Sqrt(re * re + im * im)) / 2);
                double imag = (im / Math.Abs(im)) * Math.Sqrt((-re + Math.Sqrt(re * re + im * im)) / 2);
                return new Element(real, imag);
            }
        }

        public Element SIGN()
        {
            return new Element(re / this.ABS(), im / this.ABS());
        }

        public double ABS()
        {
            return Math.Sqrt(re * re + im * im);
        }

        public override string ToString()
        {
            if (im == 0)
            {
                return string.Format("{0}", Math.Round(re, 5));
            }
            else
            {
                return string.Format("{0} + i * {1}", Math.Round(re, 5), Math.Round(im, 5));
            }
        }

        public bool IsNull()
        {
            return isNull;
        }

        public override int GetHashCode()
        {
            return Convert.ToInt32(re + im);
        }
        public override bool Equals(object obj)
        {
            if(obj == null || !(obj is Element))
            {
                return false;
            }
            return this == (Element)obj;
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
        public static bool operator >=(Element a, Element b)
        {
            if(a.re >= b.re)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        public static bool operator <=(Element a, Element b)
        {
            if (a.re <= b.re)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        public static bool operator >(Element a, Element b)
        {
            if(a.re > b.re)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        public static bool operator <(Element a, Element b)
        {
            if (a.re < b.re)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    }
    
    /// <summary>
    /// Matrix m with complex or real elements.
    /// </summary>
    public class Matrix/**/
    {
        protected int rows;
        protected int cols;
        protected bool isQuadratic;
        protected Element[,] elements;
        protected Element det;
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
            det = Element.NULL;
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
                                    elements[row, col] = Element.one;
                                }
                                else
                                {
                                    elements[row, col] = Element.one;
                                }
                            }
                        }
                    }
                    break;
            }
            det = Element.NULL;
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
            det = Element.NULL;
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
            if(det.isNull)
            {
                Determinant d = new Determinant(this);
                det = d.det;
            }
            Element[,] inverse = new Element[rows, cols];
            for(int row = 1; row <= rows; row++)
            {
                for(int col = 1; col <= cols; col++)
                {
                    Element subdeterminant = this.SubMatrix(row, col).GetDeterminant();
                    Element adjunct;
                    if((row + col) % 2 == 0)
                    {
                        adjunct = subdeterminant;
                    }
                    else
                    {
                        adjunct = -(subdeterminant);
                    }
                    inverse[row - 1, col - 1] = adjunct / det;
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
        /// calculate eigenvalues of Matrix.
        /// </summary>
        /// <returns></returns>
        public Eigenvalue[] Eigenvalues()
        {
            EigenvalueProblem evp = new EigenvalueProblem(this);
            return evp.GetEigenvalues();
        }

        /// <summary>
        /// specify if Matrix m is symmetric. 
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
        /// specify if Matrix m is skewsymmetric. 
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
        /// specify if Matrix m is complex. 
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
        /// specify if Matrix m is orthogonal. 
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
        /// specify if Matrix m is quadratic.
        /// </summary>
        /// <returns></returns>
        public bool IsQuadratic()
        {
            return isQuadratic;
        }

        /// <summary>
        /// specify if Matrix m is regular. 
        /// returns true if det(m) != 0.
        /// </summary>
        /// <returns></returns>
        public bool IsRegular()
        {
            if(isQuadratic)
            {
                if(det.isNull)
                {
                    Determinant d = new Determinant(this);
                    det = d.det;
                }
                if(det != Element.zero)
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
        /// specify if Matrix m is hermitic. 
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
        /// specify if Matrix m is skewhermitic. 
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
        /// specify if Matrix m is unitary. 
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

        public Element[,] GetElements()
        {
            return elements;
        }

        public Element GetDeterminant()
        {
            if(det.isNull)
            {
                Determinant d = new Determinant(this);
                det = d.det;
            }
            return det;
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

        public override string ToString()
        {
            string ret = "";
            for(int row = 0; row<rows; row++)
            {
                string r = "row" + (row + 1) + ":\t";
                for(int col = 0; col<cols; col++)
                {
                    r += " " + elements[row, col] + " |";
                }
                if(row == 0)
                {
                    ret += r;
                }
                else
                {
                    ret += Environment.NewLine + r;
                }
            }
            return ret;
        }

        public override int GetHashCode()
        {
            int sum = 0;
            foreach(var element in elements)
            {
                sum += element.GetHashCode();
            }
            return sum;
        }
        public override bool Equals(object obj)
        {
            if(obj == null || !(obj is Matrix))
            {
                return false;
            }
            return this == (Matrix)obj;
        }
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
            for(int row = 0; row<A.rows; row++)
            {
                for(int col = 0; col<A.cols; col++)
                {
                    if(A.elements[row, col] != B.elements[row, col])
                    {
                        return false;
                    }
                }
            }
            return true;
        }
        public static bool operator !=(Matrix A, Matrix B)
        {
            if (A.rows != B.rows || A.cols != B.cols)
            {
                return true;
            }
            for (int row = 0; row < A.rows; row++)
            {
                for (int col = 0; col < A.cols; col++)
                {
                    if (A.elements[row, col] != B.elements[row, col])
                    {
                        return true;
                    }
                }
            }
            return false;
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

        /// <summary>
        /// calculate determinant of Matrix m.
        /// </summary>
        /// <param name="matrix"></param>
        public Determinant(Matrix matrix)
        {
            if (!matrix.IsQuadratic())
            {
                throw new ArgumentException("cannot calculate determinant of non quadratic matrix.");
            }
            Tuple<int, Element[,], int[]> LUPdec = LUPDecompose(matrix.GetElements(), matrix.GetRows(), new Element(double.Epsilon));
            if (LUPdec.Item1 == 1)
            {
                det = LUPDeterminant(LUPdec.Item2, LUPdec.Item3, matrix.GetCols());
            }
            else
            {
                det = Calculate(matrix.GetElements(), matrix.GetRows());
            }
        }

        private Element LUPDeterminant(Element[,] A, int[] P, int n)
        {
            Element d = A[0, 0];

            for(int i = 1; i < n; i++)
            {
                d *= A[i, i];
            }

            if((P[n] - n) % 2 == 0)
            {
                return d;
            }
            else
            {
                return -d;
            }
        }

        private Tuple<int, Element[,], int[]> LUPDecompose(Element[,] A, int n, Element Tol)
        {
            int i, j, k, imax;
            Element maxA, absA;
            Element[] ptr = new Element[n];

            int[] P = new int[n+1];
            Element[,] decompose = new Element[n, n];
            
            for(i = 0; i < n; i++)
            {
                for(j = 0; j < n; j++)
                {
                    decompose[i, j] = A[i, j];
                }
            }

            for(i = 0; i <= n; i++)
            {
                P[i] = i;
            }

            for(i = 0; i < n; i++)
            {
                maxA = Element.zero;
                imax = i;

                for(k = i; k < n; k++)
                {
                    absA = new Element(decompose[k, i].ABS());
                    if(absA > maxA)
                    {
                        maxA = absA;
                        imax = k;
                    }
                }

                if (maxA < Tol) return new Tuple<int, Element[,], int[]>(0, decompose, P);

                if(imax != i)
                {
                    j = P[i];
                    P[i] = P[imax];
                    P[imax] = j;

                    for(k = 0; k < n; k++)
                    {
                        ptr[k] = decompose[i, k];
                    }
                    for (k = 0; k < n; k++)
                    {
                        decompose[i, k] = decompose[imax, k];
                    }
                    for (k = 0; k < n; k++)
                    {
                        decompose[imax, k] = ptr[k];
                    }

                    P[n]++;
                }

                for(j = i + 1; j < n; j++)
                {
                    decompose[j, i] /= decompose[i, i];

                    for(k = i + 1; k < n; k++)
                    {
                        decompose[j, k] -= decompose[j, i] * decompose[i, k];
                    }
                }
            }
            return new Tuple<int, Element[,], int[]>(1, decompose, P);
        }

        private Element[,] Sub(Element[,] m, int i, int k, int n)
        {
            Element[,] sub = new Element[n - 1, n - 1];
            int newrow = 0;
            int newcol;
            i--; k--;
            for (int row = 0; row < n; row++)
            {
                if (row != i)
                {
                    newcol = 0;
                    for (int col = 0; col < n; col++)
                    {
                        if (col != k)
                        {
                            sub[newrow, newcol] = m[row, col];
                            newcol++;
                        }
                    }
                    newrow++;
                }
            }
            return sub;
        }

        private Element Calculate(Element[,] m, int n)
        {
                if(n == 2)
                {
                    return (m[0, 0] * m[1, 1] - m[0, 1] * m[1, 0]);
                }
                else
                {
                    Element sum = new Element(0, 0);
                    List<Element> dets = new List<Element>();
                    for(int i = 1; i <= n; i++)
                    { 
                        dets.Add(Calculate(Sub(m, 1, i, n), n-1)); 
                    }
                    //Parallel.For(1, n + 1, index => dets.Add(Calculate(Sub(m, 1, index, n), n - 1)));
                    for (int col = 1; col <= n; col++)
                    {
                        Element sub = dets[col-1];
                        if ((col + 1) % 2 == 0)
                        {
                            sum += m[0, col-1] * sub;
                        }
                        else
                        {
                            sum += -1 * m[0, col-1] * sub;
                        }
                    }
                    return sum;
                }
        }
    }

    /// <summary>
    /// Characteristic polynomial of Matrix m.
    /// </summary>
    public class CharacteristicPolynomial
    {
        protected Element[] coefficients;
        protected Tuple<Element, Element>[,] characteristicmatrix;

        /// <summary>
        /// calculate characteristic polynomial of Matrix m.
        /// </summary>
        /// <param name="m"></param>
        public CharacteristicPolynomial(Matrix m)
        {
            if(!m.IsQuadratic())
            {
                throw new ArgumentException("cannot calculate characteristic polynomial of non quadratic matrix.");
            }
            int n = m.GetRows();
            CalculateCharacteristicMatrix(m, n);
            coefficients = CharacteristicPolynomialRecursion(characteristicmatrix, n);
        }

        private void CalculateCharacteristicMatrix(Matrix m, int n)
        {
            characteristicmatrix = new Tuple<Element, Element>[n, n];
            for(int i = 0; i<n; i++)
            {
                for(int j = 0; j<n; j++)
                {
                    Element element = m.GetElement(i + 1, j + 1);
                    if(i == j)
                    {
                        characteristicmatrix[i, j] = new Tuple<Element, Element>(element, new Element(-1));
                    }
                    else
                    {
                        characteristicmatrix[i, j] = new Tuple<Element, Element>(element, Element.zero);
                    }
                }
            }
        }

        private Element[] CharacteristicPolynomialRecursion(Tuple<Element, Element>[,] sub, int n)
        {
            Element[] ret;
            if (n == 2)
            {
                if (sub[0, 0].Item2 != Element.zero && sub[1, 1].Item2 != Element.zero)
                {
                    ret = new Element[3];
                    ret[0] = sub[0, 0].Item1 * sub[1, 1].Item1 - sub[0, 1].Item1 * sub[1, 0].Item1;
                    ret[1] = -(sub[0, 0].Item1) - sub[1, 1].Item1;
                    ret[2] = Element.one;
                }
                else if(sub[0, 0].Item2 != Element.zero)
                {
                    ret = new Element[2];
                    ret[0] = sub[0, 0].Item1 * sub[1, 1].Item1 - sub[0, 1].Item1 * sub[1, 0].Item1;
                    ret[1] = -(sub[1, 1].Item1);
                }
                else if (sub[1, 1].Item2 != Element.zero)
                {
                    ret = new Element[2];
                    ret[0] = sub[0, 0].Item1 * sub[1, 1].Item1 - sub[0, 1].Item1 * sub[1, 0].Item1;
                    ret[1] = -(sub[0, 0].Item1);
                }
                else if (sub[1, 0].Item2 != Element.zero)
                {
                    ret = new Element[2];
                    ret[0] = sub[0, 0].Item1 * sub[1, 1].Item1 - sub[0, 1].Item1 * sub[1, 0].Item1;
                    ret[1] = -(sub[0, 1].Item1);
                }
                else if (sub[0, 1].Item2 != Element.zero)
                {
                    ret = new Element[2];
                    ret[0] = sub[0, 0].Item1 * sub[1, 1].Item1 - sub[0, 1].Item1 * sub[1, 0].Item1;
                    ret[1] = -(sub[1, 0].Item1);
                }
                else
                {
                    ret = new Element[1];
                    ret[0] = sub[0, 0].Item1 * sub[1, 1].Item1 - sub[0, 1].Item1 * sub[1, 0].Item1;
                }
            }
            else
            {
                List<Element[]> subpolynomial = new List<Element[]>();
                for(int col = 0; col<n; col++)
                {
                    Tuple<Element, Element>[,] submatrix = CalculateSubMatrix(sub, n, col);
                    Tuple<Element, Element> adj;
                    if ((col + 1) % 2 == 1)
                    {
                        adj = sub[0, col];
                    }
                    else
                    {
                        adj = new Tuple<Element, Element>(-(sub[0, col].Item1), -(sub[0,col].Item2));
                    }
                    subpolynomial.Add(CalculateAdjunctPolynomial(CharacteristicPolynomialRecursion(submatrix, n - 1), adj));
                }
                ret = CalculatePolynomialSum(subpolynomial);
            }
            return ret;
        }

        private Element[] CalculatePolynomialSum(List<Element[]> subs)
        {
            int n = int.MinValue;
            foreach(var element in subs)
            {
                if(element.Length > n)
                {
                    n = element.Length;
                }
            }
            Element[] ret = new Element[n];
            for(int i = 0; i<n; i++)
            {
                ret[i] = new Element(0);
                foreach(var element in subs)
                {
                    if(i < element.Length)
                    {
                        ret[i] += element[i];
                    }
                }
            }
            return ret;
        }

        private Tuple<Element, Element>[,] CalculateSubMatrix(Tuple<Element, Element>[,] sub, int n, int column)
        {
            Tuple<Element, Element>[,] ret = new Tuple<Element, Element>[n - 1, n - 1];
            int newrow = 0;
            int newcol;
            for (int row = 0; row < n; row++)
            {
                if (row != 0)
                {
                    newcol = 0;
                    for (int col = 0; col < n; col++)
                    {
                        if (col != column)
                        {
                            ret[newrow, newcol] = sub[row, col];
                            newcol++;
                        }
                    }
                    newrow++;
                }
            }
            return ret;
        }

        private Element[] CalculateAdjunctPolynomial(Element[] sub, Tuple<Element, Element> adj)
        {
            int n = sub.Length;
            if (adj.Item2 != Element.zero)
            {
                Element[] ret = new Element[n + 1];
                ret[0] = sub[0] * adj.Item1;
                ret[n] = sub[n - 1] * adj.Item2;
                for(int i = 1; i<n; i++)
                {
                    ret[i] = sub[i] * adj.Item1 + sub[i - 1] * adj.Item2;
                }
                return ret;
            }
            else
            {
                Element[] ret = new Element[n];
                for(int i = 0; i<n; i++)
                {
                    ret[i] = sub[i] * adj.Item1;
                }
                return ret;
            }
        }

        /// <summary>
        /// get characteristic polynomial of Matrix m.
        /// returnspolynomial as array of elements P(x) = arr[0] + arr[1] * x + ....
        /// </summary>
        /// <returns></returns>
        public Element[] GetCharacteristicPolynomial()
        {
            return coefficients;
        }

    }

    /// <summary>
    /// Hessenbergtransform of Matrix m.
    /// </summary>
    public class HessenbergTransform
    {
        protected Matrix hessenbergtransform;

        public HessenbergTransform(Matrix m)
        {
            if(!m.IsQuadratic())
            {
                throw new ArgumentException("cannot calculate HessenbergTransform of non quadratic matrix.");
            }
            hessenbergtransform = CalculateHessenbergTransform(m);
        }

        private Matrix CalculateHessenbergTransform(Matrix m)
        {
            int n = m.GetRows();
            Element[,] hbtransform = new Element[n, n];
            for(int i = 1; i <= n; i++)
            {
                for(int j = 1; j <= n; j++)
                {
                    hbtransform[i - 1, j - 1] = m.GetElement(i, j);
                }
            }
            for(int j = 1; j <= n-2; j++)
            {
                for(int i = j + 2; i <= n; i++)
                {
                    Element aij = hbtransform[i - 1, j - 1];
                    Element aj1j = hbtransform[j, j - 1];
                    if(aij != Element.zero)
                    {
                        Element w, c, s;
                        if(aj1j.ABS() < double.Epsilon * aij.ABS())
                        {
                            w = -m.GetElement(i, j);
                            c = Element.zero;
                            s = Element.one;
                        }
                        else
                        {
                            Element det = aj1j * aj1j + aij * aij;
                            w = aj1j.SIGN() * det.SquareRoot();
                            c = aj1j / w;
                            s = -(aij / w);
                        }
                        hbtransform[j, j - 1] = w;
                        hbtransform[i - 1, j - 1] = Element.zero;
                        for(int k = j + 1; k <= n; k++)
                        {
                            Element h = c * hbtransform[j, k - 1] - s * hbtransform[i - 1, k - 1];
                            hbtransform[i - 1, k - 1] = s * hbtransform[j, k - 1] + c * hbtransform[i - 1, k - 1];
                            hbtransform[j, k - 1] = h;
                        }
                        for(int k = 1; k <= n; k++)
                        {
                            Element h = c * hbtransform[k - 1, j] - s * hbtransform[k - 1, i - 1];
                            hbtransform[k - 1, i - 1] = s * hbtransform[k - 1, j] + c * hbtransform[k - 1, i - 1];
                            hbtransform[k - 1, j] = h; 
                        }
                    }
                }
            }
            return new Matrix(hbtransform);
        }

        public Matrix GetHessenbergTransform()
        {
            return hessenbergtransform;
        }
    }

    public class RotationMatrix
    {
        protected Matrix rotationmatrix;
        public RotationMatrix(int p, int q, int n, Element x1, Element x2)
        {
            Element sin, cos;
            if(x1 == Element.zero)
            {
                sin = Element.one;
                cos = Element.zero;
            }
            else
            {
                Element tan = x2/x1;
                sin = tan / (1.0 + tan*tan).SquareRoot();
                cos = 1.0 / (1.0 + tan*tan).SquareRoot();
            }
            Element[,] rot = new Element[n,n];
            for(int i = 0; i < n; i++)
            {
                for(int j = 0; j < n; j++)
                {
                    if(i == j)
                    {
                        rot[i,j] = Element.one;
                    }
                    else
                    {
                        rot[i,j] = Element.zero;
                    }
                }
            }
            rot[p,q] = sin;
            rot[q,p] = - sin;
            rot[p,p] = cos;
            rot[q,q] = cos;
            rotationmatrix = new Matrix(rot);
        }

        public Matrix GetRotationMatrix()
        {
            return rotationmatrix;
        }
    }
    
    public class QRTransform
    {
        protected Element[,] qrtransform;

        public QRTransform(HessenbergTransform H)
        {
            int n = H.GetHessenbergTransform().GetRows();
            qrtransform = new Element[n, n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    qrtransform[i, j] = H.GetHessenbergTransform().GetElements()[i, j];
                }
            }
            //QR(qrtransform, n);
        }

        private void QRDoubleStep(Element[,] A, int n)
        {
            Element[,] C = new Element[2,2];
            C[0,0] = A[n-2,n-2];
            C[1,1] = A[n-1,n-1];
            C[1,0] = A[n-1,n-2];
            C[0,1] = A[n-2,n-1];
            Tuple<Element, Element> shifts = CalcEigenvals22(C);

            Element s = shifts.Item1 + shifts.Item2;
            Element t = shifts.Item1 * shifts.Item2;

            Element x1 = A[0,0]*A[0,0] + A[0,1]*A[1,0] - s * A[0,0] + t;
            Element x2 = A[1,0] * (A[0,0] + A[1,1] - s);
            Element x3 = A[1,0] * A[2,1];

            RotationMatrix U1 = new RotationMatrix(1, 2, n, x1, x2);

            Matrix B1 = U1.GetRotationMatrix().Transpose() * new Matrix(A) * U1.GetRotationMatrix();
            Element[,] EB1 = B1.GetElements();

            x1 = EB1[0,0] * EB1[0,0] + EB1[0,1]*EB1[1,0] - s * EB1[0,0] + t;

            RotationMatrix U2 = new RotationMatrix(1, 3, n, x1, x3);

            Matrix B = U2.GetRotationMatrix().Transpose() * new Matrix(A) * U2.GetRotationMatrix();
            Element[,] EB = B.GetElements();

            for(int i = 1; i < n-2; i++)
            {
                if(i == n - 3)
                {
                    RotationMatrix Ux = new RotationMatrix(i, i+1, n, EB[i,i], EB[i+1,i]);
                    B = Ux.GetRotationMatrix().Transpose() * B * Ux.GetRotationMatrix();
                }
                else
                {
                    RotationMatrix Ux1 = new RotationMatrix(i, i+1, n, EB[i,i], EB[i+1,i]);
                    RotationMatrix Ux2 = new RotationMatrix(i, i+2, n, EB[i,i], EB[i+2,i]);
                    B = Ux1.GetRotationMatrix().Transpose() * B * Ux1.GetRotationMatrix();
                    B = Ux2.GetRotationMatrix().Transpose() * B * Ux2.GetRotationMatrix();
                }
            }
        }

        private void QRStep(Element[,] A, Element shift, int n)
        {
            A[0,0] = A[0,0] - shift;
            Element w, c = Element.zero, s = Element.zero, ct = Element.zero, st = Element.zero;
            for(int i = 0; i < n; i++)
            {
                if(i < n-1)
                {
                    if(A[i,i].ABS() < double.Epsilon * A[i+1, i].ABS())
                    {
                        w = new Element(A[i+1, i].ABS());
                        c = Element.zero;
                        s = A[i+1, i].SIGN();
                    }
                    else
                    {
                        w = (A[i,i]*A[i,i] + A[i+1,i]*A[i+1,i]).SquareRoot();
                        c = A[i,i] / w;
                        s = - (A[i+1,i] / w);
                    }
                    for(int j = i+1; j < n; j++)
                    {
                        Element g = c * A[i,j] - s * A[i+1,j];
                        A[i+1, j] = s * A[i,j] + c * A[i+1,j];
                        A[i,j] = g;
                    }
                }
                
                if(i > 0)
                {
                    for(int j = 0; j < n; j++)
                    {
                        Element g = ct * A[j,i-1] - st * A[j,i];
                        A[j,i] = st * A[j,i-1] + ct * A[j,i];
                        A[j,i-1] = g;
                    }
                    A[i-1,i-1] = A[i-1,i-1] + shift;
                }

                ct = c;
                st = s;
            }
            A[n-1,n-1] = A[n-1,n-1] + shift;
        }

        private Tuple<Element, Element> CalcEigenvals22(Element[,] A)
        {
            Element b = -A[0,0] - A[1,1];
            Element c = A[0,0] * A[1,1] - A[0,1] * A[1,0];
            
            Element x = 0.5 * b;
            Element y = new Element(b.re*b.re - b.im*b.im-4*c.re, 2*b.im*b.re - 4*c.im);

            Element ev1 = x + y.SquareRoot();
            Element ev2 = x - y.SquareRoot();

            return new Tuple<Element, Element>(ev1, ev2);
        }

        public Matrix GetQRTransform()
        {
            return new Matrix(qrtransform);
        }
    }

    public class Eigenvalue
    {
        protected Element eigenvalue;
        protected Element[] eigenvector;

        public Eigenvalue(Matrix m, Element ev)
        {
            if(!m.IsQuadratic())
            {
                throw new ArgumentException("cannot calculate eigenvalues of non quadratic matrix.");
            }
            eigenvalue = ev;
            Element[,] evam = CalculateAugmentedMatrix(m, ev);
            eigenvector = Gauss(evam);
        }

        private Element[,] CalculateAugmentedMatrix(Matrix m, Element ev)
        {
            int n = m.GetRows();
            Element[,] evam = new Element[n, n+1];
            for(int i = 0; i<n; i++)
            {
                for(int j = 0; j<n; j++)
                {
                    if(i == j)
                    {
                        evam[i, j] = m.GetElement(i + 1, j + 1) - ev;
                    }
                    else
                    {
                        evam[i, j] = m.GetElement(i + 1, j + 1);
                    }
                }
                evam[i, n] = Element.zero;
            }
            return evam;
        }

        private Element[] Gauss(Element[,] am)
        {
            for(int line1 = 0; line1 < am.GetLength(0); line1++)
            {

                if(am[line1, line1] == Element.zero)
                {
                    for(int line2 = line1 + 1; line2 < am.GetLength(0); line2++)
                    {
                        if(am[line2, line1] != Element.zero)
                        {
                            for (int x = 0; x < am.GetLength(1); x++)
                            {
                                Element temp = am[x, line1];
                                am[x, line1] = am[x, line2];
                                am[x, line2] = temp;
                            }
                            break;
                        }
                    }
                }
                for(int line2 = line1 + 1; line2 < am.GetLength(0); line2++)
                {
                    if(am[line2, line1] == Element.zero)
                    {
                        continue;
                    }
                    Element factor = am[line2, line1] / am[line1, line1];
                    for(int x = line1; x < am.GetLength(1); x++)
                    {
                        am[line2, x] -= factor * am[line1, x];
                    }
                }
            }

            Element[] result = new Element[am.GetLength(0)];
            for(int line = am.GetLength(0) - 1; line >= 0; line--)
            {
                result[line] = am[line, am.GetLength(1) - 1];
                for(int x = line + 1; x < result.Length; x++)
                {
                    result[line] -= result[x] * am[line, x];
                }
                result[line] /= am[line, line];
            }
            return result;
        }

        public Element GetEigenvalue()
        {
            return eigenvalue;
        }

        public Element[] GetEigenvector()
        {
            return eigenvector;
        }
    }

    /// <summary>
    /// Eigenvalue problem of Matrix m.
    /// </summary>
    public class EigenvalueProblem
    {
        protected Eigenvalue[] eigenvalues;

        /// <summary>
        /// calculate eigenvalues and their eigenvectors for Matrix m.
        /// 1. calculate characteristic polynomial.
        /// 2. solve for eigenvalues with BairstowAlgorithm.
        /// </summary>
        /// <param name="m"></param>
        public EigenvalueProblem(Matrix m)
        {
            if (!m.IsQuadratic())
            {
                throw new ArgumentException("cannot calculate characteristic polynomial of non quadratic matrix.");
            }
            HessenbergTransform hbt = new HessenbergTransform(m);
            Matrix H = hbt.GetHessenbergTransform();
        }

        //
        //*
        //Bairstow Algorithm to solve characteristic polynomial
        //returns list of eigenvalues
        private List<Element> Bairstow(Element[] A)
        {
            List<Element> roots = new List<Element>();

            while(A.Length > 2)
            {
                int n = A.Length;
                Element u = A[n - 2] / A[n - 1];
                Element v = A[n - 3] / A[n - 1];

                Tuple<Element, Element, Element[]> quadfactor = QuadFactor(A, u, v);
                Tuple<Element, Element> quadraticroots = SolveQuadraticEquation(Element.one, quadfactor.Item1, quadfactor.Item2);
                roots.Add(quadraticroots.Item1);
                roots.Add(quadraticroots.Item2);
                A = quadfactor.Item3;
            }

            if(A.Length == 2)
            {
                roots.Add(-(A[0]) / A[1]);
            }

            return roots;
        }

        private Tuple<Element, Element, Element[]> QuadFactor(Element[] A, Element u, Element v)
        {
            const int ITERMAX = 1000;
            Element CONVTOL = new Element(0.0000000000000001);

            int n = A.Length - 1;
            Element c = Element.one;
            Element d = Element.one;
            int iter = 0;

            Element[] B = new Element[n + 1];

            while (c * c + d * d >= CONVTOL && iter < ITERMAX)
            {
                Element[] F = new Element[n + 1];
                B = new Element[n + 1];
                for(int i = 0; i<n+1; i++)
                {
                    F[i] = Element.zero;
                    B[i] = Element.zero;
                }
                for(int i = n-2; i >= 0; i--)
                {
                    B[i] = A[i + 2] - u * B[i + 1] - v * B[i + 2];
                    F[i] = B[i + 2] - u * F[i + 1] - v * F[i + 2];
                }
                c = A[1] - u * B[0] - v * B[1];
                d = A[0] - v * B[0];
                Element g = B[1] - u * F[0] - v * F[1];
                Element h = B[0] - v * F[0];
                Element det = v * g * g + h * (h - u * g);
                u -= (-h * c + g * d) / det;
                v -= (-g * v * c + (g * u - h) * d) / det;
                iter++;
            }

            Element[] ret = new Element[n - 1];
            for(int i = 0; i<n-1; i++)
            {
                ret[i] = B[i];
            }

            return new Tuple<Element, Element, Element[]>(u, v, ret);
        }

        private Element EvaluatePolynomial(Element[] A, Element cmplx)
        {
            Element result = Element.zero;
            int power = 0;
            foreach(var element in A)
            {
                result += element * Power(cmplx, power);
                power++;
            }
            return result;
        }

        private Tuple<Element, Element> SolveQuadraticEquation(Element a, Element b, Element c)
        {
            Element det = b * b - 4 * a * c;
            Element root1 = (-b + SquareRoot(det)) / (2.0 * a);
            Element root2 = -(b / a + root1);
            return new Tuple<Element, Element>(root1, root2);
        }

        private Element Power(Element cmplx, int power)
        {
            if(power == 0)
            {
                return Element.one;
            }
            else
            {
                Element result = cmplx;
                for(int i = 1; i<power; i++)
                {
                    result = result * cmplx;
                }
                return result;
            }
        }

        private Element SquareRoot(Element cmplx)
        {
            if(cmplx.im == 0)
            {
                return new Element(Math.Sqrt(cmplx.re));
            }
            else
            {
                double re = Math.Sqrt((cmplx.re + Math.Sqrt(cmplx.re * cmplx.re + cmplx.im * cmplx.im)) / 2);
                double im = (cmplx.im / Math.Abs(cmplx.im)) * Math.Sqrt((-cmplx.re + Math.Sqrt(cmplx.re * cmplx.re + cmplx.im * cmplx.im)) / 2);
                return new Element(re, im);
            }
        }
        //*
        //

        public Eigenvalue[] GetEigenvalues()
        {
            return eigenvalues;
        }
    }
}
