using System;
using System.Collections.Generic;

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

        public static Element zero = new Element(0, 0);
        public static Element one = new Element(1);

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
                        if((col + 1) % 2 == 0)
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
            Element[] ret = new Element[n + 1];
            if (n == 2)
            {
                if (sub[0, 0].Item2 != Element.zero && sub[1, 1].Item2 != Element.zero)
                {
                    ret[0] = sub[0, 0].Item1 * sub[1, 1].Item1 - sub[0, 1].Item1 * sub[1, 0].Item1;
                    ret[1] = -(sub[0, 0].Item1) - sub[1, 1].Item1;
                    ret[2] = Element.one;
                }
                if (sub[0, 1].Item2 != Element.zero && sub[1, 0].Item2 != Element.zero)
                {
                    ret[0] = sub[0, 0].Item1 * sub[1, 1].Item1 - sub[0, 1].Item1 * sub[1, 0].Item1;
                    ret[1] = sub[0, 1].Item1 + sub[1, 0].Item1;
                    ret[2] = -(Element.one);
                }
                else if(sub[0, 0].Item2 != Element.zero)
                {
                    ret[0] = sub[0, 0].Item1 * sub[1, 1].Item1 - sub[0, 1].Item1 * sub[1, 0].Item1;
                    ret[1] = -(sub[1, 1].Item1);
                    ret[2] = Element.zero;
                }
                else if (sub[1, 1].Item2 != Element.zero)
                {
                    ret[0] = sub[0, 0].Item1 * sub[1, 1].Item1 - sub[0, 1].Item1 * sub[1, 0].Item1;
                    ret[1] = -(sub[0, 0].Item1);
                    ret[2] = Element.zero;
                }
                else if (sub[1, 0].Item2 != Element.zero)
                {
                    ret[0] = sub[0, 0].Item1 * sub[1, 1].Item1 - sub[0, 1].Item1 * sub[1, 0].Item1;
                    ret[1] = -(sub[0, 1].Item1);
                    ret[2] = Element.zero;
                }
                else if (sub[0, 1].Item2 != Element.zero)
                {
                    ret[0] = sub[0, 0].Item1 * sub[1, 1].Item1 - sub[0, 1].Item1 * sub[1, 0].Item1;
                    ret[1] = -(sub[1, 0].Item1);
                    ret[2] = Element.zero;
                }
                else
                {
                    ret[0] = sub[0, 0].Item1 * sub[1, 1].Item1 - sub[0, 1].Item1 * sub[1, 0].Item1;
                    ret[1] = Element.zero;
                    ret[2] = Element.zero;
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

    public class Eigenvalue
    {
        protected Element eigenvalue;
        protected List<Matrix> eigenvectors;

        public Eigenvalue(Matrix m, CharacteristicPolynomial cp, Element val)
        {
            eigenvalue = val;
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
            CharacteristicPolynomial cpcalc = new CharacteristicPolynomial(m);
            Element[] cp = cpcalc.GetCharacteristicPolynomial();
            List<Element> ev = Bairstow(cp);
            eigenvalues = new Eigenvalue[ev.Count];
            for(int i = 0; i<ev.Count; i++)
            {
                eigenvalues[i] = new Eigenvalue(m, cpcalc, ev[i]);
            }
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
            const int ITERMAX = 100;
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
