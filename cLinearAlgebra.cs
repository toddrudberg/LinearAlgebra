using System;
using System.Collections.Generic;
using System.Collections;
using System.Text;
using System.Xml;
using System.Xml.Serialization;
using Electroimpact;
using Electroimpact.LinearAlgebra;

namespace Electroimpact
{

  namespace LinearAlgebra
  {

    public static class ExtensionMethods
    {
      public static double DegreesToRadians(this double Degrees)
      {
        return (Degrees * Math.PI / 180.0);
      }
      public static double RadiansToDegrees(this double Radians)
      {
        return (Radians * 180.0 / Math.PI);
      }

      /// <summary>
      /// Calculates the coefficients of a polynomial_degree powered polynomial
      /// 
      /// </summary>
      /// <param name="x">The X values</param>
      /// <param name="y">The Y values </param>
      /// <param name="polynomial_degree">the order or polynomial degree</param>
      /// <returns>an array of as long as the order, highest order coefficient first</returns>
      public static double[] LeastSquaresPolynomial(double[] x, double[] y, int polynomial_degree)
      {
        int numcoefs = polynomial_degree + 1;
        double[,] X = new double[x.Length, numcoefs];
        double[,] ys = new double[y.Length, 1];
        for (int ii = 0; ii < y.Length; ii++)
          ys[ii, 0] = y[ii];


        for (int row = 0; row < X.GetLength(0); row++)//xval
        {
          for (int col = 0; col < X.GetLength(1); col++)//power
          {
            X[row, col] = Math.Pow(x[row], col);  //put the lowest power in the first col. 0 to degree
          }
        }

        caMatrix am = new caMatrix(X);
        double[,] X_transpose = am.Transpose();

        double[,] xtx = Electroimpact.LinearAlgebra.caMatrix.DotProduct(X_transpose, X);

        Electroimpact.LinearAlgebra.cSquareMatrix M = new Electroimpact.LinearAlgebra.cSquareMatrix(xtx);
        double[,] xtx_inverse = M.Inverse();

        double[,] xtx_inerse_xt = Electroimpact.LinearAlgebra.caMatrix.DotProduct(xtx_inverse, X_transpose);

        double[,] coeffs = Electroimpact.LinearAlgebra.caMatrix.DotProduct(xtx_inerse_xt, ys);

        double[] retvals = new double[coeffs.Length];
        for (int ii = 0; ii < coeffs.GetLength(0); ii++)
          retvals[ii] = coeffs[ii,0];

        return retvals;
      }
    }

    #region Simple 6dof structure

    [Serializable]
    [XmlRoot]
    public class c6DOF_Args
    {
      public double X;
      public double Y;
      public double Z;
      /// <summary>
      /// rX is in Degrees because its readable.  This will bite me someday.
      /// </summary>
      public double rX;
      /// <summary>
      /// rY is in Degrees because its readable.  This will bite me someday.
      /// </summary>
      public double rY;
      /// <summary>
      /// rZ is in Degrees because its readable.  This will bite me someday.
      /// </summary>
      public double rZ;

      private const double d2r = Math.PI / 180.0;

      public c6DOF_Args()
      {
        X = Y = Z = rX = rY = rZ = 0.0;
      }

      public c6DOF_Args(c6DOF_Args input)
      {
        X = input.X;
        Y = input.Y;
        Z = input.Z;
        rX = input.rX;
        rY = input.rY;
        rZ = input.rZ;

      }

      public c6DOF_Args(double Xin, double Yin, double Zin, double rXrad, double rYrad, double rZrad)
      {
        X = Xin;
        Y = Yin;
        Z = Zin;
        rX = rXrad.RadiansToDegrees();
        rY = rYrad.RadiansToDegrees();
        rZ = rZrad.RadiansToDegrees();

      }

      public double rX_radians
      {
        get
        {
          return rX * d2r;
        }
        set
        {
          this.rX = value.RadiansToDegrees();
        }
      }
      public double rY_radians
      {
        get
        {
          return rY * d2r;
        }
        set
        {
          this.rY = value.RadiansToDegrees();
        }
      }
      public double rZ_radians
      {
        get
        {
          return rZ * d2r;
        }
        set
        {
          this.rZ = value.RadiansToDegrees();
        }
      }
    }
    #endregion

    #region Transform Class
    /// <summary>
    /// Transforms is a collection of 4x4 matricies.  
    /// These matricies can be formed by using DHMatrix or RussMatrix or
    /// whatever method you prefer.
    /// 
    /// Calling DotProduct(), returns the dot product of all the matricies in this system
    /// in the order they were added to the collection.
    /// 
    /// This Class also has a CrossProduct(double[] a, double[] b) method.  
    /// This returns a x b and is unrelated to the collection of matricies.
    /// </summary>
    public class Transforms
    {
      #region Members
      private System.Collections.ArrayList myTransforms = new ArrayList();
      private cMatrix myResult = new cMatrix();
      #endregion

      #region Constructor
      /// <summary>
      /// Default constructor
      /// </summary>
      public Transforms()
      {
      }
      #endregion

      #region Methods
      public void Clear()
      {
        this.myTransforms.Clear();
      }
      /// <summary>
      /// Cross Product is unrelated to the collection of matricies in this 
      /// class and merely returns the i, j, k scalers of the a x b function.
      /// </summary>
      /// <param name="a">a one dimensional double array of three elements.
      /// a[0] = i scaler
      /// a[1] = j scaler
      /// a[2] = k scaler</param>
      /// <param name="b">similar to a.</param>
      /// <returns>a one dimensional double array of three elements.
      /// [0] = i scaler
      /// [1] = j scaler
      /// [2] = k scaler.
      /// This function does not normalize the vector</returns>
      public double[] CrossProduct(double[] a, double[] b)
      {
        double vi, vj, vk;

        if (a.Length < 3 || b.Length < 3)
        {
          System.Exception ex = new Exception("a and b must be 3 element doubles");
          throw (ex);
        }

        vi = a[1] * b[2] - a[2] * b[1];
        vj = -(a[0] * b[2] - a[2] * b[0]);
        vk = a[0] * b[1] - a[1] * b[0];

        double[] dog = new double[3];

        dog[0] = vi;
        dog[1] = vj;
        dog[2] = vk;
        return dog;
      }
      /// <summary>
      /// Normalizes a one dimensional matrix a
      /// </summary>
      /// <param name="a">a one dimensional double matrix</param>
      /// <returns>a one dimensional double matrix with same number of elements as a</returns>
      public double[] Normalize(double[] a)
      {
        double hold = 0;
        double[] rtn = new double[a.Length];
        for (int ii = 0; ii < a.Length; ii++)
        {
          hold += a[ii] * a[ii];
        }
        hold = Math.Sqrt(hold);
        for (int ii = 0; ii < a.Length; ii++)
        {
          rtn[ii] = a[ii] / hold;
        }
        return rtn;
      }
      /// <summary>
      /// Calculates the dot product of the collection of matricies in this 
      /// class and returns the 4x4 result.
      /// </summary>
      /// <returns>Returns a double[4,4] matrix</returns>
      public double[,] DotProduct(uint start, uint finish)
      {
        double[,] Hold = new double[4, 4];
        double[,] Result = new double[4, 4];
        double[,] Transform = new double[4, 4];

        if (start > this.myTransforms.Count - 1 || finish > this.myTransforms.Count - 1 || start >= finish)
        {
          System.Exception ex = new Exception("screw up in DotProduct.");
        }

        Copy4by4(ref Result, (double[,])this.myTransforms[(int)start]);

        for (uint ii = start + 1; ii <= finish; ii++)
        {
          Copy4by4(ref Transform, (double[,])this.myTransforms[(int)ii]);
          //for (int row = 0; row < 4; row++)
          //{
          //  for (int col = 0; col < 4; col++)
          //  {
          //    Hold[row, col] =
          //      Result[row, 0] * Transform[0, col] +
          //      Result[row, 1] * Transform[1, col] +
          //      Result[row, 2] * Transform[2, col] +
          //      Result[row, 3] * Transform[3, col];
          //  }
          //}
          Hold = DotProductLHT(Result, Transform);
          Copy4by4(ref Result, Hold);
        }
        myResult = new cMatrix(Result);
        return Result;
      }
      /// <summary>
      /// Computes dotproduct of member 4x4 matrices
      /// </summary>
      /// <returns>double[4,4]</returns>
      public double[,] DotProduct()
      {
        double[,] Hold = new double[4, 4];
        double[,] Result = new double[4, 4];
        double[,] Transform = new double[4, 4];

        if (this.myTransforms.Count > 0)
        {
          //Result = (double[,])this.myTransforms[0];
          Copy4by4(ref Result, (double[,])this.myTransforms[0]);
        }
        else
        {
          System.Exception ex = new Exception("There are no Transforms to multiply.");
          throw ex;
        }
        for (int ii = 1; ii < this.myTransforms.Count; ii++)
        {
          //Transform = (double[,])this.myTransforms[ii];
          Copy4by4(ref Transform, (double[,])this.myTransforms[ii]);
          //Hold = DotProduct(Result, Transform);
          Hold = DotProductLHT(Result, Transform);
          Copy4by4(ref Result, Hold);
        }
        myResult = new cMatrix(Result);
        return Result;
      }
      /// <summary>
      /// use for any 4x4 transform
      /// </summary>
      /// <param name="left"></param>
      /// <param name="right"></param>
      /// <returns></returns>
      Double[,] DotProduct(double[,] left, double[,] right)
      {
        double[,] Result = new double[4, 4];
        for (int row = 0; row < 4; row++)
        {
          for (int col = 0; col < 4; col++)
          {
            Result[row, col] =
              left[row, 0] * right[0, col] +
              left[row, 1] * right[1, col] +
              left[row, 2] * right[2, col] +
              left[row, 3] * right[3, col];
          }
        }
        return Result;
      }
      /// <summary>
      /// Use only with LHT transforms
      /// </summary>
      /// <param name="left"></param>
      /// <param name="right"></param>
      /// <returns></returns>
      public double[,] DotProductLHT(double[,] left, double[,] right)
      {
        double[,] Result = new double[4,4];
        //1st row: (3 mult 2 add) * 4 + 1 add = 12 mult, 9 add total
        Result[0,0] = left[0,0] * right[0,0] + left[0,1] * right[1,0] + left[0,2] * right[2,0]; //right[3,0] is always 0
        Result[0,1] = left[0,0] * right[0,1] + left[0,1] * right[1,1] + left[0,2] * right[2,1]; //right[3,1] is always 0
        Result[0,2] = left[0,0] * right[0,2] + left[0,1] * right[1,2] + left[0,2] * right[2,2]; //right[3,2] is always 0
        Result[0,3] = left[0,0] * right[0,3] + left[0,1] * right[1,3] + left[0,2] * right[2,3] + left[0,3];  //right[3,3] is always 1

        //2nd row:  (3 mult 2 add) * 4 + 1 add = 12 mult, 9 add total
        Result[1,0] = left[1,0] * right[0,0] + left[1,1] * right[1,0] + left[1,2] * right[2,0]; //right[3,0] is always 0
        Result[1,1] = left[1,0] * right[0,1] + left[1,1] * right[1,1] + left[1,2] * right[2,1]; //right[3,1] is always 0
        Result[1,2] = left[1,0] * right[0,2] + left[1,1] * right[1,2] + left[1,2] * right[2,2]; //right[3,2] is always 0
        Result[1,3] = left[1,0] * right[0,3] + left[1,1] * right[1,3] + left[1,2] * right[2,3] + left[1,3];             //right[3,3] is always 1


        //3rd row: (3 mult 2 add) * 4 + 1 add = 12 mult, 9 add total
        Result[2,0] = left[2,0] * right[0,0] + left[2,1] * right[1,0] + left[2,2] * right[2,0]; //right[3,0] is always 0
        Result[2,1] = left[2,0] * right[0,1] + left[2,1] * right[1,1] + left[2,2] * right[2,1]; //right[3,1] is always 0
        Result[2,2] = left[2,0] * right[0,2] + left[2,1] * right[1,2] + left[2,2] * right[2,2]; //right[3,2] is always 0
        Result[2,3] = left[2,0] * right[0,3] + left[2,1] * right[1,3] + left[2,2] * right[2,3] + left[2,3];             //right[3,3] is always 1

        //4th row: no math
        Result[3,0] = 0.0;
        Result[3,1] = 0.0;
        Result[3,2] = 0.0;
        Result[3,3] = 1.0;

        //total cost is 36 mult, 27 add

        return Result;
      }


      /// <summary>
      /// Adds a 4x4 matrix to the collection of transforms in this class.
      /// </summary>
      /// <param name="iXform">a 4x4 double matrix</param>
      public void AddTransform(double[,] iXform)
      {
        if (iXform.GetLength(0) == 4 && iXform.GetLength(1) == 4)
        {
          this.myTransforms.Add(iXform);
        }
        else
        {
          System.Exception ex = new Exception("must be a 4x4 matrix");
          throw (ex);
        }
      }

      public cMatrix myStateMatrix
      {
        get
        {
          return myResult;
        }
      }

      /// <summary>
      /// Copies the 4x4 matrix iXform into the location indicated by Element
      /// which is an index into the collection of transforms in this class.
      /// </summary>
      /// <param name="Element">Integer element to modify.</param>
      /// <param name="iXform">double[4,4] to change it to.</param>
      public void ModifyTransform(int Element, double[,] iXform)
      {
        if (iXform.GetLength(0) == 4 && iXform.GetLength(1) == 4)
        {
          if (Element < this.myTransforms.Count)
            this.myTransforms[Element] = iXform;
        }
      }

      private void Copy4by4(ref double[,] left, double[,] right)
      {
        for (int row = 0; row < 4; row++)
        {
          for (int col = 0; col < 4; col++)
          {
            left[row, col] = right[row, col];
          }
        }
      }

      #endregion

      #region Properties
      /// <summary>
      /// Returns the number of Transforms in this system.
      /// </summary>
      public int Count
      {
        get { return this.myTransforms.Count; }
      }

      public double X
      {
        get { return this.myResult.X; }
      }
      public double Y
      {
        get { return this.myResult.Y; }
      }
      public double Z
      {
        get { return this.myResult.Z; }
      }
      public double A
      {
        get { return this.myResult.A; }
      }

      public double[] rXrYrZFixedAngle
      {
        get { return this.myResult.rXrYrZFixedAngle; }
      }

      public double[,] rZrYrX
      {
        get { return this.myResult.rZrYrX; }
      }

      public double[,] rXrYrZ
      {
        get { return this.myResult.rXrYrZ; }
      }

      public double[,] rYrXrZ
      {
        get { return this.myResult.rYrXrZ; }
      }

      public double B
      {
        get { return this.myResult.B; }
      }

      public double[] C
      {
        get { return this.myResult.C; }
      }

      #endregion
    }
    #endregion

    #region a Matrix

    public class caMatrix
    {
      private double[,] myM;

      public caMatrix(double[,] input)
      {
        myM = input;
      }

      //public double[,] Transpose()
      //{
        
      //}

      public double[,] Transpose()
      {
        int nrows = myM.GetLength(0);
        int ncols = myM.GetLength(1);

        double[,] t = new double[ncols, nrows];

        for (int row = 0; row < t.GetLength(0); row++)
        {
          for (int col = 0; col < t.GetLength(1); col++)
          {
            t[row, col] = myM[col, row];
          }
        }
        return t;
      }

      public static double[,] DotProduct(double[,] left, double[,] right)
      {

        int rows = left.GetLength(0);
        int cols = right.GetLength(1);

        double[,] Result = new double[rows, cols];

        int left_cols = left.GetLength(1);
        int right_rows = right.GetLength(0);

        if (left_cols != right_rows)
        {
          System.Windows.Forms.MessageBox.Show("dimensions not compatible");
        }

        for (int row = 0; row < rows; row++)
        {
          for (int col = 0; col < cols; col++)
          {

            double result = 0;
            for (int elem = 0; elem < left_cols; elem++)
            {
              result += left[row, elem] * right[elem, col];
            }

            Result[row, col] = result;
          }
        }
        return Result;
      }
    }



    #endregion


    #region nxn Matrix Class
    public class cSquareMatrix
    {
      private double[,] myM;
      private int size;
      private const double d2r = Math.PI / 180.0;
      private const double r2d = 180.0 / Math.PI;

      public cSquareMatrix(double[,] Matrix)
      {
        size = Matrix.GetLength(0);
        myM = new double[size, size];
        myM = CopyMatrix(Matrix);
      }

      public cSquareMatrix(double[,] Matrix, int Size)
      {
        size = Size;
        myM = new double[size, size];
        myM = CopyMatrix(Matrix);
      }

      /// <summary>
      /// Dot Product of this matrix with single dimensional n element array of doubles.
      /// </summary>
      /// <param name="Rvector">single dimensional matrix that maches size of this matrix class of doubles</param>
      /// <returns>Single dimensional n element array of doubles</returns>
      /// 
      public double[] DotMe(double[] Rvector)
      {
        double[,] Hold = new double[size, size];
        double[,] Result = new double[size, size];
        double[,] Transform = new double[size, size];
        double[] dog = new double[size];

        for (int ii = 0; ii < size; ii++)
        {
          dog[ii] = 0;
          for (int jj = 0; jj < size; jj++)
            dog[ii] += myM[ii, jj] * Rvector[jj];
        }
        return dog;
      }

      public double Determinate()
      {
        double det = 1;
        double[,] copy = CopyMatrix();

        for (int diag = 0; diag < size; diag++)
        {
          double element = copy[diag, diag];
          if (element == 0)
          {
            int jj = diag;
            while (jj < size && copy[jj, jj] == 0)
            {
              jj++;
            }
            if (jj == size)
              return 0.0;
            for (int ii = diag; ii < size; ii++)
            {
              double hold = copy[ii, jj];
              copy[ii, jj] = copy[ii, diag];
              copy[ii, diag] = hold;
            }
            det = -det;
          }
          else
          {
            det = det * element;
            if (diag < size - 1)
            {
              int jj = diag + 1;
              for (int row = jj; row < size; row++)
              {
                for (int col = jj; col < size; col++)
                  copy[row, col] = copy[row, col] - copy[row, diag] * (copy[diag, col] / element);
              }
            }
          }
        }
        return det;
      }

      public double[,] Inverse2x2()
      {
        double det = Determinate();
        double[,] m = new double[2, 2];

        m[0, 0] = myM[1, 1] / det;
        m[0, 1] = -myM[0, 1] / det;
        m[1, 0] = -myM[1, 0] / det;
        m[1, 1] = myM[0, 0] / det;
        return m;
      }

      public double[,] Inverse()
      {
        double[,] copy = CopyMatrix();
        double[,] id = IdentityMatrix();
        double[,] b = new double[size, size];

        for (int diag = 0; diag < size; diag++)
        {
          if (Math.Abs(copy[diag, diag]) < 0.0000000001)
          {
            for (int col = diag + 1; col < size; col++)
            {
              if (col != diag)
              {
                if (Math.Abs(copy[col, diag]) > 0.0000000001)
                {
                  for (int row = 0; row < size; row++)
                  {
                    copy[row, diag] = copy[row, diag] + copy[row, col];
                    id[row, diag] = id[row, diag] + id[row, col];
                  }
                  break;
                }
              }
              else
              {
                System.Exception ex = new Exception("Mother Fucker!");
                throw ex;
              }
            }
          }
          double N = 1.0 / copy[diag, diag];
          for (int row = 0; row < size; row++) //normalize the diagonal
          {
            copy[row, diag] = N * copy[row, diag];
            id[row, diag] = N * id[row, diag];
          }
          for (int col = 0; col < size; col++)
          {
            if (col != diag)
            {
              double a = copy[diag, col];
              for (int row = 0; row < size; row++)
              {
                copy[row, col] = copy[row, col] - a * copy[row, diag];
                id[row, col] = id[row, col] - a * id[row, diag];
              }
            }
          }
        }
        return id;
      }

      private double[,] CopyMatrix()
      {
        return CopyMatrix(myM);
      }

      private double[,] CopyMatrix(double[,] Matrix)
      {
        double[,] copy = new double[size, size];
        for (int row = 0; row < size; row++)
        {
          for (int col = 0; col < size; col++)
            copy[row, col] = Matrix[row, col];
        }
        return copy;
      }

      private double[,] IdentityMatrix()
      {
        return IdentityMatrix(size);
      }
      private double[,] IdentityMatrix(int Size)
      {
        double[,] r = new double[Size, Size];
        for (int row = 0; row < Size; row++)
        {
          for (int col = 0; col < Size; col++)
          {
            r[row, col] = 0;
          }
          r[row, row] = 1;
        }
        return r;
      }
    }
    #endregion

    #region Matrix Manipulator Class
    /// <summary>
    /// cMatrix is a 4x4 LHT matrix manipulator.
    /// </summary>
    public class cMatrix
    {
      int size = 4;
      private double[,] myM = new double[4, 4];
      private const double d2r = Math.PI / 180.0;
      private const double SmallDouble = 1.0e-15;
      /// <summary>
      /// Rutuns the angle of the vector composed of x and y as measured from the positive x axis going counter clockwise.
      /// </summary>
      /// <param name="x">component value of the vector from the axis which the angle is to be measured</param>
      /// <param name="y">component value fo the vector from the axis orthoganal to the above axis</param>
      /// <returns>double -- the angle.</returns>
      public double aTan2(double x, double y)
      {
        double r = Math.Sqrt(x * x + y * y);

        //if x is very small relative to r we are on or near an asymtote of the tan function.
        if (Math.Abs(x) < .0000001 * r)
        {
          if (y >= 0)
            return 90.0;
          else
            return 270.0;

        }
        double a = Math.Atan(y / x).RadiansToDegrees();
        if (x > 0.0)
        {
          a = a < 0 ? a + 360 : a;
          return a;
        }
        else
        {
          return 180 + a;
        }
      }

      /// <summary>
      /// Dot Product of this matrix with single dimensional 4 element array of doubles.
      /// </summary>
      /// <param name="Rvector">single dimensional 4 element array of doubles</param>
      /// <returns>Single dimensional 4 element array of doubles</returns>
      /// 
      public double[] DotMe(double[] Rin)
      {
        double[,] Hold = new double[4, 4];
        double[,] Result = new double[4, 4];
        double[,] Transform = new double[4, 4];
        double[] dog = new double[4];
        double[] Rvector;
        if (Rin.Length == 3)
          Rvector = new double[] { Rin[0], Rin[1], Rin[2], 1.0 };
        else
          Rvector = new double[] { Rin[0], Rin[1], Rin[2], Rin[3] };

        for (int ii = 0; ii < 4; ii++)
        {
          dog[ii] = myM[ii, 0] * Rvector[0] +
                    myM[ii, 1] * Rvector[1] +
                    myM[ii, 2] * Rvector[2] +
                    myM[ii, 3] * Rvector[3];
        }
        return dog;
      }

      public double[] DotMe(Vector R)
      {
        return DotMe(R.R);

      }
      private void Copy4by4(ref double[,] left, double[,] right)
      {
        for (int row = 0; row < 4; row++)
        {
          for (int col = 0; col < 4; col++)
          {
            left[row, col] = right[row, col];
          }
        }
      }

      /// <summary>
      /// Multiplies the matrix class matrix x MatrixIn and returns the result (the dot product).
      /// </summary>
      /// <param name="MatrixIn">a 4x4 double, there is no error checking!</param>
      /// <returns>a 4x4 double</returns>
      public double[,] DotMe(double[,] MatrixIn)
      {
        double[,] Result = new double[4, 4];

        for (int row = 0; row < 4; row++)
        {
          for (int col = 0; col < 4; col++)
          {
            Result[row, col] =
              myM[row, 0] * MatrixIn[0, col] +
              myM[row, 1] * MatrixIn[1, col] +
              myM[row, 2] * MatrixIn[2, col] +
              myM[row, 3] * MatrixIn[3, col];
          }
        }
        return Result;
      }

      /// <summary>
      /// Warning! assumed LHT matrix
      /// Multiplies the matrix class matrix x MatrixIn and returns the result (the dot product).
      /// </summary>
      /// <param name="MatrixIn">a 4x4 double, there is no error checking!</param>
      /// <returns>a 4x4 double</returns>
      public double[,] DotMeLHT(double[,] right)
      {
        {
          double[,] Result = new double[4,4];
          //1st row: (3 mult 2 add) * 4 + 1 add = 12 mult, 9 add total
          Result[0,0] = myM[0,0] * right[0,0] + myM[0,1] * right[1,0] + myM[0,2] * right[2,0]; //right[3,0] is always 0
          Result[0,1] = myM[0,0] * right[0,1] + myM[0,1] * right[1,1] + myM[0,2] * right[2,1]; //right[3,1] is always 0
          Result[0,2] = myM[0,0] * right[0,2] + myM[0,1] * right[1,2] + myM[0,2] * right[2,2]; //right[3,2] is always 0
          Result[0,3] = myM[0,0] * right[0,3] + myM[0,1] * right[1,3] + myM[0,2] * right[2,3] + myM[0,3];  //right[3,3] is always 1

          //2nd row:  (3 mult 2 add) * 4 + 1 add = 12 mult, 9 add total
          Result[1,0] = myM[1,0] * right[0,0] + myM[1,1] * right[1,0] + myM[1,2] * right[2,0]; //right[3,0] is always 0
          Result[1,1] = myM[1,0] * right[0,1] + myM[1,1] * right[1,1] + myM[1,2] * right[2,1]; //right[3,1] is always 0
          Result[1,2] = myM[1,0] * right[0,2] + myM[1,1] * right[1,2] + myM[1,2] * right[2,2]; //right[3,2] is always 0
          Result[1,3] = myM[1,0] * right[0,3] + myM[1,1] * right[1,3] + myM[1,2] * right[2,3] + myM[1,3];             //right[3,3] is always 1


          //3rd row: (3 mult 2 add) * 4 + 1 add = 12 mult, 9 add total
          Result[2,0] = myM[2,0] * right[0,0] + myM[2,1] * right[1,0] + myM[2,2] * right[2,0]; //right[3,0] is always 0
          Result[2,1] = myM[2,0] * right[0,1] + myM[2,1] * right[1,1] + myM[2,2] * right[2,1]; //right[3,1] is always 0
          Result[2,2] = myM[2,0] * right[0,2] + myM[2,1] * right[1,2] + myM[2,2] * right[2,2]; //right[3,2] is always 0
          Result[2,3] = myM[2,0] * right[0,3] + myM[2,1] * right[1,3] + myM[2,2] * right[2,3] + myM[2,3];             //right[3,3] is always 1

          //4th row: no math
          Result[3,0] = 0.0;
          Result[3,1] = 0.0;
          Result[3,2] = 0.0;
          Result[3,3] = 1.0;

          //total cost is 36 mult, 27 add

          return Result;
        }

        //double[,] Result = new double[4, 4];

        //for (int row = 0; row < 4; row++)
        //{
        //  for (int col = 0; col < 4; col++)
        //  {
        //    Result[row, col] =
        //      myM[row, 0] * MatrixIn[0, col] +
        //      myM[row, 1] * MatrixIn[1, col] +
        //      myM[row, 2] * MatrixIn[2, col] +
        //      myM[row, 3] * MatrixIn[3, col];
        //  }
        //}
        //return Result;
      }

      /// <summary>
      /// Constructor requires a 4x4 matrix for input.
      /// </summary>
      /// <param name="M">a 4x4 double array</param>
      public cMatrix(double[,] M)
      {
        if (M.GetLength(0) != 4 || M.GetLength(1) != 4)
          throw (new System.Exception("Problem in cMatrix Constructor."));
        this.myM = M;
      }
      /// <summary>
      /// Matrix will be the identity matrix.
      /// </summary>
      public cMatrix()
      {
        for (int row = 0; row < 4; row++)
        {
          for (int col = 0; col < 4; col++)
          {
            if (row == col)
              this.myM[row, col] = 1;
            else
              this.myM[row, col] = 0;
          }
        }
      }
      /// <summary>
      /// Only works for standard transfer matrix.
      /// 
      /// </summary>

      public void InvertMe()
      {
        myM = Inverse();
      }
      private double[,] CopyMatrix()
      {
        return CopyMatrix(myM);
      }
      private double[,] CopyMatrix(double[,] Matrix)
      {
        double[,] copy = new double[size, size];
        for (int row = 0; row < size; row++)
        {
          for (int col = 0; col < size; col++)
            copy[row, col] = Matrix[row, col];
        }
        return copy;
      }
      private double[,] IdentityMatrix()
      {
        return IdentityMatrix(size);
      }
      private double[,] IdentityMatrix(int Size)
      {
        double[,] r = new double[Size, Size];
        for (int row = 0; row < Size; row++)
        {
          for (int col = 0; col < Size; col++)
          {
            r[row, col] = 0;
          }
          r[row, row] = 1;
        }
        return r;
      }
      /// <summary>
      /// Returns the 4x4 of the inverse matrix.  This works only for homegeneous matrices.
      /// </summary>
      /// <returns></returns>
      public double[,] InverseHomogeneous()
      {
        cMatrix Mprime = new cMatrix(this.TransposeOrientation());
        double[] Rprime = this.R;
        Rprime = Mprime.DotMe(Rprime);
        Rprime[0] *= -1.0;
        Rprime[1] *= -1.0;
        Rprime[2] *= -1.0;
        Mprime.R = Rprime;
        return Mprime.GetMatrix;
      }

      public double[,] Inverse()
      {
        double[,] copy = CopyMatrix();
        double[,] id = IdentityMatrix();
        double[,] b = new double[size, size];

        for (int diag = 0; diag < size; diag++)
        {
          if (Math.Abs(copy[diag, diag]) < 0.0000000001)
          {
            for (int col = diag + 1; col < size; col++)
            {
              if (col != diag)
              {
                if (Math.Abs(copy[col, diag]) > 0.0000000001)
                {
                  for (int row = 0; row < size; row++)
                  {
                    copy[row, diag] = copy[row, diag] + copy[row, col];
                    id[row, diag] = id[row, diag] + id[row, col];
                  }
                  break;
                }
              }
              else
              {
                System.Exception ex = new Exception("Mother Fucker!");
                throw ex;
              }
            }
          }
          double N = 1.0 / copy[diag, diag];
          for (int row = 0; row < size; row++) //normalize the diagonal
          {
            copy[row, diag] = N * copy[row, diag];
            id[row, diag] = N * id[row, diag];
          }
          for (int col = 0; col < size; col++)
          {
            if (col != diag)
            {
              double a = copy[diag, col];
              for (int row = 0; row < size; row++)
              {
                copy[row, col] = copy[row, col] - a * copy[row, diag];
                id[row, col] = id[row, col] - a * id[row, diag];
              }
            }
          }
        }
        return id;
      }

      public double[,] TransposeOrientation()
      {
        double[,] m = new double[4, 4];
        m[0, 0] = vi[0];
        m[0, 1] = vi[1];
        m[0, 2] = vi[2];
        m[1, 0] = vj[0];
        m[1, 1] = vj[1];
        m[1, 2] = vj[2];
        m[2, 0] = vk[0];
        m[2, 1] = vk[1];
        m[2, 2] = vk[2];

        //R is no longer probably interesting.
        m[0, 3] = 0;
        m[1, 3] = 0;
        m[2, 3] = 0;

        //The bottom row should be 0001
        m[3, 0] = 0;
        m[3, 1] = 0;
        m[3, 2] = 0;
        m[3, 3] = 1;
        return m;
      }

      /// <summary>
      /// get the full transpose of this 4x4.  Solid works format is the transpose of the rest of the known universe's for instance.
      /// </summary>
      /// <returns></returns>
      public double[,] GetTranspose()
      {
        double[,] m = new double[4, 4];
        m[0, 0] = vi[0];
        m[0, 1] = vi[1];
        m[0, 2] = vi[2];
        m[1, 0] = vj[0];
        m[1, 1] = vj[1];
        m[1, 2] = vj[2];
        m[2, 0] = vk[0];
        m[2, 1] = vk[1];
        m[2, 2] = vk[2];

        m[3, 0] = vR[0];
        m[3, 1] = vR[1];
        m[3, 2] = vR[2];

        m[0, 3] = 0;
        m[1, 3] = 0;
        m[2, 3] = 0;
        m[3, 3] = 1;
        return m;
      }

      public double VectorDot(double[] left, double[] right)
      {
        double ret = 0;
        ret = left[0] * right[0] + left[1] * right[1] + left[2] * right[2];
        return ret;
      }
      public double[,] GetMatrix2()
      {
        return this.GetMatrix;
      }
      public double[,] GetMatrix
      {
        get { return myM; }
      }
      /// <summary>
      /// Returns string representing 4x4 array.  Units are unitless and mm.
      /// </summary>
      /// <returns></returns>
      public override string ToString()
      {
        return ProduceString(false);
      }

      /// <summary>
      /// Returns string representing 4x4 array.  Units are unitless and mm.
      /// </summary>
      /// <returns></returns>
      public string ToString(bool ConvertToInch)
      {
        return ProduceString(true);
      }

      private string ProduceString(bool ConvertToInch)
      {
        string sz = "";
        for (int row = 0; row < 4; row++)
        {
          //sz += "[ ";
          for (int col = 0; col < 3; col++)
          {
            sz += " " + this.myM[row, col].ToString("f6") + "\t";
          }
          if (ConvertToInch && row < 3)
            sz += " " + (this.myM[row, 3] / 25.4).ToString("f6") + "\t";
          else
            sz += " " + this.myM[row, 3].ToString("f6") + "\t";
          sz += "\n";
        }
        return sz;
      }
      /// <summary>
      /// Returns the 0,3 element of array or the tp.x
      /// </summary>
      public double X
      {
        get { return this.myM[0, 3]; }
      }
      /// <summary>
      /// return 1,3 element of array or the tp.y
      /// </summary>
      public double Y
      {
        get { return this.myM[1, 3]; }
      }
      /// <summary>
      /// returns the 2,3 element of array or the tp.z
      /// </summary>
      public double Z
      {
        get { return this.myM[2, 3]; }
      }
      /// <summary>
      /// return one dimensional 3 element array representing i.
      /// </summary>
      public double[] i
      {
        get
        {
          double[] dog = new double[3];
          dog[0] = this.myM[0, 0];
          dog[1] = this.myM[1, 0];
          dog[2] = this.myM[2, 0];
          dog = this.Normalize(dog);
          return dog;
        }
        set
        {
          if (value.Length > 2)
          {
            myM[0, 0] = value[0];
            myM[1, 0] = value[1];
            myM[2, 0] = value[2];
            myM[3, 0] = 0;
          }
        }
      }
      public Vector vi
      {
        get
        {
          return new Vector(this.i);
        }
      }
      /// <summary>
      /// return one dimensional 3 element array representing j.
      /// </summary>
      public double[] j
      {
        get
        {
          double[] dog = new double[3];
          dog[0] = this.myM[0, 1];
          dog[1] = this.myM[1, 1];
          dog[2] = this.myM[2, 1];
          dog = this.Normalize(dog);
          return dog;
        }
        set
        {
          if (value.Length > 2)
          {
            myM[0, 1] = value[0];
            myM[1, 1] = value[1];
            myM[2, 1] = value[2];
            myM[3, 1] = 0;
          }
        }
      }
      public Vector vj
      {
        get
        {
          return new Vector(this.j);
        }
      }
      /// <summary>
      /// return one dimensional 3 element array representing k.
      /// </summary>
      public double[] k
      {
        get
        {
          double[] dog = new double[3];
          dog[0] = this.myM[0, 2];
          dog[1] = this.myM[1, 2];
          dog[2] = this.myM[2, 2];
          dog = this.Normalize(dog);
          return dog;
        }
        set
        {
          if (value.Length > 2)
          {
            myM[0, 2] = value[0];
            myM[1, 2] = value[1];
            myM[2, 2] = value[2];
            myM[3, 2] = 0;
          }
        }
      }
      public Vector vk
      {
        get
        {
          return new Vector(this.k);
        }
      }

      /// <summary>
      /// Returns a two element double array, where: <br/>
      /// A[0] = A in terms of 360 axis <br/>
      /// A[1] = A in terms of a +/- axis.
      /// Calculates the Atp Based on the k vector.  Units are degrees.
      /// </summary>
      public double A
      {
        get
        {
          return -Math.Atan2(k[1], k[2]).RadiansToDegrees();
        }
      }

      /// <summary>
      /// rXrYrZFixedAngle returns a double array. In DEGREES.
      /// This is NOT the same as the old school ABC.
      /// </summary>
      public double[] rXrYrZFixedAngle //rotation about X, Y and then Z
      {
        get
        {
          double[] ret = new double[3];

          double rY = Math.Atan2(-i[2], Math.Sqrt(j[2] * j[2] + k[2] * k[2]));
          double rZ = Math.Atan2(i[1], i[0]);
          double rX = Math.Atan2(j[2], k[2]);

          //double rY = Math.Atan2(-i[2], Math.Sqrt(i[1] * i[1] + i[2] * i[2]));
          //double rZ = Math.Atan2(j[2], k[2]);
          //double rX = Math.Atan2(i[0], i[1]);


          ret[0] = rX.RadiansToDegrees();
          ret[1] = rY.RadiansToDegrees();
          ret[2] = rZ.RadiansToDegrees();





          return ret;
        }
      }


      /// <summary>
      /// rXrYrz returns a 3 x 2 matrix.  Column one is one solution, Column two is the other.
      /// Utilize your axis limits to determine the correct solution for this case..
      /// Returns two solutions, both are in degrees.  
      /// [0,0] is rX solution 1 and [0,1] is rX solution 2
      /// </summary>
      /// 

      public double[,] rXrYrZ //rotation about X, Y and then Z
      {
        get
        {
          double[,] ret = new double[3, 2];

          //First Solution
          double rY = Math.Atan2(k[0], Math.Sqrt(k[1] * k[1] + k[2] * k[2]));
          double rX = Math.Atan2(-k[1], k[2]);
          double rZ = Math.Atan2(-j[0], i[0]);

          ret[0, 0] = rX.RadiansToDegrees();
          ret[1, 0] = rY.RadiansToDegrees();
          ret[2, 0] = rZ.RadiansToDegrees();

          //Second Solution
          rY = Math.Atan2(k[0], -Math.Sqrt(k[1] * k[1] + k[2] * k[2]));
          double cosb = Math.Sign(Math.Cos(rY));

          rX = Math.Atan2(-k[1] / cosb, k[2] / cosb);
          rZ = Math.Atan2(-j[0] / cosb, i[0] / cosb);

          ret[0, 1] = rX.RadiansToDegrees();
          ret[1, 1] = rY.RadiansToDegrees();
          ret[2, 1] = rZ.RadiansToDegrees();

          return ret;
        }
      }

      public double[,] rXrYrZtest //rotation about X, Y and then Z
      {
        get
        {
          double[,] ret = new double[3, 2];

          double rY = JoshAtan2(k[0], Math.Sqrt(k[1] * k[1] + k[2] * k[2]));
          rY = Math.Asin(k[0]);
          double cosb = Math.Cos(rY);
          double rX = 0.0;
          double rZ = 0.0;
          double rzp = 0.0;

          double smallNumber = .00001;

          if (Math.Abs(cosb) < smallNumber)
            cosb = Math.Sign(cosb) * smallNumber;
          rX = JoshAtan2(-k[1] / cosb, k[2] / cosb);
          rX = Math.Asin(-k[1] / cosb);

          double a = Math.Cos(rX);
          double b = Math.Sin(rX) * Math.Sin(rY);
          double c = Math.Sin(rX);
          double d = Math.Cos(rX) * Math.Sin(rY);
          double e = Math.Cos(rX);
          double f = Math.Sin(rX) * Math.Sin(rY);
          double g = Math.Sin(rX);
          double h = Math.Cos(rX) * Math.Sin(rY);

          double cosrz = (h * j[1] + f * j[2]) / (h * e + f * g);
          double sinrz = (i[2] * b + i[1] * d) / (d * a + b * c);
          rzp = Math.Asin((i[2] * b + i[1] * d) / (d * a + b * c));
          rZ = JoshAtan2(-j[0] / cosb, i[0] / cosb);
          //rZ = Math.Atan2(-j[0] / cosb, i[0] / cosb);
          //rZ = Math.Asin(-j[0] / cosb);
          Console.WriteLine(sinrz);
          Console.WriteLine(cosrz);
          Console.WriteLine(rZ * 180.0 / Math.PI);
          Console.WriteLine(rzp * 180.0 / Math.PI);

          ret[0, 0] = rX.RadiansToDegrees();
          ret[1, 0] = rY.RadiansToDegrees();
          ret[2, 0] = rZ.RadiansToDegrees();

          rY = JoshAtan2(k[0], -Math.Sqrt(k[1] * k[1] + k[2] * k[2]));
          rY = Math.Asin(k[0]);
          cosb = Math.Cos(rY);

          if (Math.Abs(cosb) < smallNumber)
            cosb = Math.Sign(cosb) * smallNumber;
          rX = JoshAtan2(-k[1] / cosb, k[2] / cosb);
          rX = Math.Asin(-k[1] / cosb);
          rZ = JoshAtan2(-j[0] / cosb, i[0] / cosb);
          //rZ = Math.Asin(-j[0] / cosb);
          ret[0, 1] = rX.RadiansToDegrees();
          ret[1, 1] = rY.RadiansToDegrees();
          ret[2, 1] = rZ.RadiansToDegrees();

          return ret;
        }
      }

      public double JoshAtan2(double Y, double X)
      {
        double alpha = Math.Atan(Y / X);
        double answer = alpha;
        if (X < 0.0)
          answer = Math.PI + alpha;
        else if (Y < 0.0)
          answer = 2 * Math.PI - alpha;
        return answer;
      }
      /// <summary>
      /// returns values in degrees.  both solutions, rx is[0,0] and [0,1] etc
      /// </summary>
      public double[,] rYrXrZ
      {
        get
        {
          double[,] ret = new double[3, 2];

          double rX = Math.Atan2(-k[1], Math.Sqrt(i[1] * i[1] + j[1] * j[1]));
          double cosa = Math.Cos(rX);
          double rY = 0.0;
          double rZ = 0.0;

          double smallNumber = .00001;

          if (Math.Abs(cosa) < smallNumber)
            cosa = Math.Sign(cosa) * smallNumber;

          rY = Math.Atan2(k[0] / cosa, k[2] / cosa);
          rZ = Math.Atan2(i[1] / cosa, j[1] / cosa);
          ret[0, 0] = rX.RadiansToDegrees();
          ret[1, 0] = rY.RadiansToDegrees();
          ret[2, 0] = rZ.RadiansToDegrees();

          rX = Math.Atan2(-k[1], -Math.Sqrt(i[1] * i[1] + j[1] * j[1]));
          cosa = Math.Cos(rX);
          if (Math.Abs(cosa) < smallNumber)
            cosa = Math.Sign(cosa) * smallNumber;
          rY = Math.Atan2(k[0] / cosa, k[2] / cosa);
          rZ = Math.Atan2(i[1] / cosa, j[1] / cosa);
          ret[0, 1] = rX.RadiansToDegrees();
          ret[1, 1] = rY.RadiansToDegrees();
          ret[2, 1] = rZ.RadiansToDegrees();
          return ret;
        }
      }
      /// <summary>
      /// returns values in degrees.  both solutions, rx is[0,0] and [0,1] etc
      /// </summary>
      public double[,] rYrXrZbusted
      {
        get
        {
          double[,] ret = new double[3, 2];

          double rX = Math.Atan2(-k[1], Math.Sqrt(i[1] * i[1] + j[1] * j[1]));
          double cosa = Math.Cos(rX);
          double rY = 0.0;
          double rZ = 0.0;

          double smallNumber = .00001;

          if (Math.Abs(cosa) < smallNumber)
            cosa = Math.Sign(cosa) * smallNumber;

          rY = Math.Atan2(k[0] / cosa, k[2] / cosa);
          rZ = Math.Atan2(i[1] / cosa, j[1] / cosa);
          ret[0, 0] = rX.RadiansToDegrees();
          ret[1, 0] = rY.RadiansToDegrees();
          ret[2, 0] = rZ.RadiansToDegrees();

          rX = Math.Atan2(-k[1], -Math.Sqrt(i[1] * i[1] + j[1] * j[1]));
          cosa = Math.Cos(rX);
          if (Math.Abs(cosa) < smallNumber)
            cosa = Math.Sign(cosa) * smallNumber;
          rY = Math.Atan2(k[0] / cosa, k[2] / cosa);
          rZ = Math.Atan2(i[1] / cosa, j[1] / cosa);
          ret[0, 1] = rX.RadiansToDegrees();
          ret[1, 1] = rY.RadiansToDegrees();
          ret[2, 1] = rZ.RadiansToDegrees();
          return ret;
        }
      }

      public double[,] rZrYrX
      {
        get
        {
          double[,] ret = new double[3, 2];

          double rY = Math.Atan2(-i[2], Math.Sqrt(i[1] * i[1] + i[0] * i[0]));
          double cosb = Math.Cos(rY);
          double rX = 0.0;
          double rZ = 0.0;

          double smallNumber = .00001;

          if (Math.Abs(cosb) < smallNumber)
            cosb = Math.Sign(cosb) * smallNumber;
          rX = Math.Atan2(j[2] / cosb, k[2] / cosb);
          rZ = Math.Atan2(i[1] / cosb, i[0] / cosb);

          ret[0, 0] = rX.RadiansToDegrees();
          ret[1, 0] = rY.RadiansToDegrees();
          ret[2, 0] = rZ.RadiansToDegrees();

          rY = Math.Atan2(-i[2], Math.Sqrt(j[2] * j[2] + k[2] * k[2]));
          cosb = Math.Cos(rY);

          if (Math.Abs(cosb) < smallNumber)
            cosb = Math.Sign(cosb) * smallNumber;
          rX = Math.Atan2(j[2] / cosb, k[2] / cosb);
          rZ = Math.Atan2(i[1] / cosb, i[0] / cosb);
          ret[0, 1] = rX.RadiansToDegrees();
          ret[1, 1] = rY.RadiansToDegrees();
          ret[2, 1] = rZ.RadiansToDegrees();

          return ret;
        }
      }

      /// <summary>
      /// Returns a two element double array, where: <br/>
      /// B[0] = B in terms of 360 axis <br/>
      /// B[1] = B in terms of a +/- axis.
      /// Calculates the Btp Based on the k vector.  Units are degrees.
      /// </summary>
      public double B
      {
        get
        {
          return Math.Atan2(k[0], k[2]).RadiansToDegrees();
        }
      }

      /// <summary>
      /// Returns a two element double array, where: <br/>
      /// C[0] = C in terms of 360 axis <br/>
      /// C[1] = C in terms of a +/- axis.
      /// Calculates the Btp Based on the i vector.  Units are degrees.
      /// </summary>
      public double[] C
      {
        get
        {
          double[] kk = this.i;
          double dx = kk[0];
          double dy = kk[1];
          double[] ret = this.FindAngle(dy, dx);
          ret[0] = 360 - ret[0];
          ret[1] = -ret[1];
          return ret;
        }
      }

      private double[] FindAngle(double dy, double dz)
      {
        double r = Math.Sqrt(dz * dz + dy * dy);
        double[] ret = new double[2];
        if (r < 0.00001 && r > -0.00001)
        {
          ret[0] = ret[1] = 0;
          return ret;
        }
        if (dz >= 0)
        {
          double alpha = Math.Acos(dz / r);
          if (dy >= 0)
          {
            ret[0] = 360 - alpha.RadiansToDegrees();
            ret[1] = -alpha.RadiansToDegrees();
            return ret;
          }
          else
          {
            ret[0] = ret[1] = alpha.RadiansToDegrees();
            return ret;
          }
        }
        else
        {
          double alpha = Math.Acos(-dz / r);
          if (dy >= 0)
          {
            ret[0] = 180 + alpha.RadiansToDegrees();
            ret[1] = alpha.RadiansToDegrees() - 180;
            return ret;
          }
          else
          {
            ret[0] = ret[1] = 180 - alpha.RadiansToDegrees();
            return ret;
          }
        }
      }
      /// <summary>
      /// return one dimensional 3 element array representing R.
      /// </summary>
      public double[] R
      {
        get
        {
          double[] dog = new double[4];
          dog[0] = this.myM[0, 3];
          dog[1] = this.myM[1, 3];
          dog[2] = this.myM[2, 3];
          dog[3] = 1.0;
          return dog;
        }
        set
        {
          if (value.Length > 2)
          {
            myM[0, 3] = value[0];
            myM[1, 3] = value[1];
            myM[2, 3] = value[2];
            myM[3, 3] = 1;
          }
        }
      }
      /// <summary>
      /// Same as R, but in Russ's cool vector class.
      /// </summary>
      public Vector vR
      {
        get
        {
          return new Vector(R);
        }
        set
        {
          this.myM[0, 3] = value[0];
          this.myM[1, 3] = value[1];
          this.myM[2, 3] = value[2];
        }
      }
      /// <summary>
      /// Normalizes a one dimensional matrix a
      /// </summary>
      /// <param name="a">a one dimensional double matrix</param>
      /// <returns>a one dimensional double matrix with same number of elements as a</returns>
      public double[] Normalize(double[] a)
      {
        double hold = 0;
        double[] rtn = new double[a.Length];
        for (int ii = 0; ii < a.Length; ii++)
        {
          hold += a[ii] * a[ii];
        }
        hold = Math.Sqrt(hold);
        for (int ii = 0; ii < a.Length; ii++)
        {
          rtn[ii] = a[ii] / hold;
        }
        return rtn;
      }
    }
    #endregion

    #region Vector Class
    /// <summary>
    /// Copied straight from Russ's vector class.
    /// </summary>
    [Serializable]
    public class Vector
    {
      private static double d2r = Math.PI / 180.0;
      //private double[] V = new double[3];
      public double a;
      public double b;
      public double c;

      //below is for special error functions.
      public ePrimaryAxis PrimaryAxis = ePrimaryAxis.i;


      public enum ePrimaryAxis
      {
        i,
        j,
        k,
      }

      public double phi_deg
      {
        get
        {
          double d2r = Math.PI / 180.0;
          double ret = 0;
          switch (PrimaryAxis)
          {
            case ePrimaryAxis.i:
              ret = Math.Atan2(b, c) / d2r;
              break;
            case ePrimaryAxis.j:
              ret = Math.Atan2(a, c) / d2r;
              break;
            case ePrimaryAxis.k:
              ret = Math.Atan2(b, a) / d2r;
              break;
          }
          return ret;
        }
      }
      public double re
      {
        get
        {
          double ret = 0;
          double re2;
          switch (PrimaryAxis)
          {
            case ePrimaryAxis.i:
              re2 = 1 - a * a;
              if (re2 > 0)
                ret = Math.Sqrt(re2);
              break;
            case ePrimaryAxis.j:
              re2 = 1 - b * b;
              if (re2 > 0)
                ret = Math.Sqrt(re2);
              break;
            case ePrimaryAxis.k:
              re2 = 1 - c * c;
              if (re2 > 0)
                ret = Math.Sqrt(re2);
              break;
          }
          return ret;
        }
      }
      public Vector(double phi_deg, double re, ePrimaryAxis pa)
      {
        this.PrimaryAxis = pa;
        if (re >= 1.0)
          return;
        double phi = phi_deg * Math.PI / 180.0;
        double recosphi = re * Math.Cos(phi);
        double resinphi = re * Math.Sin(phi);
        double remainder = Math.Sqrt(1 - re * re);

        switch (pa)
        {
          case ePrimaryAxis.i:
            this.a = remainder;
            this.b = resinphi;
            this.c = recosphi;
            break;
          case ePrimaryAxis.j:
            this.a = resinphi;
            this.b = remainder;
            this.c = recosphi;
            break;
          case ePrimaryAxis.k:
            this.a = recosphi;
            this.b = resinphi;
            this.c = remainder;
            break;
        }
      }
      public Vector(double a, double b, double c)
      {
        this.a = a;
        this.b = b;
        this.c = c;


        //V[0] = this.a;
        //V[1] = this.b;
        //V[2] = this.c;
      }
      public Vector()
      {
        //V[0] = this.a;
        //V[1] = this.b;
        //V[2] = this.c;
      }
      public Vector(double[] input)
      {
        if (input.Length >= 3)
        {
          this.a = input[0];
          this.b = input[1];
          this.c = input[2];


          //V[0] = this.a;
          //V[1] = this.b;
          //V[2] = this.c;
        }
        else
        {
          throw new Exception("public Vector(double[] input)...error here");
        }

      }
      [XmlIgnore]
      public double Magnitude
      {
        get { return Math.Sqrt(a * a + b * b + c * c); }
      }
      public double i
      {
        get
        {
          double s = this.a;
          return s;
        }
        set
        {
          this.a = value;
        }
      }
      public double j
      {
        get
        {
          double s = this.b;
          return s;
        }
        set
        {
          this.b = value;
        }
      }
      public double k
      {
        get
        {
          double s = this.c;
          return s;
        }
        set
        {
          this.c = value;
        }
      }

      public static Vector operator *(double s, Vector a)
      {
        double[] doggy = new double[3];
        doggy[0] = a[0] * s;
        doggy[1] = a[1] * s;
        doggy[2] = a[2] * s;

        return new Vector(doggy);
      }
      public static Vector operator /(Vector a, double s)
      {
        double[] doggy = new double[3];
        doggy[0] = a[0] / s;
        doggy[1] = a[1] / s;
        doggy[2] = a[2] / s;

        return new Vector(doggy);
      }
      public static Vector operator +(Vector a, Vector b)
      {
        Vector Peanut = new Vector();
        Peanut[0] = a[0] + b[0];
        Peanut[1] = a[1] + b[1];
        Peanut[2] = a[2] + b[2];
        return Peanut;
      }
      public static Vector operator -(Vector a, Vector b)
      {
        Vector Peanut = new Vector();
        Peanut[0] = a[0] - b[0];
        Peanut[1] = a[1] - b[1];
        Peanut[2] = a[2] - b[2];
        return Peanut;
      }
      public static double dot(Vector a, Vector b)
      {
        double s;
        s = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        return s;
      }
      public static Vector multiply_by_scaler(Vector a, double Mag)
      {
        Vector s;
        s = new Vector(a.i * Mag, a.j * Mag, a.k * Mag);
        return s;
      }

      public Vector Cross(Vector dog)
      {
        double[] V = { a, b, c };
        Vector Peanut = new Vector(
                      V[1] * dog[2] - dog[1] * V[2],
                      V[2] * dog[0] - dog[2] * V[0],
                      V[0] * dog[1] - dog[0] * V[1]);
        return Peanut;
      }
      public double this[int index]
      {
        get
        {
          double[] V = { a, b, c };
          return V[index];
        }
        set
        {
          double[] V = { a, b, c };
          V[index] = value;
          a = V[0];
          b = V[1];
          c = V[2];
        }
      }
      [XmlIgnore]
      public double[] R
      {
        get
        {
          double[] r = { a, b, c, 1.0 };
          return r;
        }
        set
        {
          a = value[0];
          b = value[1];
          c = value[2];
          //V[3] = 1.0;
        }
      }

      public void NormalizeMe()
      {
        double m = Magnitude;
        a /= m;
        b /= m;
        c /= m;
      }

      /// <summary>
      /// Normalizes the vector.  
      /// Returns rx and ry as if this is the K vector.
      /// </summary>
      /// <returns>array of two elements, rx and ry respectively in degrees</returns>
      public double[] GetRxRy()
      {
        return GetRxRy(this.i, this.j, this.k);
      }

      public static double[] GetRxRy(double i, double j, double k)
      {
        Vector vv = new Vector(i, j, k);
        vv.NormalizeMe();
        double rx = -Math.Atan2(vv.j, vv.k) / d2r;
        double zp = Math.Sqrt(vv.j * vv.j + vv.k * vv.k);

        double ry = Math.Atan2(vv.i, zp) / d2r;

        return new double[] { rx, ry };      
      }
      /// <summary>
      /// Returns the included angle between this vector and RayIn.  Absolute value in degrees.
      /// </summary>
      /// <param name="RayIn"></param>
      /// <returns>Absolute Value in Degrees</returns>
      public double IncludedAngle(Vector RayIn)
      {
        Electroimpact.LinearAlgebra.Vector n = this.Cross(RayIn);
        double theta = (Math.Asin(n.Magnitude / (this.Magnitude * RayIn.Magnitude))) / d2r;
        return Math.Abs(theta);
      }

      /// <summary>
      /// returns signed angle in degrees, 2d only.  Uses j and k components
      /// 
      /// </summary>
      /// <param name="RayIn"></param>
      /// <returns></returns>
      public double SignedAngle2D(Vector RayIn)
      {
        double crossed = j * RayIn.k - k * RayIn.j;
        double angle = Math.Atan2(crossed, dot(this, RayIn)) * 180.0 / Math.PI;
        return angle;
      }

      public override string ToString()
      {
        return ToString(1.0);
      }
      public string ToString(double scaler)
      {
        double[] V = { a, b, c };
        return (V[0] * scaler).ToString("F3") + " " + (V[1] * scaler).ToString("F3") + " " + (V[2] * scaler).ToString("F3") + " " + (Magnitude * scaler).ToString("F3");
      }
    }
    #endregion

    #region Standard Denavit-Hartenberg Matrix
    /// <summary>
    /// The Denavit-Hartenberg Matrix used in all of ToddR's kinematics derivations.
    /// Please see Todd's Redbook for more information on this or "Robotic Engineering an Integrated Approach".
    /// <br /><br />
    /// Derivation of arguments for the DH matrix must follow the following procedure:<br />
    /// 1) Find Theta (rotation about Z)<br />
    /// 2) Find d (translate along Z)<br />
    /// 3) Find a (translate along X)<br />
    /// 4) Find Alpha (rotation about X)<br />
    /// </summary>
    public class DHMatrix
    {
      #region Members
      protected double[,] mySelf = new double[4, 4];
      protected double myTheta;
      protected double myAlpha;
      protected double myA;
      protected double myD;

      public bool RecalcReqd = true;
      public bool RecalcTheta = true;
      public bool RecalcAlpha = true;
      public bool RecalcA = true;
      public bool RecalcD = true;
      private bool AlwaysRecalc;

      public class cDHStrings
      {
        public string Theta;
        public string d;
        public string a;
        public string Alpha;
      }

      public cDHStrings myStrings = new cDHStrings();
      #endregion

      #region Methods

      public DHMatrix(cDHStrings FromStrings, Electroimpact.StringCalc.cStringCalc Calculator, bool AlwaysRecalculate)
      {
        this.myStrings = FromStrings;
        this.mySelf = this.formDHMatrix(Calculator);
        AlwaysRecalc = AlwaysRecalculate;
      }
      public void ForceRecalculation()
      {
        RecalcReqd = true;
        RecalcTheta = true;
        RecalcAlpha = true;
        RecalcA = true;
        RecalcD = true;
      }
      #endregion

      #region Properties
      /// <summary>
      /// Returns properly formed Denavit - Hartenberg Matrix calculated from Theta, Alpha, a and d 
      /// which is a 4x4 double array.
      /// </summary>
      public double[,] getDH(Electroimpact.StringCalc.cStringCalc Calculator)
      {
        if (RecalcReqd || AlwaysRecalc)
          return this.formDHMatrix(Calculator);
        else
          return this.Self;
      }
      #endregion

      #region private Methods
      //virtual protected void ReformMe()
      //{
      //  this.mySelf = formDHMatrix(this.myTheta, this.myAlpha, this.myA, this.myD);
      //}
      protected double[,] formDHMatrix(Electroimpact.StringCalc.cStringCalc Calculator)// double Theta, double Alpha, double a, double d)
      {
        if (AlwaysRecalc)
        {
          RecalcTheta = true;
          RecalcAlpha = true;
          RecalcA = true;
          RecalcD = true;
        }
        if (RecalcTheta)
          this.myTheta = Calculator.Degrees ? Calculator.SimpleCalc(myStrings.Theta) * Math.PI / 180.0 : Calculator.SimpleCalc(myStrings.Theta);
        if (RecalcAlpha)
          this.myAlpha = Calculator.Degrees ? Calculator.SimpleCalc(myStrings.Alpha) * Math.PI / 180.0 : Calculator.SimpleCalc(myStrings.Alpha);
        if (RecalcD)
          this.myD = Calculator.SimpleCalc(myStrings.d);
        if (RecalcA)
          this.myA = Calculator.SimpleCalc(myStrings.a);

        double[,] dog = new double[4, 4];
        double sinAlpha = Math.Sin(this.myAlpha);
        double cosAlpha = Math.Cos(this.myAlpha);
        double sinTheta = Math.Sin(this.myTheta);
        double cosTheta = Math.Cos(this.myTheta);

        dog[0, 0] = cosTheta;
        dog[0, 1] = -cosAlpha * sinTheta;
        dog[0, 2] = sinAlpha * sinTheta;
        dog[0, 3] = this.myA * cosTheta;

        dog[1, 0] = sinTheta;
        dog[1, 1] = cosAlpha * cosTheta;
        dog[1, 2] = -cosTheta * sinAlpha;
        dog[1, 3] = this.myA * sinTheta;

        dog[2, 0] = 0;
        dog[2, 1] = sinAlpha;
        dog[2, 2] = cosAlpha;
        dog[2, 3] = this.myD;

        dog[3, 0] = 0;
        dog[3, 1] = 0;
        dog[3, 2] = 0;
        dog[3, 3] = 1;

        this.mySelf = dog;

        RecalcReqd = false;
        RecalcA = false;
        RecalcTheta = false;
        RecalcAlpha = false;
        RecalcD = false;

        return dog;
      }

      protected double[,] Self
      {
        get { return this.mySelf; }
      }
      #endregion
    }
    #endregion

    #region Standard 6dof string enabled matrix...arguments X,Y,Z,rX,rY,rZ
    /// <summary>
    /// The Denavit-Hartenberg Matrix used in all of ToddR's kinematics derivations.
    /// Please see Todd's Redbook for more information on this or "Robotic Engineering an Integrated Approach".
    /// <br /><br />
    /// Derivation of arguments for the DH matrix must follow the following procedure:<br />
    /// 1) Find Theta (rotation about Z)<br />
    /// 2) Find d (translate along Z)<br />
    /// 3) Find a (translate along X)<br />
    /// 4) Find Alpha (rotation about X)<br />
    /// </summary>
    public class string6DOF
    {
      #region Members
      protected double[,] mySelf = new double[4, 4];
      protected double myX;
      protected double myY;
      protected double myZ;
      protected double myrX;
      protected double myrY;
      protected double myrZ;
      protected Electroimpact.StringCalc.cStringCalc myCalculator;

      public bool RecalcReqd = true;
      public bool RecalcX = true;
      public bool RecalcY = true;
      public bool RecalcZ = true;
      public bool RecalcrX = true;
      public bool RecalcrY = true;
      public bool RecalcrZ = true;
      private bool AlwaysRecalc;

      public class cArgStrings
      {
        public string X;
        public string Y;
        public string Z;
        public string rX;
        public string rY;
        public string rZ;
      }

      public cArgStrings myStrings = new cArgStrings();
      #endregion

      #region Methods

      public string6DOF(cArgStrings FromStrings, Electroimpact.StringCalc.cStringCalc Calculator, bool AlwaysRecalculate)
      {
        this.myStrings = FromStrings;
        this.mySelf = this.formMatrix(Calculator);
        AlwaysRecalc = AlwaysRecalculate;
        myCalculator = Calculator;
      }
      public void ForceRecalculation()
      {
        RecalcReqd = true;
        RecalcX = true;
        RecalcY = true;
        RecalcZ = true;
        RecalcrX = true;
        RecalcrY = true;
        RecalcrZ = true;
      }

      public double[,] DotMe(double[,] Min)
      {
        c6dof me = new c6dof(formMatrix(myCalculator));
        return me.DotMe(Min);
      }

      #endregion

      #region Properties
      /// <summary>
      /// Returns properly formed Denavit - Hartenberg Matrix calculated from Theta, Alpha, a and d 
      /// which is a 4x4 double array.
      /// </summary>
      public double[,] getDH(Electroimpact.StringCalc.cStringCalc Calculator)
      {
        if (RecalcReqd || AlwaysRecalc)
          return this.formMatrix(Calculator);
        else
          return this.Self;
      }
      #endregion

      #region private Methods
      //virtual protected void ReformMe()
      //{
      //  this.mySelf = formDHMatrix(this.myTheta, this.myAlpha, this.myA, this.myD);
      //}
      protected double[,] formMatrix(Electroimpact.StringCalc.cStringCalc Calculator)// double Theta, double Alpha, double a, double d)
      {
        if (AlwaysRecalc)
        {
          RecalcX = true;
          RecalcY = true;
          RecalcZ = true;
          RecalcrX = true;
          RecalcrY = true;
          RecalcrZ = true;
        }

        if (RecalcX)
          this.myX = Calculator.SimpleCalc(myStrings.X);
        if (RecalcY)
          this.myY = Calculator.SimpleCalc(myStrings.Y);
        if (RecalcZ)
          this.myZ = Calculator.SimpleCalc(myStrings.Z);
        if (RecalcrX)
          this.myrX = Calculator.Degrees ? Calculator.SimpleCalc(myStrings.rX).DegreesToRadians() : Calculator.SimpleCalc(myStrings.rX);
        if (RecalcrY)
          this.myrY = Calculator.Degrees ? Calculator.SimpleCalc(myStrings.rY).DegreesToRadians() : Calculator.SimpleCalc(myStrings.rY);
        if (RecalcrZ)
          this.myrZ = Calculator.Degrees ? Calculator.SimpleCalc(myStrings.rZ).DegreesToRadians() : Calculator.SimpleCalc(myStrings.rZ);

        c6dof m = new c6dof(this.myX, this.myY, this.myZ, this.myrX, this.myrY, this.myrZ);

        this.mySelf = m.GetMatrix();

        RecalcReqd = false;
        RecalcX = false;
        RecalcY = false;
        RecalcZ = false;
        RecalcrX = false;
        RecalcrY = false;
        RecalcrZ = false;

        return m.GetMatrix();
      }

      protected double[,] Self
      {
        get { return this.mySelf; }
      }
      #endregion
    }
    #endregion

    #region Single Axis Rotation and Translation Matricies
    /// <summary>
    /// A matrix class describing a rotation of coordinates about X.
    /// </summary>
    public class RxMatrix
    {
      #region Members
      double phi = 0;
      double[,] mySelf = new double[4, 4];
      #endregion

      #region Constructor
      /// <summary>
      /// Default constructor.
      /// </summary>
      public RxMatrix()
      {
        phi = 0;
        this.mySelf = formMatrix(this.phi);
      }
      /// <summary>
      /// Specific Constructor.
      /// </summary>
      /// <param name="Phi">double rotation about X amount in Radians.</param>
      public RxMatrix(double Phi)
      {
        phi = Phi;
        this.mySelf = formMatrix(this.phi);
      }
      #endregion

      #region Properties
      /// <summary>
      /// returns the 4x4 matrix contained in this class.
      /// </summary>
      public double[,] GetMatrix
      {
        get { return this.mySelf; }
      }
      /// <summary>
      /// gets and sets the amount of rotation about axis.  This number should be in radians.
      /// </summary>
      public double Phi
      {
        get { return this.phi; }
        set
        {
          this.phi = value;
          this.ReformMe();
        }
      }
      #endregion

      #region Methods
      public double[,] DotMe(double[,] Input)
      {
        cMatrix m = new cMatrix(this.mySelf);
        return m.DotMe(Input);
      }

      public double[] DotMe(double[] Input)
      {
        cMatrix m = new cMatrix(this.mySelf);
        return m.DotMe(Input);
      }

      public override string ToString()
      {
        cMatrix m = new cMatrix(this.mySelf);
        return m.ToString();
      }

      public double[] GetIJK()
      {
        Vector v = new Vector(0, 0, 1);

        double[] ret = this.DotMe(v.R);
        return ret;
      }
      #endregion

      #region private Methods
      private void ReformMe()
      {
        this.mySelf = formMatrix(this.phi);
      }

      private double[,] formMatrix(double Phi)
      {
        double[,] dog = new double[4, 4];
        double sinPhi = Math.Sin(Phi);
        double cosPhi = Math.Cos(Phi);

        dog[0, 0] = 1;
        dog[0, 1] = 0;
        dog[0, 2] = 0;
        dog[0, 3] = 0;

        dog[1, 0] = 0;
        dog[1, 1] = cosPhi;
        dog[1, 2] = -sinPhi;
        dog[1, 3] = 0;

        dog[2, 0] = 0;
        dog[2, 1] = sinPhi;
        dog[2, 2] = cosPhi;
        dog[2, 3] = 0;

        dog[3, 0] = 0;
        dog[3, 1] = 0;
        dog[3, 2] = 0;
        dog[3, 3] = 1;

        return dog;
      }
      #endregion
    }

    /// <summary>
    /// A matrix class describing a rotation of coordinates about Y.
    /// </summary>
    public class RyMatrix
    {
      #region Members
      double phi = 0;
      double[,] mySelf = new double[4, 4];
      #endregion

      #region Constructor
      /// <summary>
      /// Default constructor.
      /// </summary>
      public RyMatrix()
      {
        phi = 0;
        this.mySelf = formMatrix(this.phi);
      }
      /// <summary>
      /// Specific Constructor.
      /// </summary>
      /// <param name="Phi">double rotation about Y amount.</param>
      public RyMatrix(double Phi)
      {
        phi = Phi;
        this.mySelf = formMatrix(this.phi);
      }
      #endregion

      #region Properties
      /// <summary>
      /// returns the 4x4 matrix contained in this class.
      /// </summary>
      public double[,] GetMatrix
      {
        get { return this.mySelf; }
      }
      /// <summary>
      /// gets and sets the amount of rotation about axis
      /// </summary>
      public double Phi
      {
        get { return this.phi; }
        set
        {
          this.phi = value;
          this.ReformMe();
        }
      }
      #endregion
      #region Methods
      public double[,] DotMe(double[,] Input)
      {
        cMatrix m = new cMatrix(this.mySelf);
        return m.DotMe(Input);
      }

      public double[] DotMe(double[] Input)
      {
        cMatrix m = new cMatrix(this.mySelf);
        return m.DotMe(Input);
      }

      public override string ToString()
      {
        cMatrix m = new cMatrix(this.mySelf);
        return m.ToString();
      }
      #endregion
      #region private Methods
      private void ReformMe()
      {
        this.mySelf = formMatrix(this.phi);
      }

      private double[,] formMatrix(double Phi)
      {
        double[,] dog = new double[4, 4];
        double sinPhi = Math.Sin(Phi);
        double cosPhi = Math.Cos(Phi);

        dog[0, 0] = cosPhi;
        dog[0, 1] = 0;
        dog[0, 2] = sinPhi;
        dog[0, 3] = 0;

        dog[1, 0] = 0;
        dog[1, 1] = 1;
        dog[1, 2] = 0;
        dog[1, 3] = 0;

        dog[2, 0] = -sinPhi;
        dog[2, 1] = 0;
        dog[2, 2] = cosPhi;
        dog[2, 3] = 0;

        dog[3, 0] = 0;
        dog[3, 1] = 0;
        dog[3, 2] = 0;
        dog[3, 3] = 1;

        return dog;
      }
      #endregion

      public double[,] Inverse()
      {
        Electroimpact.LinearAlgebra.cMatrix m = new cMatrix(GetMatrix);
        return m.Inverse();
      }

    }

    /// <summary>
    /// A matrix class describing a rotation of coordinates about Z.
    /// </summary>
    public class RzMatrix
    {
      #region Members
      double phi = 0;
      double[,] mySelf = new double[4, 4];
      #endregion

      #region Constructor
      /// <summary>
      /// Default constructor.
      /// </summary>
      public RzMatrix()
      {
        phi = 0;
        this.mySelf = formMatrix(this.phi);
      }
      /// <summary>
      /// Specific Constructor.
      /// </summary>
      /// <param name="Phi">double rotation about Z in radians.</param>
      public RzMatrix(double Phi)
      {
        phi = Phi;
        this.mySelf = formMatrix(this.phi);
      }
      #endregion

      #region Properties
      /// <summary>
      /// returns the 4x4 matrix contained in this class.
      /// </summary>
      public double[,] GetMatrix
      {
        get { return this.mySelf; }
      }
      /// <summary>
      /// gets and sets the amount of rotation about axis
      /// </summary>
      public double Phi
      {
        get { return this.phi; }
        set
        {
          this.phi = value;
          this.ReformMe();
        }
      }
      #endregion

      #region Methods
      public double[,] DotMe(double[,] Input)
      {
        cMatrix m = new cMatrix(this.mySelf);
        return m.DotMe(Input);
      }

      public double[] DotMe(double[] Input)
      {
        cMatrix m = new cMatrix(this.mySelf);
        return m.DotMe(Input);
      }

      public override string ToString()
      {
        cMatrix m = new cMatrix(this.mySelf);
        return m.ToString();
      }
      #endregion

      #region private Methods
      private void ReformMe()
      {
        this.mySelf = formMatrix(this.phi);
      }

      private double[,] formMatrix(double Phi)
      {
        double[,] dog = new double[4, 4];
        double sinPhi = Math.Sin(Phi);
        double cosPhi = Math.Cos(Phi);

        dog[0, 0] = cosPhi;
        dog[0, 1] = -sinPhi;
        dog[0, 2] = 0;
        dog[0, 3] = 0;

        dog[1, 0] = sinPhi;
        dog[1, 1] = cosPhi;
        dog[1, 2] = 0;
        dog[1, 3] = 0;

        dog[2, 0] = 0;
        dog[2, 1] = 0;
        dog[2, 2] = 1;
        dog[2, 3] = 0;

        dog[3, 0] = 0;
        dog[3, 1] = 0;
        dog[3, 2] = 0;
        dog[3, 3] = 1;

        return dog;
      }
      #endregion

      public double[,] Inverse()
      {
        cMatrix m = new cMatrix(GetMatrix);
        return m.Inverse();
      }
    }

    public class RrMatrix
    {
      #region Members
      double x = 0;
      double y = 0;
      double z = 0;
      double[,] mySelf = new double[4, 4];
      #endregion

      #region Constructor
      /// <summary>
      /// Default constructor.
      /// </summary>
      public RrMatrix()
      {
        this.mySelf = formMatrix(this.x, this.y, this.z);
      }
      /// <summary>
      /// Translation Matrix Constructor
      /// </summary>
      /// <param name="X">double value representing translation along the X axis.</param>
      /// <param name="Y">double value representing translation along the Y axis.</param>
      /// <param name="Z">double value representing translation along the Z axis.</param>
      public RrMatrix(double X, double Y, double Z)
      {
        this.x = X;
        this.y = Y;
        this.z = Z;
        this.mySelf = formMatrix(this.x, this.y, this.z);
      }
      #endregion

      #region Properties
      /// <summary>
      /// returns the 4x4 matrix contained in this class.
      /// </summary>
      public double[,] GetMatrix
      {
        get { return this.mySelf; }
      }
      /// <summary>
      /// double value representing translation along the X axis.
      /// </summary>
      public double X
      {
        get { return this.x; }
        set { this.x = value; ReformMe(); }
      }
      /// <summary>
      /// double value representing translation along the Y axis.
      /// </summary>
      public double Y
      {
        get { return this.y; }
        set { this.y = value; ReformMe(); }
      }
      /// <summary>
      /// double value representing translation along the Z axis.
      /// </summary>
      public double Z
      {
        get { return this.z; }
        set { this.z = value; ReformMe(); }
      }
      #endregion

      #region private Methods
      private void ReformMe()
      {
        this.mySelf = formMatrix(this.x, this.y, this.z);
      }

      private double[,] formMatrix(double x, double y, double z)
      {
        double[,] dog = new double[4, 4];

        dog[0, 0] = 1;
        dog[0, 1] = 0;
        dog[0, 2] = 0;
        dog[0, 3] = x;

        dog[1, 0] = 0;
        dog[1, 1] = 1;
        dog[1, 2] = 0;
        dog[1, 3] = y;

        dog[2, 0] = 0;
        dog[2, 1] = 0;
        dog[2, 2] = 1;
        dog[2, 3] = z;

        dog[3, 0] = 0;
        dog[3, 1] = 0;
        dog[3, 2] = 0;
        dog[3, 3] = 1;

        return dog;
      }
      #endregion
    }
    #endregion



    public class c6dof_robot : c6dof
    {
      protected override double[,] BuildMatrix()
      {

        double[,] matrix = new double[4, 4];

        double sinrx = Math.Sin(args.rX_radians);
        double sinry = Math.Sin(args.rY_radians);
        double sinrz = Math.Sin(args.rZ_radians);
        double cosrx = Math.Cos(args.rX_radians);
        double cosry = Math.Cos(args.rY_radians);
        double cosrz = Math.Cos(args.rZ_radians);

        matrix[0, 0] = cosry * cosrz;
        matrix[1, 0] = cosry * sinrz;
        matrix[2, 0] = -sinry;
        matrix[3, 0] = 0.0;

        matrix[0, 1] = cosrz * sinrx * sinry - cosrx * sinrz;
        matrix[1, 1] = cosrx * cosrz + sinrx * sinry * sinrz;
        matrix[2, 1] = cosry * sinrx;
        matrix[3, 1] = 0.0;

        matrix[0, 2] = sinrx * sinrz + cosrx * cosrz * sinry;
        matrix[1, 2] = cosrx * sinry * sinrz - cosrz * sinrx;
        matrix[2, 2] = cosrx * cosry;
        matrix[3, 2] = 0.0;

        matrix[0, 3] = args.X;
        matrix[1, 3] = args.Y;
        matrix[2, 3] = args.Z;
        matrix[3, 3] = 1.0;

        return matrix;
      }

      #region constructors
      public c6dof_robot(c6DOF_Args args_in)
      {
        //rx = new RxMatrix(args_in.rX_radians);
        //ry = new RyMatrix(args_in.rY_radians);
        //rz = new RzMatrix(args_in.rZ_radians);
        ////        rr = new RrMatrix(X, Y, Z);
        //rr = new RrMatrix(args_in.X, args_in.Y, args_in.Z);

        this.args = new c6DOF_Args(args_in);
      }

      public c6dof_robot(double X, double Y, double Z, double rXradians, double rYradians, double rZradians)
      {
        //rx = new RxMatrix(rXradians);
        //ry = new RyMatrix(rYradians);
        //rz = new RzMatrix(rZradians);
        //rr = new RrMatrix(X, Y, Z);

        this.args = new c6DOF_Args(X, Y, Z, rXradians, rYradians, rZradians);
      }
      /// <summary>
      /// Base constructor
      /// </summary>
      public c6dof_robot()
      {
        //rx = new RxMatrix(0);
        //ry = new RyMatrix(0);
        //rz = new RzMatrix(0);
        //rr = new RrMatrix(0, 0, 0);

        this.args = new c6DOF_Args();
      }
      /// <summary>
      /// Copy Constructor
      /// </summary>
      /// <param name="Input">c6dof class</param>
      public c6dof_robot(c6dof Input)
      {
        //rx = new RxMatrix(Input.rX);
        //ry = new RyMatrix(Input.rY);
        //rz = new RzMatrix(Input.rZ);
        //rr = new RrMatrix(Input.X, Input.Y, Input.Z);

        this.args = new c6DOF_Args(Input.X, Input.Y, Input.Z, Input.rX, Input.rY, Input.rZ);
      }
      public c6dof_robot(double[,] Input)
      {
        Electroimpact.LinearAlgebra.cMatrix m = new cMatrix(Input);
        //rx = new RxMatrix();
        //ry = new RyMatrix();
        //rz = new RzMatrix();
        //rr = new RrMatrix();
        //rX = m.rXrYrZ[0, 0].DegreesToRadians();
        //rY = m.rXrYrZ[1, 0].DegreesToRadians();
        //rZ = m.rXrYrZ[2, 0].DegreesToRadians();
        //X = m.X;
        //Y = m.Y;
        //Z = m.Z;

        this.args = new c6DOF_Args(m.X, m.Y, m.Z, m.rZrYrX[0, 0].DegreesToRadians(), m.rZrYrX[1, 0].DegreesToRadians(), m.rZrYrX[2, 0].DegreesToRadians());
      }
      #endregion


      public override void Make6DoF(double[,] StateMatrix)
      {
        Electroimpact.LinearAlgebra.cMatrix m = new cMatrix(StateMatrix);
        X = m.X;
        Y = m.Y;
        Z = m.Z;
        rX = m.rZrYrX[0, 0].DegreesToRadians();
        rY = m.rZrYrX[1, 0].DegreesToRadians();
        rZ = m.rZrYrX[2, 0].DegreesToRadians();
      }
    }


    #region 6dof matrix
    /// <summary>
    /// Displace linear, rotate about X, then Y then Z
    /// </summary>
    public class c6dof
    {

      protected c6DOF_Args args = new c6DOF_Args();

      //input values in radians please

      #region constructors
      public c6dof(c6DOF_Args args_in)
      {
        //rx = new RxMatrix(args_in.rX_radians);
        //ry = new RyMatrix(args_in.rY_radians);
        //rz = new RzMatrix(args_in.rZ_radians);
        ////        rr = new RrMatrix(X, Y, Z);
        //rr = new RrMatrix(args_in.X, args_in.Y, args_in.Z);

        this.args = new c6DOF_Args(args_in);
      }

      public c6dof(double X, double Y, double Z, double rXradians, double rYradians, double rZradians)
      {
        //rx = new RxMatrix(rXradians);
        //ry = new RyMatrix(rYradians);
        //rz = new RzMatrix(rZradians);
        //rr = new RrMatrix(X, Y, Z);

        this.args = new c6DOF_Args(X, Y, Z, rXradians, rYradians, rZradians);
      }
      /// <summary>
      /// Base constructor
      /// </summary>
      public c6dof()
      {
        //rx = new RxMatrix(0);
        //ry = new RyMatrix(0);
        //rz = new RzMatrix(0);
        //rr = new RrMatrix(0, 0, 0);

        this.args = new c6DOF_Args();
      }
      /// <summary>
      /// Copy Constructor
      /// </summary>
      /// <param name="Input">c6dof class</param>
      public c6dof(c6dof Input)
      {
        //rx = new RxMatrix(Input.rX);
        //ry = new RyMatrix(Input.rY);
        //rz = new RzMatrix(Input.rZ);
        //rr = new RrMatrix(Input.X, Input.Y, Input.Z);

        this.args = new c6DOF_Args(Input.X, Input.Y, Input.Z, Input.rX, Input.rY, Input.rZ);
      }
      public c6dof(double[,] Input)
      {
        Electroimpact.LinearAlgebra.cMatrix m = new cMatrix(Input);
        //rx = new RxMatrix();
        //ry = new RyMatrix();
        //rz = new RzMatrix();
        //rr = new RrMatrix();
        //rX = m.rXrYrZ[0, 0].DegreesToRadians();
        //rY = m.rXrYrZ[1, 0].DegreesToRadians();
        //rZ = m.rXrYrZ[2, 0].DegreesToRadians();
        //X = m.X;
        //Y = m.Y;
        //Z = m.Z;

        this.args = new c6DOF_Args(m.X, m.Y, m.Z, m.rXrYrZ[0,0].DegreesToRadians(), m.rXrYrZ[1,0].DegreesToRadians(), m.rXrYrZ[2,0].DegreesToRadians());

      }
      #endregion


      protected virtual double[,] BuildMatrix()
      {

        double[,] matrix = new double[4, 4];

        double sinrx = Math.Sin(args.rX_radians);
        double sinry = Math.Sin(args.rY_radians);
        double sinrz = Math.Sin(args.rZ_radians);
        double cosrx = Math.Cos(args.rX_radians);
        double cosry = Math.Cos(args.rY_radians);
        double cosrz = Math.Cos(args.rZ_radians);

        matrix[0, 0] = cosry * cosrz;
        matrix[1, 0] = cosrx * sinrz + cosrz * sinrx * sinry;
        matrix[2, 0] = sinrx * sinrz - cosrx * cosrz * sinry;
        matrix[3,0] = 0.0;

        matrix[0, 1] = -cosry * sinrz;
        matrix[1, 1] = cosrx * cosrz - sinrx * sinry * sinrz;
        matrix[2, 1] = cosrz * sinrx + cosrx * sinry * sinrz;
        matrix[3, 1] = 0.0;

        matrix[0, 2] = sinry;
        matrix[1, 2] = -cosry * sinrx;
        matrix[2, 2] = cosrx * cosry;
        matrix[3, 2] = 0.0;

        matrix[0, 3] = args.X;
        matrix[1, 3] = args.Y;
        matrix[2, 3] = args.Z;
        matrix[3, 3] = 1.0;

        return matrix;
      }

      protected virtual double[,] BuildMatrix_rZrYrX(double X_in, double Y_in, double Z_in, double rXradians, double rYradians, double rZradians)
      {
        double[,] matrix = new double[4, 4];

        double sinrx = Math.Sin(rXradians);//args.rX_radians);
        double sinry = Math.Sin(rYradians);//args.rY_radians);
        double sinrz = Math.Sin(rZradians);//args.rZ_radians);
        double cosrx = Math.Cos(rXradians);//args.rX_radians);
        double cosry = Math.Cos(rYradians);//args.rY_radians);
        double cosrz = Math.Cos(rZradians);//args.rZ_radians);

        matrix[0, 0] = cosrz * cosry;
        matrix[1, 0] = sinrz * cosry; 
        matrix[2, 0] = -sinry;
        matrix[3, 0] = 0;

        matrix[0, 1] = -sinrz*cosrx + cosrz*sinry*sinrx; 
        matrix[1, 1] = cosrz * cosrx + sinrz * sinry * sinrx;
        matrix[2, 1] = cosry * sinrx; 
        matrix[3, 1] = 0;

        matrix[0, 2] = sinrz * sinrx + cosrz * sinry * cosrx;
        matrix[1, 2] = -cosrz*sinrx + sinrz*sinry*cosrx; 
        matrix[2, 2] = cosry * cosrx; 
        matrix[3, 2] = 0;

        matrix[0, 3] = X_in;// args.X;
        matrix[1, 3] = Y_in;// args.Y;
        matrix[2, 3] = Z_in;// args.Z;
        matrix[3, 3] = 1.0;

        return matrix;


      }

      public override string ToString()
      {
        string ret = "";
        ret += "X: " + args.X.ToString("F3") + "\n";
        ret += "Y: " + args.Y.ToString("F3") + "\n";
        ret += "Z: " + args.Z.ToString("F3") + "\n";
        ret += "rX: " + (args.rX).ToString("F3") + "\n";
        ret += "rY: " + (args.rY).ToString("F3") + "\n";
        ret += "rZ: " + (args.rZ).ToString("F3");
        return ret;
      }

      /// <summary>
      /// Gets the matrix based on the arguments below.  Pass in radians please.
      /// </summary>
      /// <param name="X">linear distance</param>
      /// <param name="Y">linear distance</param>
      /// <param name="Z">linear distance</param>
      /// <param name="rX">radians</param>
      /// <param name="rY">radians</param>
      /// <param name="rZ">radians</param>
      /// <returns></returns>
      public double[,] GetMatrix(double X, double Y, double Z, double rXradians, double rYradians, double rZradians)
      {
        args = new c6DOF_Args(X, Y, Z, rXradians, rYradians, rZradians);
        return BuildMatrix();
      }
      /// <summary>
      /// Returns a state matrix represending this transform
      /// </summary>
      /// <returns></returns>
      public double[,] GetMatrix()
      {
        return BuildMatrix();
      }

      /// <summary>
      /// Gets the matrix based on the arguments below.  Pass in radians please.
      /// </summary>
      /// <param name="X">linear distance</param>
      /// <param name="Y">linear distance</param>
      /// <param name="Z">linear distance</param>
      /// <param name="rX">radians</param>
      /// <param name="rY">radians</param>
      /// <param name="rZ">radians</param>
      /// <returns></returns>
      public void MakeMatrix_rZrYrX(double X_in, double Y_in, double Z_in, double rXradians, double rYradians, double rZradians)
      {
        //args = new c6DOF_Args(X, Y, Z, rXradians, rYradians, rZradians);

        double[,] temp = BuildMatrix_rZrYrX(X_in, Y_in, Z_in, rXradians, rYradians, rZradians);
        Electroimpact.LinearAlgebra.cMatrix m = new cMatrix(temp);
        X = m.X;
        Y = m.Y;
        Z = m.Z;
        rX = m.rXrYrZ[0, 0].DegreesToRadians();
        rY = m.rXrYrZ[1, 0].DegreesToRadians();
        rZ = m.rXrYrZ[2, 0].DegreesToRadians();

        args = new c6DOF_Args(X,Y,Z,rX,rY,rZ); //set the arguments the matrix is based on

        //double Output = {X,Y,Z,rX,rY,rZ);
        

        return; 
      }     

      public cMatrix myMatixClass
      {
        get
        {
          return new cMatrix(GetMatrix());
        }
      }

      /// <summary>
      /// Returns only the roatation portion of the state matrix representing this transform
      /// </summary>
      /// <returns></returns>
      public double[,] GetOrientationMatrix()
      {
        double[,] ret = BuildMatrix();
        ret[0, 3] = 0.0;
        ret[1, 3] = 0.0;
        ret[2, 3] = 0.0;
        return ret;
      }
      /// <summary>
      /// a double value
      /// </summary>
      public double X
      {
        get { return args.X; }
        set { args.X = value; }
      }
      /// <summary>
      /// a double value
      /// </summary>
      public double Y
      {
        get { return args.Y; }
        set { args.Y = value; }
      }
      /// <summary>
      /// a double value
      /// </summary>
      public double Z
      {
        get { return args.Z; }
        set { args.Z = value; }
      }
      /// <summary>
      /// Sets or returns a value in radians
      /// </summary>
      public double rX
      {
        get { return args.rX_radians; }
        set { args.rX_radians = value; }
      }
      /// <summary>
      /// Sets or returns a value in radians
      /// </summary>      
      public double rY
      {
        get { return args.rY_radians; }
        set { args.rY_radians = value; }
      }
      /// <summary>
      /// Sets or returns a value in radians
      /// </summary>      
      public double rZ
      {
        get { return args.rZ_radians; }
        set { args.rZ_radians = value; }
      }

      /// <summary>
      /// Copy function.  I haven't overrided =, == etc.  Too much work..use this.
      /// </summary>
      /// <param name="sixDOF"></param>
      public void Copy(c6dof sixDOF)
      {
        X = sixDOF.X;
        Y = sixDOF.Y;
        Z = sixDOF.Z;
        rX = sixDOF.rX;
        rY = sixDOF.rY;
        rZ = sixDOF.rZ;
        BuildMatrix();
      }
      /// <summary>
      /// This function converts a 4x4 state matrix into the appropriate 6dof Euler matrix w/ rx,ry,rz rotations.
      /// </summary>
      /// <param name="StateMatrix">array, 4x4 double</param>
      public virtual void Make6DoF(double[,] StateMatrix)
      {
        Electroimpact.LinearAlgebra.cMatrix m = new cMatrix(StateMatrix);
        X = m.X;
        Y = m.Y;
        Z = m.Z;
        rX = m.rXrYrZ[0, 0].DegreesToRadians();
        rY = m.rXrYrZ[1, 0].DegreesToRadians();
        rZ = m.rXrYrZ[2, 0].DegreesToRadians();
      }
      /// <summary>
      /// This function changes the rXrYrZ parameters of the c6dof to the order of rx,ry,rz rotations.
      /// </summary>
      public virtual void Make6DoF_rXrYrZ()
      {
        double[,] temp = this.GetMatrix();

        Electroimpact.LinearAlgebra.cMatrix m = new cMatrix(temp);
        X = m.X;
        Y = m.Y;
        Z = m.Z;
        rX = m.rXrYrZ[0, 0].DegreesToRadians();
        rY = m.rXrYrZ[1, 0].DegreesToRadians();
        rZ = m.rXrYrZ[2, 0].DegreesToRadians();
      }

      /// <summary>
      /// This function converts a 4x4 state matrix into the appropriate 6dof Euler matrix w/ rz,ry,rx rotations.
      /// </summary>
      /// <param name="StateMatrix">array, 4x4 double</param>
      public virtual void Make6DoF_rZrYrX(double[,] StateMatrix)
      {
        Electroimpact.LinearAlgebra.cMatrix m = new cMatrix(StateMatrix);

        X = m.X;
        Y = m.Y;
        Z = m.Z;
        rX = m.rZrYrX[0, 0].DegreesToRadians();
        rY = m.rZrYrX[1, 0].DegreesToRadians();
        rZ = m.rZrYrX[2, 0].DegreesToRadians();
      }
      /// <summary>
      /// pull out values of rZrYrX, in that order, from the existing rotation matrix without modifying it.  .
      /// </summary>
      public virtual double[] Get_rZrYrX()
      {
        double[,] temp = this.GetMatrix();
        Electroimpact.LinearAlgebra.cMatrix m = new cMatrix(temp);
        double [] output = {m.rZrYrX[0, 0].DegreesToRadians(), m.rZrYrX[1, 0].DegreesToRadians(), m.rZrYrX[2, 0].DegreesToRadians()};
        return output;
      }
      
      /// <summary>
      /// Inverts this Euler Transform
      /// </summary>
      public void InvertMe()
      {
        Electroimpact.LinearAlgebra.cMatrix m = new cMatrix(this.GetMatrix());
        m = new cMatrix(m.InverseHomogeneous());
        X = m.X;
        Y = m.Y;
        Z = m.Z;
        rX = m.rXrYrZ[0, 0].DegreesToRadians();
        rY = m.rXrYrZ[1, 0].DegreesToRadians();
        rZ = m.rXrYrZ[2, 0].DegreesToRadians();
      }

      /// <summary>
      /// Returns the inverse of this Transform in the form of a 4x4 state matrix
      /// </summary>
      /// <returns>4x4 array of doubles</returns>
      public double[,] Inverse()
      {
        Electroimpact.LinearAlgebra.cMatrix m = new cMatrix(this.GetMatrix());
        return m.InverseHomogeneous();
      }

      /// <summary>
      /// Return Dot Product of this point
      /// </summary>
      /// <param name="p"></param>
      /// <returns></returns>
      public double[] DotMe(double[] p)
      {
        double[] ret = new double[] { p[0], p[1], p[2], 1 };
        cMatrix m = new cMatrix(this.GetMatrix());
        //p[3] = 1.0;//just in case abuser forgets.
        return m.DotMe(p);
      }
      /// <summary>
      /// returns dot product of this transfrom with the following state matrix
      /// </summary>
      /// <param name="p"></param>
      /// <returns></returns>
      public double[,] DotMe(double[,] p)
      {
        cMatrix m = new cMatrix(this.GetMatrix());
        return m.DotMeLHT(p);
      }

      public double A
      {
        get
        {
          return myMatixClass.A;
        }
      }

      public double B
      {
        get
        {
          return myMatixClass.B;
        }
      }
    }

    /// <summary>
    /// Displace linear, rotate about Y, then X then Z
    /// </summary>
    public class c6dof_rYrXrZ
    {

      protected c6DOF_Args args = new c6DOF_Args();

      //input values in radians please

      #region constructors
      public c6dof_rYrXrZ(c6DOF_Args args_in)
      {
        //rx = new RxMatrix(args_in.rX_radians);
        //ry = new RyMatrix(args_in.rY_radians);
        //rz = new RzMatrix(args_in.rZ_radians);
        ////        rr = new RrMatrix(X, Y, Z);
        //rr = new RrMatrix(args_in.X, args_in.Y, args_in.Z);

        this.args = new c6DOF_Args(args_in);
        BuildMatrix();
      }

      public c6dof_rYrXrZ(double X, double Y, double Z, double rXradians, double rYradians, double rZradians)
      {
        //rx = new RxMatrix(rXradians);
        //ry = new RyMatrix(rYradians);
        //rz = new RzMatrix(rZradians);
        //rr = new RrMatrix(X, Y, Z);

        this.args = new c6DOF_Args(X, Y, Z, rXradians, rYradians, rZradians);
        BuildMatrix();
      }
      /// <summary>
      /// Base constructor
      /// </summary>
      public c6dof_rYrXrZ()
      {
        //rx = new RxMatrix(0);
        //ry = new RyMatrix(0);
        //rz = new RzMatrix(0);
        //rr = new RrMatrix(0, 0, 0);

        this.args = new c6DOF_Args();
        BuildMatrix();
      }
      /// <summary>
      /// Copy Constructor
      /// </summary>
      /// <param name="Input">c6dof class</param>
      public c6dof_rYrXrZ(c6dof_rYrXrZ Input)
      {
        //rx = new RxMatrix(Input.rX);
        //ry = new RyMatrix(Input.rY);
        //rz = new RzMatrix(Input.rZ);
        //rr = new RrMatrix(Input.X, Input.Y, Input.Z);

        this.args = new c6DOF_Args(Input.X, Input.Y, Input.Z, Input.rX, Input.rY, Input.rZ);
        BuildMatrix();
      }
      public c6dof_rYrXrZ(double[,] Input)
      {
        Electroimpact.LinearAlgebra.cMatrix m = new cMatrix(Input);
        //rx = new RxMatrix();
        //ry = new RyMatrix();
        //rz = new RzMatrix();
        //rr = new RrMatrix();
        //rX = m.rXrYrZ[0, 0].DegreesToRadians();
        //rY = m.rXrYrZ[1, 0].DegreesToRadians();
        //rZ = m.rXrYrZ[2, 0].DegreesToRadians();
        //X = m.X;
        //Y = m.Y;
        //Z = m.Z;

        this.args = new c6DOF_Args(m.X, m.Y, m.Z, m.rYrXrZ[0, 0].DegreesToRadians(), m.rYrXrZ[1, 0].DegreesToRadians(), m.rYrXrZ[2, 0].DegreesToRadians());
        BuildMatrix();
      }
      #endregion

      /// <summary>
      /// Builds a matrix based on ryrxry format
      /// </summary>
      /// <returns></returns>
      protected virtual double[,] BuildMatrix()
      {

        double[,] matrix = new double[4, 4];

        double sinrx = Math.Sin(args.rX_radians);
        double sinry = Math.Sin(args.rY_radians);
        double sinrz = Math.Sin(args.rZ_radians);
        double cosrx = Math.Cos(args.rX_radians);
        double cosry = Math.Cos(args.rY_radians);
        double cosrz = Math.Cos(args.rZ_radians);

        matrix[0, 0] = cosry * cosrz + sinrx * sinry * sinrz;
        matrix[1, 0] = cosrx * sinrz;
        matrix[2, 0] = sinrx * cosry * sinrz - sinry * cosrz;
        matrix[3, 0] = 0.0;

        matrix[0, 1] = sinrx * sinry * cosrz - cosry * sinrz;
        matrix[1, 1] = cosrx * cosrz;
        matrix[2, 1] = sinry * sinrz + sinrx * cosry * cosrz;
        matrix[3, 1] = 0.0;

        matrix[0, 2] = cosrx * sinry;
        matrix[1, 2] = -sinrx;
        matrix[2, 2] = cosrx * cosry;
        matrix[3, 2] = 0.0;

        matrix[0, 3] = args.X;
        matrix[1, 3] = args.Y;
        matrix[2, 3] = args.Z;
        matrix[3, 3] = 1.0;

        return matrix;
      }
      /// <summary>
      /// Gets the matrix based on the arguments below.  Pass in radians please.
      /// </summary>
      /// <param name="X">linear distance</param>
      /// <param name="Y">linear distance</param>
      /// <param name="Z">linear distance</param>
      /// <param name="rX">radians</param>
      /// <param name="rY">radians</param>
      /// <param name="rZ">radians</param>
      /// <returns></returns>
      public double[,] GetMatrix(double X, double Y, double Z, double rXradians, double rYradians, double rZradians)
      {
        args = new c6DOF_Args(X, Y, Z, rXradians, rYradians, rZradians);
        return BuildMatrix();
      }
      /// <summary>
      /// Returns a state matrix represending this transform
      /// </summary>
      /// <returns></returns>
      public double[,] GetMatrix()
      {
        return BuildMatrix();
      }

      public cMatrix myMatixClass
      {
        get
        {
          return new cMatrix(GetMatrix());
        }
      }

      /// <summary>
      /// Returns only the roatation portion of the state matrix representing this transform
      /// </summary>
      /// <returns></returns>
      public double[,] GetOrientationMatrix()
      {
        double[,] ret = BuildMatrix();
        ret[0, 3] = 0.0;
        ret[1, 3] = 0.0;
        ret[2, 3] = 0.0;
        return ret;
      }
      /// <summary>
      /// a double value
      /// </summary>
      public double X
      {
        get { return args.X; }
        set { args.X = value; }
      }
      /// <summary>
      /// a double value
      /// </summary>
      public double Y
      {
        get { return args.Y; }
        set { args.Y = value; }
      }
      /// <summary>
      /// a double value
      /// </summary>
      public double Z
      {
        get { return args.Z; }
        set { args.Z = value; }
      }
      /// <summary>
      /// Sets or returns a value in radians
      /// </summary>
      public double rX
      {
        get { return args.rX_radians; }
        set { args.rX_radians = value; }
      }
      /// <summary>
      /// Sets or returns a value in radians
      /// </summary>      
      public double rY
      {
        get { return args.rY_radians; }
        set { args.rY_radians = value; }
      }
      /// <summary>
      /// Sets or returns a value in radians
      /// </summary>      
      public double rZ
      {
        get { return args.rZ_radians; }
        set { args.rZ_radians = value; }
      }

      /// <summary>
      /// Copy function.  I haven't overrided =, == etc.  Too much work..use this.
      /// </summary>
      /// <param name="sixDOF"></param>
      public void Copy(c6dof_rYrXrZ sixDOF)
      {
        X = sixDOF.X;
        Y = sixDOF.Y;
        Z = sixDOF.Z;
        rX = sixDOF.rX;
        rY = sixDOF.rY;
        rZ = sixDOF.rZ;
        BuildMatrix();
      }
      /// <summary>
      /// This function converts a 4x4 state matrix into the appropriate 6dof Euler matrix w/ rx,ry,rz rotations.
      /// </summary>
      /// <param name="StateMatrix">array, 4x4 double</param>
      public virtual void Make6DoF(double[,] StateMatrix)
      {
        Electroimpact.LinearAlgebra.cMatrix m = new cMatrix(StateMatrix);
        X = m.X;
        Y = m.Y;
        Z = m.Z;
        rX = m.rYrXrZ[0, 0].DegreesToRadians();
        rY = m.rYrXrZ[1, 0].DegreesToRadians();
        rZ = m.rYrXrZ[2, 0].DegreesToRadians();
      }
      /// <summary>
      /// This function changes the rXrYrZ parameters of the c6dof to the order of rx,ry,rz rotations.
      /// </summary>
      //public virtual void Make6DoF_rXrYrZ()
      //{
      //  double[,] temp = this.GetMatrix();

      //  Electroimpact.LinearAlgebra.cMatrix m = new cMatrix(temp);
      //  X = m.X;
      //  Y = m.Y;
      //  Z = m.Z;
      //  rX = m.rXrYrZ[0, 0].DegreesToRadians();
      //  rY = m.rXrYrZ[1, 0].DegreesToRadians();
      //  rZ = m.rXrYrZ[2, 0].DegreesToRadians();
      //}

      /// <summary>
      /// This function converts a 4x4 state matrix into the appropriate 6dof Euler matrix w/ rz,ry,rx rotations.
      /// </summary>
      /// <param name="StateMatrix">array, 4x4 double</param>
      public virtual void Make6DoF_rZrYrX(double[,] StateMatrix)
      {
        Electroimpact.LinearAlgebra.cMatrix m = new cMatrix(StateMatrix);

        X = m.X;
        Y = m.Y;
        Z = m.Z;
        rX = m.rZrYrX[0, 0].DegreesToRadians();
        rY = m.rZrYrX[1, 0].DegreesToRadians();
        rZ = m.rZrYrX[2, 0].DegreesToRadians();
      }
      /// <summary>
      /// pull out values of rZrYrX, in that order, from the existing rotation matrix without modifying it.  .
      /// </summary>
      public virtual double[] Get_rZrYrX()
      {
        double[,] temp = this.GetMatrix();
        Electroimpact.LinearAlgebra.cMatrix m = new cMatrix(temp);
        double[] output = { m.rZrYrX[0, 0].DegreesToRadians(), m.rZrYrX[1, 0].DegreesToRadians(), m.rZrYrX[2, 0].DegreesToRadians() };
        return output;
      }

      /// <summary>
      /// Inverts this Euler Transform
      /// </summary>
      public void InvertMe()
      {
        Electroimpact.LinearAlgebra.cMatrix m = new cMatrix(this.GetMatrix());
        m = new cMatrix(m.InverseHomogeneous());
        X = m.X;
        Y = m.Y;
        Z = m.Z;
        rX = m.rYrXrZ[0, 0].DegreesToRadians();
        rY = m.rYrXrZ[1, 0].DegreesToRadians();
        rZ = m.rYrXrZ[2, 0].DegreesToRadians();
      }

      /// <summary>
      /// Returns the inverse of this Transform in the form of a 4x4 state matrix
      /// </summary>
      /// <returns>4x4 array of doubles</returns>
      public double[,] Inverse()
      {
        Electroimpact.LinearAlgebra.cMatrix m = new cMatrix(this.GetMatrix());
        return m.InverseHomogeneous();
      }

      /// <summary>
      /// Return Dot Product of this point
      /// </summary>
      /// <param name="p"></param>
      /// <returns></returns>
      public double[] DotMe(double[] p)
      {
        double[] ret = new double[] { p[0], p[1], p[2], 1 };
        cMatrix m = new cMatrix(this.GetMatrix());
        //p[3] = 1.0;//just in case abuser forgets.
        return m.DotMe(p);
      }
      /// <summary>
      /// returns dot product of this transfrom with the following state matrix
      /// </summary>
      /// <param name="p"></param>
      /// <returns></returns>
      public double[,] DotMe(double[,] p)
      {
        cMatrix m = new cMatrix(this.GetMatrix());
        return m.DotMeLHT(p);
      }
    }

    /// <summary>
    /// rotate about X, Y, then Z then translate.
    /// </summary>
    public class c6dofFixedAngle
    {
      #region Members

      private RxMatrix rx;
      private RyMatrix ry;
      private RzMatrix rz;
      private RrMatrix rr;
      //input values in radians please
      public c6dofFixedAngle(double X, double Y, double Z, double rX, double rY, double rZ)
      {
        rx = new RxMatrix(rX);
        ry = new RyMatrix(rY);
        rz = new RzMatrix(rZ);
        rr = new RrMatrix(X, Y, Z);
      }
      public c6dofFixedAngle()
      {
        rx = new RxMatrix(0);
        ry = new RyMatrix(0);
        rz = new RzMatrix(0);
        rr = new RrMatrix(0, 0, 0);
      }
      public c6dofFixedAngle(double[,] Input)
      {
        double d2r = Math.PI / 180.0;
        Electroimpact.LinearAlgebra.cMatrix m = new cMatrix(Input);
        rx = new RxMatrix();
        ry = new RyMatrix();
        rz = new RzMatrix();
        rr = new RrMatrix();
        rX = m.rXrYrZ[0, 0] * d2r;
        rY = m.rXrYrZ[1, 0] * d2r;
        rZ = m.rXrYrZ[2, 0] * d2r;



        X = m.X;
        Y = m.Y;
        Z = m.Z;
      }
      #endregion

      private double[,] BuildMatrix()
      {
        Transforms t = new Transforms();

        t.AddTransform(rr.GetMatrix);
        t.AddTransform(rz.GetMatrix);
        t.AddTransform(ry.GetMatrix);
        t.AddTransform(rx.GetMatrix);
        return t.DotProduct();
      }
      /// <summary>
      /// Gets the matrix based on the arguments below.  Pass in radians please.
      /// </summary>
      /// <param name="X">linear distance</param>
      /// <param name="Y">linear distance</param>
      /// <param name="Z">linear distance</param>
      /// <param name="rX">radians</param>
      /// <param name="rY">radians</param>
      /// <param name="rZ">radians</param>
      /// <returns></returns>
      public double[,] GetMatrix(double X, double Y, double Z, double rX, double rY, double rZ)
      {
        rx.Phi = rX;
        ry.Phi = rY;
        rz.Phi = rZ;
        rr.X = X;
        rr.Y = Y;
        rr.Z = Z;
        return BuildMatrix();
      }
      public double[,] GetMatrix()
      {
        return BuildMatrix();
      }

      public double[,] GetOrientationMatrix()
      {
        double[,] ret = BuildMatrix();
        ret[0, 3] = 0.0;
        ret[1, 3] = 0.0;
        ret[2, 3] = 0.0;
        return ret;
      }
      public cMatrix myMatixClass
      {
        get
        {
          return new cMatrix(GetMatrix());
        }
      }
      public double X
      {
        get { return rr.X; }
        set { rr.X = value; }
      }
      public double Y
      {
        get { return rr.Y; }
        set { rr.Y = value; }
      }
      public double Z
      {
        get { return rr.Z; }
        set { rr.Z = value; }
      }
      public double rX
      {
        get { return rx.Phi; }
        set { rx.Phi = value; }
      }
      public double rY
      {
        get { return ry.Phi; }
        set { ry.Phi = value; }
      }
      public double rZ
      {
        get { return rz.Phi; }
        set { rz.Phi = value; }
      }
    }

    #endregion

  }
  namespace Machine
  {
    #region Axis Class
    /// <summary>
    /// An axis of a CNC or Robotic Machine and all the stuff I like to use to make it work right.
    /// </summary>
    public class cAxis
    {
      #region Members and Constuctor
      private string name;
      private bool IsRotary;
      private bool IsHidden = false;

      private double scalefactor = 1.0;
      private double offset = 0.0;
      private double myPos = 0;
      private const double DegToRad = Math.PI / 180.0;
      private bool using_Degrees = false;
      private bool rollover = false;

      /// <summary>
      /// Constructor requires a Name and to know if this axis is Rotary.
      /// </summary>
      /// <param name="Name">string Unique Axis name must match setup file</param>
      /// <param name="Rotary">bool</param>
      public cAxis(string Name, bool Rotary, bool UseDegrees)
      {
        name = Name;
        this.IsRotary = Rotary;
        this.using_Degrees = UseDegrees;
      }
      public cAxis(string Name, bool Rotary)
      {
        name = Name;
        this.IsRotary = Rotary;
        //this.using_Degrees = UseDegrees;
      }
      public cAxis(string Name, bool Rotary, bool UseDegrees, bool Rollover)
      {
        name = Name;
        this.IsRotary = Rotary;
        this.using_Degrees = UseDegrees;
        this.rollover = Rollover;
      }


      #endregion

      #region Properties

      /// <summary>
      /// Defined in constructor, this is read only afterwards.
      /// </summary>
      public string Name
      {
        get { return this.name; }
      }
      /// <summary>
      /// double value.  The axis is position will be multiplied by this number when used in calculations.
      /// <example>ScaleFactor = 1.006 <br/> 
      /// base.value returns (this.myPos - Offset) * ScaleFactor</example>
      /// </summary>
      public double ScaleFactor
      {
        get { return this.scalefactor; }
        set { this.scalefactor = value; }
      }
      /// <summary>
      /// double value.  The Offset will be subtracted from the axis position before multiplied by ScaleFactor when used in calculations.
      /// <example>Offset = 100.0 <br/> 
      /// base.value returns (this.myPos + Offset) * ScaleFactor</example>
      /// </summary>
      public double Offset
      {
        get { return this.offset; }
        set { this.offset = value; }
      }

      public bool Rollover
      {
        get { return this.rollover; }
        set { this.rollover = value; }
      }

      /// <summary>
      /// 
      /// </summary>
      public bool Hidden
      {
        get { return this.IsHidden; }
        set { this.IsHidden = value; }
      }
      /// <summary>
      /// If this is a rotary axis, the value must be in degrees!  This function returns a value in degrees.
      /// Use the property "Value" if you want scaled, offsetted position in radians for caliculations.
      /// </summary>
      public double Position
      {
        get
        {
          double pos = this.IsRotary && !using_Degrees ? myPos / DegToRad : myPos;
          return pos;
        }
        set
        {
          this.myPos = IsRotary && !using_Degrees ? value * DegToRad : value;
        }
      }
      /// <summary>
      /// Used in calculations not for display!!!!
      /// <example>myArg = Value where value is base.myPos * base.ScaleFactor - base.Offset</example>
      /// </summary>
      public double Value
      {
        get
        {
          return (this.myPos - Offset) * ScaleFactor + this.myPos;
        }
      }
      /// <summary>
      /// true causes this axis to behave like a rotary axis.  This parameter can only be set in the constructor.
      /// </summary>
      public bool Rotary
      {
        get { return this.IsRotary; }
      }
      /// <summary>
      /// Gets or sets the amount of compensation required for this axis.
      /// This is not added into the value function right now...Possibly later.
      /// </summary>
      #endregion
    }
    #endregion

    #region Basic Machine Class
    /// <summary>
    ///  An Electroimpact machine. 
    /// </summary>
    public class cMachine //: IMachine
    {
      #region Members and Constructor
      private const double SmallestIncrement = 1.0e-25;
      protected Electroimpact.StringCalc.cStringCalc myCalculator = new Electroimpact.StringCalc.cStringCalc();
      protected System.Collections.Hashtable myAxisHash = new System.Collections.Hashtable();
      protected System.Collections.ArrayList myDHs = new System.Collections.ArrayList();
      protected Electroimpact.LinearAlgebra.Transforms myXfroms = new Electroimpact.LinearAlgebra.Transforms();
      protected cCompTable myCompTable = new cCompTable();
      protected bool HookedUp = false;
      protected string FileName;
      protected bool Rebuild = true;
      protected bool comptableoverride = false;
      protected bool HasBaseXform = false;
      protected string CompTableAxis = "none";
      protected int station_count = 0;

      #region Variable to DH class
      //The following is a hashtable of arguments used.  If the argument is used in a DH, it will be assigned a number
      //corresponding to this DH.  If the argument is written to, the DH will be re calculated...otherwise it will 
      //be left alone.  Tough problem actually.  Shit.  The string may contain Xm and lXmx.  Xm would be in appropriately
      //assigned to a DH with lXmx.  Sucks.  Needs more planning.

      private class cVariableLocation
      {
        private ArrayList myVariables = new ArrayList();
        private Hashtable myAssignedVariables = new Hashtable();
        System.Text.RegularExpressions.Regex objAlphaPattern = new System.Text.RegularExpressions.Regex("[^a-zA-Z]");

        public enum eArgument
        {
          Theta,
          Alpha,
          A,
          D, //Above for DH arguments, Below for 6dof arguments
          X,
          Y,
          Z,
          rX,
          rY,
          rZ
        }

        private class cDHVariable
        {
          public int DHIndex;
          public eArgument myType = new eArgument();

          public cDHVariable(int DHidx, eArgument ArgType)
          {
            DHIndex = DHidx;
            myType = ArgType;
          }
        }

        /// <summary>
        /// Once the entire machine config file is read, hand this class the VariableNames and the DHs derived from the file.
        /// This function will assignblame to the appropriate DHs.  Later, when changing a variables value, call MarkVairable to 
        /// cause the DH to update the appropriate arguments and matrix.
        /// </summary>
        /// <param name="DHs"></param>
        public void AssignVariablesToDHs(ArrayList DHs, bool Debug)
        {
          int ii;
          for (ii = 0; ii < DHs.Count; ii++)
          {
            object o = DHs[ii];
            if (o is Electroimpact.LinearAlgebra.DHMatrix)
            {
              LinearAlgebra.DHMatrix dh = (LinearAlgebra.DHMatrix)DHs[ii];
              FindVariable(dh.myStrings.Theta, ii, eArgument.Theta, Debug);
              FindVariable(dh.myStrings.Alpha, ii, eArgument.Alpha, Debug);
              FindVariable(dh.myStrings.a, ii, eArgument.A, Debug);
              FindVariable(dh.myStrings.d, ii, eArgument.D, Debug);
            }
            else if (o is Electroimpact.LinearAlgebra.string6DOF)
            {
              LinearAlgebra.string6DOF dh = (LinearAlgebra.string6DOF)DHs[ii];
              FindVariable(dh.myStrings.X, ii, eArgument.X, Debug);
              FindVariable(dh.myStrings.Y, ii, eArgument.Y, Debug);
              FindVariable(dh.myStrings.Z, ii, eArgument.Z, Debug);
              FindVariable(dh.myStrings.rX, ii, eArgument.rX, Debug);
              FindVariable(dh.myStrings.rY, ii, eArgument.rY, Debug);
              FindVariable(dh.myStrings.rZ, ii, eArgument.rZ, Debug);
            }
          }
          string message = "";
          for (ii = 0; ii < myVariables.Count; ii++)
          {
            if (!myAssignedVariables.ContainsKey((string)myVariables[ii]))
              message += (string)myVariables[ii] + "\n";
          }
          if (message.Length > 0 && Debug)
          {
            System.Windows.Forms.MessageBox.Show("Not all variables were assigned to at least one DH.  Here is a list: \n" + message);
          }
        }

        private void FindVariable(string Argument, int DHIndex, eArgument ArgType, bool Debug)
        {
          for (int ii = 0; ii < myVariables.Count; ii++)
          {
            int idx = 0;
            string avar = (string)myVariables[ii];
            while (idx != -1)
            {
              idx = Argument.IndexOf(avar, idx);
              if (idx >= 0) //Try to find the bastard
              {
                bool checkL = false;
                bool checkR = false;

                //Check left, if at beginning of string, left is good.
                if (idx == 0)
                  checkL = true;
                else
                {
                  string l = Argument.Substring(idx - 1, 1);
                  checkL = objAlphaPattern.IsMatch(l);
                }
                if (checkL) //No reason to check Right if Left fails.
                {
                  if (idx + avar.Length == Argument.Length)
                    checkR = true;
                  else
                  {
                    string r = Argument.Substring(idx + avar.Length, 1);
                    checkR = objAlphaPattern.IsMatch(r);
                  }

                  if (checkL && checkR)
                  {
                    cDHVariable cv = new cDHVariable(DHIndex, ArgType);
                    if (myAssignedVariables.ContainsKey(avar))
                    {
                      if (Debug)
                        System.Windows.Forms.MessageBox.Show("Warning!  You have assigned " + avar + " to more than one DH.  This is OK, to kill this warning in the future, set the variable _Debug to false in your config file.");
                      object o = myAssignedVariables[avar];
                      if (o is ArrayList)
                        ((ArrayList)o).Add(cv); //this should work, however, you need to test to ensure it does.  Tested...it works!
                      else if (o is cDHVariable)
                      {
                        ArrayList al = new ArrayList();
                        al.Add(cv);
                        al.Add((cDHVariable)o);
                        myAssignedVariables.Remove(avar);
                        myAssignedVariables.Add(avar, al);
                      }
                    }
                    else
                    {
                      myAssignedVariables.Add(avar, cv);
                    }
                  }
                }
              }
              if (idx >= 0)
                idx += ((string)myVariables[ii]).Length;
            }
          }
        }

        /// <summary>
        /// Markes the appropriate DH as needs recalc
        /// </summary>
        /// <param name="Variable">A variable name</param>
        /// <param name="DHs">Array list of DH matrix classes for the machine you are working on.</param>
        public void MarkVariable(string Variable, ArrayList DHs)
        {
          if (myAssignedVariables.ContainsKey(Variable))
          {
            object o = myAssignedVariables[Variable];
            if (o is ArrayList) //happens when a variable shows up in more than one dh.
            {
              foreach (cDHVariable dhv in (ArrayList)o)
              {
                if (dhv.DHIndex < DHs.Count)
                {
                  if (DHs[dhv.DHIndex] is LinearAlgebra.DHMatrix)
                  {
                    Electroimpact.LinearAlgebra.DHMatrix dh = (Electroimpact.LinearAlgebra.DHMatrix)DHs[dhv.DHIndex];
                    dh.RecalcReqd = true;
                    switch (dhv.myType)
                    {
                      case eArgument.Theta:
                        dh.RecalcTheta = true;
                        break;
                      case eArgument.Alpha:
                        dh.RecalcAlpha = true;
                        break;
                      case eArgument.A:
                        dh.RecalcA = true;
                        break;
                      case eArgument.D:
                        dh.RecalcD = true;
                        break;
                    }
                  }
                  else if (DHs[dhv.DHIndex] is LinearAlgebra.string6DOF)
                  {
                    Electroimpact.LinearAlgebra.string6DOF dh = (Electroimpact.LinearAlgebra.string6DOF)DHs[dhv.DHIndex];
                    dh.RecalcReqd = true;
                    switch (dhv.myType)
                    {
                      case eArgument.X:
                        dh.RecalcX = true;
                        break;
                      case eArgument.Y:
                        dh.RecalcY = true;
                        break;
                      case eArgument.Z:
                        dh.RecalcZ = true;
                        break;
                      case eArgument.rX:
                        dh.RecalcrX = true;
                        break;
                      case eArgument.rY:
                        dh.RecalcrY = true;
                        break;
                      case eArgument.rZ:
                        dh.RecalcrZ = true;
                        break;
                    }
                  }
                }
              }
            }
            else
            {
              cDHVariable dhv = (cDHVariable)myAssignedVariables[Variable];
              if (dhv.DHIndex < DHs.Count)
              {
                if (DHs[dhv.DHIndex] is LinearAlgebra.DHMatrix)
                {
                  Electroimpact.LinearAlgebra.DHMatrix dh = (Electroimpact.LinearAlgebra.DHMatrix)DHs[dhv.DHIndex];
                  dh.RecalcReqd = true;
                  switch (dhv.myType)
                  {
                    case eArgument.Theta:
                      dh.RecalcTheta = true;
                      break;
                    case eArgument.Alpha:
                      dh.RecalcAlpha = true;
                      break;
                    case eArgument.A:
                      dh.RecalcA = true;
                      break;
                    case eArgument.D:
                      dh.RecalcD = true;
                      break;
                  }
                }
                else if (DHs[dhv.DHIndex] is LinearAlgebra.string6DOF)
                {
                  Electroimpact.LinearAlgebra.string6DOF dh = (Electroimpact.LinearAlgebra.string6DOF)DHs[dhv.DHIndex]; dh.RecalcReqd = true;
                  switch (dhv.myType)
                  {
                    case eArgument.X:
                      dh.RecalcX = true;
                      break;
                    case eArgument.Y:
                      dh.RecalcY = true;
                      break;
                    case eArgument.Z:
                      dh.RecalcZ = true;
                      break;
                    case eArgument.rX:
                      dh.RecalcrX = true;
                      break;
                    case eArgument.rY:
                      dh.RecalcrY = true;
                      break;
                    case eArgument.rZ:
                      dh.RecalcrZ = true;
                      break;
                  }
                }
              }
            }
          }
        }

        public void Add(string VariableName)
        {
          myVariables.Add(VariableName);
        }

        internal void Clear()
        {
          myVariables.Clear();
          myAssignedVariables.Clear();
        }
      }
      #endregion

      cVariableLocation myVariableToDH_Hash = new cVariableLocation();

      //The following is used for writing out data later on.
      private System.Collections.ArrayList myStaticVarNames = new ArrayList();
      //private System.Collections.ArrayList myOtherStaticVarNames = new ArrayList();
      private List<string> myOtherStaticVarNames = new List<string>();
      private System.Collections.ArrayList myMachineAxisNames = new ArrayList();

      /// <summary>
      /// This class requires a setup file.  Please see example located somewhere.
      /// </summary>
      /// <param name="FileName">string representing a loc and filename e.g. c:\doggy\doggystyle.config</param>
      public cMachine(string FileName, bool UseDegrees)
      {
        this.FileName = FileName;
        this.HookedUp = this.HookMeUp(FileName);
        this.myCalculator.Degrees = UseDegrees;
        //this.myCalculator.Degrees = false;
      }

      protected cMachine()
      {
      }
      #endregion

      #region Setup Functions

      #region Reading In File
      protected bool HookMeUp(string FileName)
      {
        XmlNodeReader r;
        XmlDocument document;
        bool rtn = true;
        try
        {
          document = new XmlDocument();
          document.Load(FileName);
          r = new XmlNodeReader(document);
        }
        catch (Exception ex)
        {
          System.Windows.Forms.MessageBox.Show("Config File: " + FileName.ToString() + " is missing or invalid.");
          throw (ex);
        }

        while (r.Read())
        {
          switch (r.NodeType)
          {
            case XmlNodeType.Element:
              {
                if (r.Name == "EITransform")
                {
                  while (r.Read() && rtn == true)
                  {
                    if (r.NodeType == XmlNodeType.Element)
                    {
                      if (r.Name == "variables")
                      {
                        rtn = ReadInVariables(r);
                      }
                      if (r.Name == "TransForms")
                      {
                        rtn = ReadInTransforms(r);
                      }
                    }
                  }
                }
                break;
              }
          }
        }
        r.Close();
        r = null;
        document = null;
        double debug = myCalculator.SimpleCalc("_Debug");
        bool _Debug = double.IsNaN(debug) || System.Math.Round(debug, 0) == 1;
        myVariableToDH_Hash.AssignVariablesToDHs(myDHs, _Debug);

        return rtn;
      }

      private bool ReadInVariables(XmlNodeReader r)
      {
        myVariableToDH_Hash.Clear();
        myCompTable.Clear();
        myMachineAxisNames.Clear();
        this.myStaticVarNames.Clear();

        return ReadInVariables(r, false);

      }

      private bool ReadInVariables(XmlNodeReader r, bool PartialConfigurationFile)
      {
        //Only clear variables if we are loading a whole new machine config file.
        if (!PartialConfigurationFile)
        {
          myVariableToDH_Hash.Clear();
          myCompTable.Clear();
          myMachineAxisNames.Clear();
          this.myStaticVarNames.Clear();
        }
        while (r.Read())
        {
          if (r.NodeType == XmlNodeType.EndElement)
          {
            if (r.Name != "variables")
            {
              Exception ex = new Exception("XML file format error!");
              throw (ex);
            }
            for (int ii = 0; ii < myMachineAxisNames.Count; ii++)
            {
              if (myCompTable.ContainsAxis((string)myMachineAxisNames[ii]))
              {
                System.Collections.Hashtable args = myCompTable.GetCompVariables((string)myMachineAxisNames[ii], ((cAxis)myAxisHash[(string)myMachineAxisNames[ii]]).Position);
                string[] names = new string[args.Count];
                args.Keys.CopyTo(names, 0);
                for (int jj = 0; jj < names.Length; jj++)
                {
                  myCalculator._AssignVariable(names[jj], (double)args[names[jj]]);
                }
              }
            }

            for (int ii = 0; ii < myMachineAxisNames.Count; ii++)
            {
              cAxis cax = (cAxis)myAxisHash[myMachineAxisNames[ii].ToString()];
              if (cax.Rotary && cax.Rollover)
              {
                if (myCompTable.ContainsAxis(myMachineAxisNames[ii].ToString()))
                {
                  List<cCompTable.cCompStation> ct = myCompTable.GetCompensationTable(myMachineAxisNames[ii].ToString());
                  bool GotTheZero = false;
                  bool GotThe360 = false;
                  for (int jj = 0; jj < ct.Count; jj++)
                  {
                    cCompTable.cCompStation cs = ct[jj];
                    if (Math.Abs(cs.loc) < 1.0e-6)
                      GotTheZero = true;
                    if (Math.Abs(360.0 - cs.loc) < 1.0e-6)
                      GotThe360 = true;
                  }
                  ct = myCompTable.GetCompensationTable(myMachineAxisNames[ii].ToString());
                  if (!GotTheZero)
                  {
                    cCompTable.cCompStation cs = new cCompTable.cCompStation();
                    List<cCompTable.cCompStation> newTable = new List<cCompTable.cCompStation>();
                    cs.loc = 0.0;
                    cs.X = 0.0;
                    cs.Y = 0.0;
                    cs.Z = 0.0;
                    cs.rX = 0.0;
                    cs.rY = 0.0;
                    cs.rZ = 0.0;
                    bool addedit = false;
                    for (int jj = 0; jj < ct.Count; jj++)
                    {
                      if (ct[jj].loc > 0.0 && !addedit)
                      {
                        newTable.Add(cs);
                        addedit = true;
                      }
                      newTable.Add(ct[jj]);
                    }
                    if (!addedit)
                      newTable.Add(cs);
                    myCompTable.PutCompensationTable(myMachineAxisNames[ii].ToString(), newTable);
                  }
                  ct = myCompTable.GetCompensationTable(myMachineAxisNames[ii].ToString());
                  if (!GotThe360)
                  {
                    cCompTable.cCompStation cs = new cCompTable.cCompStation();
                    List<cCompTable.cCompStation> newTable = new List<cCompTable.cCompStation>();
                    cs.loc = 360.0;
                    cs.X = 0.0;
                    cs.Y = 0.0;
                    cs.Z = 0.0;
                    cs.rX = 0.0;
                    cs.rY = 0.0;
                    cs.rZ = 0.0;
                    bool addedit = false;
                    for (int jj = 0; jj < ct.Count; jj++)
                    {
                      if (ct[jj].loc > 360.0 && !addedit)
                      {
                        newTable.Add(cs);
                        addedit = true;
                      }
                      newTable.Add(ct[jj]);
                    }
                    if (!addedit)
                      newTable.Add(cs);
                    myCompTable.PutCompensationTable(myMachineAxisNames[ii].ToString(), newTable);
                  }
                }
              }
            }
            return true;
          }
          if (r.NodeType == XmlNodeType.Element)
          {
            #region static Variables
            if (r.Name == "static")
            {
              string value = r.GetAttribute("multiply");
              if (value == null)
                value = "10000.0";
              this.myCalculator._AssignVariable("StaticVariableScaler", this.myCalculator.SimpleCalc(value));

              while (r.Read())
              {
                if (r.NodeType == XmlNodeType.EndElement)
                {
                  if (r.Name != "static")
                  {
                    Exception ex = new Exception("XML Format Error!");
                    throw (ex);
                  }
                  break;
                }
                if (r.NodeType == XmlNodeType.Element)
                {//process a kinematic variable.
                  string name = r.Name;

                  if (name == "ToolPointOffset")//we shall slowly irradicate this silly thing.
                  {
                    value = r.GetAttribute("X");
                    //this.myCalculator._AssignVariable("ToolPointOffset_X", this.myCalculator.SimpleCalc(value));
                    this.myVariableToDH_Hash.Add("BaseShift_X");
                    myCalculator._AssignVariable("BaseShift_X", this.myCalculator.SimpleCalc(value));

                    value = r.GetAttribute("Y");
                    //this.myCalculator._AssignVariable("ToolPointOffset_Y", this.myCalculator.SimpleCalc(value));
                    this.myVariableToDH_Hash.Add("BaseShift_Y");
                    this.myCalculator._AssignVariable("BaseShift_Y", this.myCalculator.SimpleCalc(value));

                    value = r.GetAttribute("Z");
                    //this.myCalculator._AssignVariable("ToolPointOffset_Z", this.myCalculator.SimpleCalc(value));
                    this.myVariableToDH_Hash.Add("BaseShift_Z");
                    this.myCalculator._AssignVariable("BaseShift_Z", this.myCalculator.SimpleCalc(value));

                    this.myVariableToDH_Hash.Add("BaseShift_rX");
                    this.myCalculator._AssignVariable("BaseShift_rX", 0.0);

                    this.myVariableToDH_Hash.Add("BaseShift_rY");
                    this.myCalculator._AssignVariable("BaseShift_rY", 0.0);

                    this.myVariableToDH_Hash.Add("BaseShift_rZ");
                    this.myCalculator._AssignVariable("BaseShift_rZ", 0.0);
                  }
                  else if (name == "BaseShift")
                  {
                    value = r.GetAttribute("X");
                    myCalculator._AssignVariable("BaseShift_X", this.myCalculator.SimpleCalc(value));
                    if (!this.myStaticVarNames.Contains("BaseShift_X"))
                    {
                      this.myStaticVarNames.Add("BaseShift_X");
                      this.myVariableToDH_Hash.Add("BaseShift_X");
                    }

                    value = r.GetAttribute("Y");
                    this.myCalculator._AssignVariable("BaseShift_Y", this.myCalculator.SimpleCalc(value));
                    if (!this.myStaticVarNames.Contains("BaseShift_Y"))
                    {
                      this.myStaticVarNames.Add("BaseShift_Y");
                      this.myVariableToDH_Hash.Add("BaseShift_Y");
                    }

                    value = r.GetAttribute("Z");
                    this.myCalculator._AssignVariable("BaseShift_Z", this.myCalculator.SimpleCalc(value));
                    if (!this.myStaticVarNames.Contains("BaseShift_Z"))
                    {
                      this.myStaticVarNames.Add("BaseShift_Z");
                      this.myVariableToDH_Hash.Add("BaseShift_Z");
                    }
                    value = r.GetAttribute("rX");
                    this.myCalculator._AssignVariable("BaseShift_rX", this.myCalculator.SimpleCalc(value));
                    if (!this.myStaticVarNames.Contains("BaseShift_rX"))
                    {
                      this.myStaticVarNames.Add("BaseShift_rX");
                      this.myVariableToDH_Hash.Add("BaseShift_rX");
                    }

                    value = r.GetAttribute("rY");
                    this.myCalculator._AssignVariable("BaseShift_rY", this.myCalculator.SimpleCalc(value));
                    if (!this.myStaticVarNames.Contains("BaseShift_rY"))
                    {
                      this.myStaticVarNames.Add("BaseShift_rY");
                      this.myVariableToDH_Hash.Add("BaseShift_rY");
                    }

                    value = r.GetAttribute("rZ");
                    this.myCalculator._AssignVariable("BaseShift_rZ", this.myCalculator.SimpleCalc(value));
                    if (!this.myStaticVarNames.Contains("BaseShift_rZ"))
                    {
                      this.myStaticVarNames.Add("BaseShift_rZ");
                      this.myVariableToDH_Hash.Add("BaseShift_rZ");
                    }

                    //this.myStaticVarNames.Add(name);
                    //this.myVariableToDH_Hash.Add(name);

                  }
                  else if (name.ToLower() == "_debug")
                  {
                    if (r.GetAttribute("value").ToLower() == "true" || r.GetAttribute("value").ToLower() == "1")
                      this.myCalculator._AssignVariable(name, 1);
                    else
                      this.myCalculator._AssignVariable(name, 0);
                  }
                  else
                  {
                    value = r.GetAttribute("value");
                    bool hide = (r.GetAttribute("Hide") != null) ? r.GetAttribute("Hide").ToLower() == "true" : false;
                    this.myCalculator._AssignVariable(name, this.myCalculator.SimpleCalc(value));
                    string note = r.GetAttribute("note");
                    if (note == null)
                    {
                      note = r.GetAttribute("Note");
                    }
                    if (note != null)
                      this.myCalculator._AssignVariableNote(name, note);

                    if (name.IndexOf("_") != 0 && !hide)
                    {
                      //Make sure we haven't already added this variable.  There is a good chance if we have only a Partial config.
                      if (!this.myStaticVarNames.Contains(name))
                      {
                        this.myStaticVarNames.Add(name);
                        this.myVariableToDH_Hash.Add(name);
                      }
                    }
                    else if (name.IndexOf("_") == 0 || hide)
                    {
                      value = r.GetAttribute("value");
                      myCalculator._AssignVariable(name, this.myCalculator.SimpleCalc(value));
                      if (!this.myOtherStaticVarNames.Contains(name))
                      {
                        this.myOtherStaticVarNames.Add(name);
                        this.myVariableToDH_Hash.Add(name); //090203 nuke comment
                      }
                    }
                  }
                }
              }
            }
            #endregion

            #region Configuring Machine Axes
            if (r.Name == "machineaxes")
            {
              string value = r.GetAttribute("multiply");
              if (value == null)
                value = "1000000.0";
              this.myCalculator._AssignVariable("MachineAxisScaler", this.myCalculator.SimpleCalc(value));

              while (r.Read())
              {
                if (r.NodeType == XmlNodeType.EndElement)
                {
                  if (r.Name != "machineaxes")
                  {
                    Exception ex = new Exception("XML Format Error!");
                    throw (ex);
                  }
                  break;
                }
                if (r.NodeType == XmlNodeType.Element)
                {
                  string name = r.Name;
                  string sScaleFactor = r.GetAttribute("scalefactor");
                  string sOffset = r.GetAttribute("offset");
                  bool bRotary = String.Compare(r.GetAttribute("rotary"), "true", true) == 0;
                  Electroimpact.Machine.cAxis ax = new Electroimpact.Machine.cAxis(name, bRotary, myCalculator.Degrees);
                  double ScaleFactor = Double.TryParse(sScaleFactor, out ScaleFactor) ? ScaleFactor : 0;
                  double Offset = Double.TryParse(sOffset, out Offset) ? Offset : 0.0;
                  bool hidden = (r.GetAttribute("Hide") != null) ? r.GetAttribute("Hide").ToLower() == "true" : false;
                  bool rollover_axis = (r.GetAttribute("rollover") != null) ? r.GetAttribute("rollover").ToLower() == "true" : false;
                  if (ax.Rotary && ax.Rollover)
                    ax.ScaleFactor = 0.0;
                  else
                    ax.ScaleFactor = ScaleFactor;
                  ax.Offset = Offset;
                  ax.Position = 0;
                  ax.Hidden = hidden;
                  ax.Rollover = rollover_axis;
                  this.myCalculator._AssignVariable(name, ax.Value);
                  //We only need to add the axis to these hash's if the thing isn't already there.  Consider Partial Config.
                  if (this.myAxisHash.ContainsKey(name))
                    this.myAxisHash[name] = ax;
                  else
                  {
                    this.myAxisHash.Add(name, ax);
                    this.myMachineAxisNames.Add(name);
                    this.myVariableToDH_Hash.Add(name);
                  }
                }
              }
            }
            #endregion

            #region Comp Tables
            if (r.Name == "comptables")
            {
              string value = r.GetAttribute("linear_multiply");
              if (value == null)
                value = "10000.0";
              this.myCalculator._AssignVariable("LinearCompScaler", this.myCalculator.SimpleCalc(value));
              value = r.GetAttribute("angle_multiply");
              if (value == null)
                value = "10000000000.0";
              this.myCalculator._AssignVariable("AngleCompScaler", this.myCalculator.SimpleCalc(value));

              while (r.Read())
              {
                if (r.NodeType == XmlNodeType.EndElement)
                {
                  if (r.Name != "comptables")
                  {
                    Exception ex = new Exception("XML Format Error!");
                    throw (ex);
                  }
                  break;
                }
                if (r.NodeType == XmlNodeType.Element)
                {
                  if (r.Name == "table")
                  {
                    ArrayList names = myCompTable.AddTable(r);
                    for (int ii = 0; ii < names.Count; ii++)
                      myVariableToDH_Hash.Add((string)names[ii]);
                  }
                }
              }
            }
            #endregion
          }
        }
        return false;
      }
      private bool ReadInTransforms(XmlNodeReader r)
      {
        if (myDHs.Count == 0)
        {
          LinearAlgebra.string6DOF.cArgStrings FromSrings = new Electroimpact.LinearAlgebra.string6DOF.cArgStrings();
          FromSrings.X = "BaseShift_X";
          FromSrings.Y = "BaseShift_Y";
          FromSrings.Z = "BaseShift_Z";
          FromSrings.rX = "BaseShift_rX";
          FromSrings.rY = "BaseShift_rY";
          FromSrings.rZ = "BaseShift_rZ";
          LinearAlgebra.string6DOF sd = new Electroimpact.LinearAlgebra.string6DOF(FromSrings, myCalculator, false);
          myDHs.Add(sd);
        }
        while (r.Read())
        {
          if (r.NodeType == XmlNodeType.EndElement)
          {
            if (r.Name != "TransForms")
            {
              Exception ex = new Exception("XML Format Error!");
              throw (ex);
            }

            for (int ii = 0; ii < myDHs.Count; ii++)
            {
              object dog = new object();
              dog = myDHs[ii];
              //Electroimpact.LinearAlgebra.DHMatrix thisdh = (Electroimpact.LinearAlgebra.DHMatrix)dog;
              //this.myXfroms.AddTransform(thisdh.getDH(this.myCalculator));

              //If you get here, Todd realizes that it would be much better to make a transform class that took care of the switching
              //back and forth.  Good problem for you to solve.  Be carefull though.  The going is tricky.
              if (dog is LinearAlgebra.DHMatrix)
                this.myXfroms.AddTransform(((LinearAlgebra.DHMatrix)dog).getDH(this.myCalculator));
              else if (dog is LinearAlgebra.string6DOF)
                this.myXfroms.AddTransform(((LinearAlgebra.string6DOF)dog).getDH(this.myCalculator));
            }

            Electroimpact.LinearAlgebra.cMatrix m = new Electroimpact.LinearAlgebra.cMatrix(this.myXfroms.DotProduct());
            //Console.WriteLine(m.ToString()); 
            return true;
          }
          if (r.NodeType == XmlNodeType.Element)
          {
            if (r.Name == "DH")
            {
              string id = r.GetAttribute("n");
              string theta = r.GetAttribute("theta");
              string alpha = r.GetAttribute("alpha");
              string a = r.GetAttribute("a");
              string d = r.GetAttribute("d");

              double dtheta = this.myCalculator.SimpleCalc(theta);
              double dalpha = this.myCalculator.SimpleCalc(alpha);
              double da = this.myCalculator.SimpleCalc(a);
              double dd = this.myCalculator.SimpleCalc(d);

              if (double.IsNaN(dtheta))
                ShowDHArgumentError(id, "theta", theta);
              if (double.IsNaN(dalpha))
                ShowDHArgumentError(id, "alpha", alpha);
              if (double.IsNaN(da))
                ShowDHArgumentError(id, "a", a);
              if (double.IsNaN(dd))
                ShowDHArgumentError(id, "d", d);

              Electroimpact.LinearAlgebra.DHMatrix.cDHStrings xf = new Electroimpact.LinearAlgebra.DHMatrix.cDHStrings();
              xf.a = a;
              xf.d = d;
              xf.Theta = theta;
              xf.Alpha = alpha;

              Electroimpact.LinearAlgebra.DHMatrix dh = new Electroimpact.LinearAlgebra.DHMatrix(xf, this.myCalculator, false);
              //dh.myStrings.Theta = theta;
              //dh.myStrings.d = d;
              //dh.myStrings.a = a;
              //dh.myStrings.Alpha = alpha;
              myDHs.Add(dh);
            }
            else if (r.Name == "sixdof")
            {
              string id = r.GetAttribute("n");
              string sX = r.GetAttribute("X");
              string sY = r.GetAttribute("Y");
              string sZ = r.GetAttribute("Z");
              string srX = r.GetAttribute("rX");
              string srY = r.GetAttribute("rY");
              string srZ = r.GetAttribute("rZ");

              double dX = this.myCalculator.SimpleCalc(sX);
              double dY = this.myCalculator.SimpleCalc(sY);
              double dZ = this.myCalculator.SimpleCalc(sZ);
              double drX = this.myCalculator.SimpleCalc(srX);
              double drY = this.myCalculator.SimpleCalc(srY);
              double drZ = this.myCalculator.SimpleCalc(srZ);

              if (double.IsNaN(dX))
                ShowDHArgumentError(id, "X", sX);
              if (double.IsNaN(dY))
                ShowDHArgumentError(id, "Y", sY);
              if (double.IsNaN(dZ))
                ShowDHArgumentError(id, "Z", sZ);
              if (double.IsNaN(drX))
                ShowDHArgumentError(id, "rX", srX);
              if (double.IsNaN(drY))
                ShowDHArgumentError(id, "rY", srY);
              if (double.IsNaN(drZ))
                ShowDHArgumentError(id, "rZ", srZ);

              Electroimpact.LinearAlgebra.string6DOF.cArgStrings xf = new Electroimpact.LinearAlgebra.string6DOF.cArgStrings();
              xf.X = sX;
              xf.Y = sY;
              xf.Z = sZ;
              xf.rX = srX;
              xf.rY = srY;
              xf.rZ = srZ;

              Electroimpact.LinearAlgebra.string6DOF dh = new Electroimpact.LinearAlgebra.string6DOF(xf, this.myCalculator, false);
              myDHs.Add(dh);
            }
          }
        }
        return false;
      }
      private void ShowDHArgumentError(string id, string argumentname, string argumentvalue)
      {
        System.Windows.Forms.MessageBox.Show("Error!  There is a problem with DH number: " + id + " in the " + argumentname + " argument. " + argumentvalue);
      }
      #endregion

      #region Writing Out File



      private void WriteStaticVariables(XmlTextWriter xw)
      {
        xw.WriteStartElement("_Debug");
        xw.WriteAttributeString("value", (Math.Round(this.myCalculator._GetVariable("_Debug")) == 1).ToString());
        xw.WriteEndElement();


        for (int ii = 0; ii < this.myStaticVarNames.Count; ii++)
        {
          if (((string)myStaticVarNames[ii]).StartsWith("BaseShift"))
            continue;
          xw.WriteStartElement(this.myStaticVarNames[ii].ToString());
          xw.WriteAttributeString("value", this.myCalculator._GetVariable((string)this.myStaticVarNames[ii]).ToString("F4"));
          ExtractNote(xw, ii, (string)this.myStaticVarNames[ii]);
          xw.WriteEndElement();
        }
        xw.WriteStartElement("BaseShift");
        {
          string value;
          value = this.myCalculator._GetVariable("BaseShift_X").ToString("F4");
          xw.WriteAttributeString("X", value);
          value = this.myCalculator._GetVariable("BaseShift_Y").ToString("F4");
          xw.WriteAttributeString("Y", value);
          value = this.myCalculator._GetVariable("BaseShift_Z").ToString("F4");
          xw.WriteAttributeString("Z", value);
          value = this.myCalculator._GetVariable("BaseShift_rX").ToString("F4");
          xw.WriteAttributeString("rX", value);
          value = this.myCalculator._GetVariable("BaseShift_rY").ToString("F4");
          xw.WriteAttributeString("rY", value);
          value = this.myCalculator._GetVariable("BaseShift_rZ").ToString("F4");
          xw.WriteAttributeString("rZ", value);
          //}

        }
        xw.WriteEndElement();
        for (int ii = 0; ii < this.myOtherStaticVarNames.Count; ii++)
        {
          xw.WriteStartElement(this.myOtherStaticVarNames[ii].ToString());
          xw.WriteAttributeString("value", this.myCalculator._GetVariable((string)this.myOtherStaticVarNames[ii]).ToString("F4"));
          ExtractNote(xw, ii, this.myOtherStaticVarNames[ii]);
          if (this.myOtherStaticVarNames[ii].ToString().IndexOf("_") != 0)
            xw.WriteAttributeString("Hide", "true");
          xw.WriteEndElement();
        }
      }

      private void ExtractNote(XmlTextWriter xw, int ii, string varname)
      {
        string note = this.myCalculator._GetVariableNote(varname);
        if (note != null)
          xw.WriteAttributeString("note", note);
      }

      private void WriteMachineAxes(XmlTextWriter xw)
      {
        for (int ii = 0; ii < this.myMachineAxisNames.Count; ii++)
        {
          string axname = (string)myMachineAxisNames[ii];
          xw.WriteStartElement(axname);
          xw.WriteAttributeString("rotary", ((cAxis)myAxisHash[axname]).Rotary.ToString());
          xw.WriteAttributeString("scalefactor", ((cAxis)myAxisHash[axname]).ScaleFactor.ToString("F6"));
          xw.WriteAttributeString("offset", ((cAxis)myAxisHash[axname]).Offset.ToString("F6"));
          if (((cAxis)myAxisHash[axname]).Rollover == true)
            xw.WriteAttributeString("rollover", "true");
          if (((cAxis)myAxisHash[axname]).Hidden)
            xw.WriteAttributeString("Hide", "true");
          xw.WriteEndElement();
        }
      }

      #endregion

      protected virtual void UpdateTransforms()
      {
        if (Rebuild)
        {
          //First update the axis values in the calculator.
          string[] keys;
          keys = new string[myAxisHash.Count];
          myAxisHash.Keys.CopyTo(keys, 0);
          for (int ii = 0; ii < keys.Length; ii++)
          {
            myCalculator._AssignVariable(keys[ii], ((cAxis)myAxisHash[keys[ii]]).Value);
            //myVariableToDH_Hash.MarkVariable(keys[ii], myDHs);  Already done in WriteAxisPosition.
            if (myCompTable.ContainsAxis(keys[ii]) && !this.comptableoverride)
            {
              System.Collections.Hashtable h = myCompTable.GetCompVariables(keys[ii], ((cAxis)myAxisHash[keys[ii]]).Position);
              string[] v = new string[h.Count];
              h.Keys.CopyTo(v, 0);
              for (int jj = 0; jj < h.Count; jj++)
              {
                if (Math.Abs(myCalculator._GetVariable(v[jj]) - (double)h[v[jj]]) > SmallestIncrement)
                {
                  myCalculator._AssignVariable(v[jj], (double)h[v[jj]]);
                  myVariableToDH_Hash.MarkVariable(v[jj], myDHs);
                }
              }
            }
          }
          this.myXfroms.Clear();
          for (int ii = 0; ii < myDHs.Count; ii++)
          {
            if (myDHs[ii] is LinearAlgebra.DHMatrix)
              this.myXfroms.AddTransform(((Electroimpact.LinearAlgebra.DHMatrix)myDHs[ii]).getDH(myCalculator));
            else if (myDHs[ii] is LinearAlgebra.string6DOF)
              this.myXfroms.AddTransform(((Electroimpact.LinearAlgebra.string6DOF)myDHs[ii]).getDH(myCalculator));
          }
          Rebuild = false;
          DotProduct();
        }
      }

      #endregion

      #region Public Properties

      /// <summary>
      /// Indicates if constructor appeared successful.
      /// </summary>
      public bool IsHookedUp
      {
        get { return this.HookedUp; }
      }

      public double[,] my4x4
      {
        get
        {
          UpdateTransforms();
          return this.myXfroms.myStateMatrix.GetMatrix;
        }
      }

      public double X
      {
        get
        {
          UpdateTransforms();
          //double xoffset = myCalculator._GetVariable("ToolPointOffset_X");
          return this.myXfroms.X;// +xoffset;
        }
      }

      public double Y
      {
        get
        {
          UpdateTransforms();
          return this.myXfroms.Y;// +myCalculator._GetVariable("ToolPointOffset_Y");
        }
      }

      public double Z
      {
        get
        {
          UpdateTransforms();
          return this.myXfroms.Z;// +myCalculator._GetVariable("ToolPointOffset_Z");
        }
      }

      /// <summary>
      /// two dimensional array of the A pos.
      /// 0 - returns 0 to 360.
      /// 1 - returns +/- 180.
      /// </summary>
      public double A
      {
        get
        {
          UpdateTransforms();
          return this.myXfroms.A;
        }
      }

      public double[] rXrYrZFixedAngle
      {
        get
        {
          UpdateTransforms();
          return this.myXfroms.rXrYrZFixedAngle;
        }
      }
      /// <summary>
      /// rXrYrz returns a 3 x 2 matrix.  Column one is one solution, Column two is the other.
      /// Utilize your axis limits to determine the correct solution for this case.
      /// this number is in degrees!  annoying!
      /// </summary>
      public double[,] rXrYrZ
      {
        get
        {
          UpdateTransforms();
          return this.myXfroms.rXrYrZ;
        }
      }

      public double[,] rYrXrZ
      {
        get
        {
          UpdateTransforms();
          return this.myXfroms.rYrXrZ;
        }
      }

      /// <summary>
      /// Standard robot notation.  Toolpoint rotated about Z then by rotated Y then rotated X
      /// </summary>
      /// 
      public double[,] rZrYrX
      {
        get
        {
          UpdateTransforms();
          return this.myXfroms.rZrYrX;
        }
      }


      /// <summary>
      /// two dimensional array of the B pos.
      /// 0 - returns 0 to 360.
      /// 1 - returns +/- 180.
      /// </summary>
      public double B
      {
        get
        {
          UpdateTransforms();
          return this.myXfroms.B;
        }
      }

      /// <summary>
      /// two dimensional array of the C pos.
      /// 0 - returns 0 to 360.
      /// 1 - returns +/- 180.
      /// </summary>
      public double[] C
      {
        get
        {
          UpdateTransforms();
          return this.myXfroms.C;
        }
      }
      /// <summary>

      public string[] DHStrings
      {

        get
        {
          int count = myDHs.Count;
          string[] dhstrings = new string[count];
          int jj = 0;
          for (int ii = 0; ii < count; ii++)
          {
            if (myDHs[ii] is LinearAlgebra.DHMatrix)
            {
              Electroimpact.LinearAlgebra.DHMatrix.cDHStrings s = ((Electroimpact.LinearAlgebra.DHMatrix)myDHs[ii]).myStrings;
              dhstrings[jj] = "Theta " + s.Theta + " d " + s.d + " a " + s.a + " alpha " + s.Alpha; jj++;
            }
          }
          return dhstrings;
        }
      }

      #endregion

      #region Public Methods
      public bool AddConfigFile(string FileName)
      {
        XmlNodeReader r;
        XmlDocument document;
        bool rtn = true;
        try
        {
          document = new XmlDocument();
          document.Load(FileName);
          r = new XmlNodeReader(document);
        }
        catch (Exception ex)
        {
          throw (ex);
        }

        while (r.Read())
        {
          switch (r.NodeType)
          {
            case XmlNodeType.Element:
              {
                if (r.Name == "EITransform")
                {
                  while (r.Read() && rtn == true)
                  {
                    if (r.NodeType == XmlNodeType.Element)
                    {
                      if (r.Name == "variables")
                      {
                        rtn = ReadInVariables(r, true);
                      }
                    }
                  }
                }
                break;
              }
          }
        }
        r.Close();
        r = null;
        document = null;

        foreach (object o in myDHs)
        {
          if (o is LinearAlgebra.string6DOF)
            ((LinearAlgebra.string6DOF)o).ForceRecalculation();
          if (o is LinearAlgebra.c6dof)
            ((LinearAlgebra.string6DOF)o).ForceRecalculation();
          if (o is Electroimpact.LinearAlgebra.DHMatrix)
            ((LinearAlgebra.DHMatrix)o).ForceRecalculation();
        }
        return rtn;
      }
      public void NullCompTables()
      {
        myCompTable.NullTableValues();
      }
      /// <summary>
      /// Modifies the Position.  This is the raw unscaled machine position of an axis in this machine.  
      /// You must call DotProduct for this to take effect!  You should modify all axis positions before calling DotProduct.
      /// </summary>
      /// <param name="name">string Axis Name (eg Am)</param>
      /// <param name="MachinePos">double position in either degrees or linear unit</param>
      /// <returns></returns>
      public bool WriteAxisPosition(string name, double MachinePos)
      {
        Rebuild = true;
        if (!this.myAxisHash.ContainsKey(name))
          return false;
        Electroimpact.Machine.cAxis ax = (Electroimpact.Machine.cAxis)this.myAxisHash[name];
        if (Math.Abs(ax.Position - MachinePos) > .0000001)
        {
          myVariableToDH_Hash.MarkVariable(name, myDHs);
          ax.Position = MachinePos;
        }
        return true;
      }
      /// <summary>
      /// Intended to be used with the solver tool...otherwise these are set during instantiation from the config file.
      /// </summary>
      /// <param name="axName">one of the axes specified in your config file</param>
      /// <param name="attributeName">one of the attributes allowed: offset, scale1, scale2, scalefactor</param>
      /// <param name="attributeValue">a double value.</param>
      public void WriteAxisAttribute(string axName, string attributeName, double attributeValue)
      {
        Rebuild = true;
        Exception ex;
        if (this.myAxisHash.ContainsKey(axName))
        {
          cAxis ax = (cAxis)this.myAxisHash[axName];
          switch (attributeName)
          {
            case "scalefactor":
              if (Math.Abs(ax.ScaleFactor - attributeValue) > SmallestIncrement)
              {
                if (ax.Rotary && ax.Rollover)
                  ax.ScaleFactor = 0.0;
                else
                  ax.ScaleFactor = attributeValue;
                myVariableToDH_Hash.MarkVariable(axName, myDHs);
              }
              return;
            default:
              ex = new Exception("Invalid attribute name in: " + System.Reflection.MethodInfo.GetCurrentMethod().Name);
              break;
          }
        }
        else
          ex = new Exception("Invalid axis name in: " + System.Reflection.MethodInfo.GetCurrentMethod().Name);
        throw ex;
      }
      /// <summary>
      /// Reads an axis type attribute
      /// </summary>
      /// <param name="axName">string -- the axis name</param>
      /// <param name="attributeName">string -- the attribute to read</param>
      /// <returns>double</returns>
      public double ReadAxisAttribute(string axName, string attributeName)
      {
        Exception ex;
        if (this.myAxisHash.ContainsKey(axName))
        {
          cAxis ax = (cAxis)myAxisHash[axName];
          switch (attributeName)
          {
            case "rotary":
              return ax.Rotary ? 1 : 0;
            case "scalefactor":
              return ax.ScaleFactor;
            case "offset":
              return ax.Offset;
            default:
              ex = new Exception("Invalid attribute name in: " + System.Reflection.MethodInfo.GetCurrentMethod().Name);
              break;
          }
        }
        else
          ex = new Exception("Invalid axis name in: " + System.Reflection.MethodInfo.GetCurrentMethod().Name);
        throw ex;
      }
      /// <summary>
      /// Reads a machine attribute, these are in the static variables section of the config file.
      /// </summary>
      /// <param name="AttribName">The name</param>
      /// <returns>a double</returns>
      public double ReadAttribute(string AttribName)
      {
        return this.myCalculator._GetVariable(AttribName);
      }
      /// <summary>
      /// Writes a machine attribute, these are in the static variables section of the config file.
      /// </summary>
      /// <param name="AttribName">string, the name</param>
      /// <param name="Value">double, the value</param>
      public void WriteAttribute(string AttribName, double Value)
      {
        Rebuild = true;
        if (!this.myCalculator._ContainsVariable(AttribName))
          System.Windows.Forms.MessageBox.Show("Variable: " + AttribName + " does not exist in the string calculator.");
        if (Math.Abs(this.myCalculator._GetVariable(AttribName) - Value) > SmallestIncrement)
        {
          myVariableToDH_Hash.MarkVariable(AttribName, myDHs);
          this.myCalculator._AssignVariable(AttribName, Value);
        }
        else
        {
          //System.Windows.Forms.MessageBox.Show("Variable " + AttribName + " not assigned a value.");
        }
      }
      /// <summary>
      /// This is the scaled, offset and full fancy comp'd axis position.  Use for calculations only!
      /// </summary>
      /// <param name="axName"></param>
      /// <returns>a double</returns>
      public double GetAxisScaledPosition(string axName)
      {
        if (this.myAxisHash.ContainsKey(axName))
        {
          cAxis ax = (cAxis)this.myAxisHash[axName];
          return ax.Value;
        }
        else
        {
          Exception ex = new Exception("Invalid axis position thrown from: " + System.Reflection.MethodInfo.GetCurrentMethod().Name);
          throw ex;
        }
      }

      public double GetAxisPostion(string axName)
      {
        if (this.myAxisHash.ContainsKey(axName))
        {
          cAxis ax = (cAxis)this.myAxisHash[axName];
          return ax.Position;
        }
        else
        {
          Exception ex = new Exception("Invalid axis position thrown from: " + System.Reflection.MethodInfo.GetCurrentMethod().Name);
          return 0.0;
          //throw ex;
        }
      }

      /// <summary>
      /// Calculates this machine's transform
      /// </summary>
      /// <returns>retruns double[4,4] representation of dotproduct with these machine axis values.</returns>
      public double[,] DotProduct()
      {
        UpdateTransforms();
        return this.myXfroms.DotProduct();
      }
      /// <summary>
      /// Returns a partial solution
      /// </summary>
      /// <param name="start">First transfrom to include in the solution (ONE BASED)</param>
      /// <param name="finish">Last transform to include in the solution</param>
      /// <returns>double[4,4] of the requested partial solution</returns>
      /// <example>cMatrix m = new cMatrix(myMachine.DotProduct(0,1)) or <br/>
      ///	double[,] dog = myMachine.DotProduct(0,1)</example>
      public double[,] DotProduct(uint start, uint finish)
      {
        UpdateTransforms();
        Rebuild = true; //in this case, we need to rebuild the system if we want to get an accurate toolpoint.
        //        return this.myXfroms.DotProduct(start - 1, finish - 1);
        // Yikes!  I added BaseShift to transforms and it is transparent.  Therefore I screwed anybody using this function.  I noticed this 2012-10-20.  TWR
        return this.myXfroms.DotProduct(start, finish);
      }
      /// <summary>
      /// Writes a configuration file representing the current state of the machine.
      /// </summary>
      /// <param name="FileOut">string representing the filename you wish to use</param>
      public void ToFile(string FileOut, bool RobotFile = false, string version = "")
      {


        bool AllowComments = false;


        if (AllowComments)
        {
          #region WithComments
          //string FileOut = "";
          string FileOutTemp = "";

          //FileOut = lblFileGenerated.Text;

          if (!System.IO.File.Exists(FileOut))
          {

            System.Windows.Forms.MessageBox.Show("Error: There is no Pre-existing Output Config File");
            return;
          }


          //get just the filename without directory
          string[] temp_ar = FileOut.Split('\\');
          string temp = temp_ar[temp_ar.Length - 1]; //the file name without directory


          //subtract the file name to get just the directory
          string DirConfig = FileOut.Replace(temp, "");
          //string DirConfig = this.lblFileGenerated.Text.Replace(temp, "");
          FileOutTemp = DirConfig + "TempConfig.xml";


          //create FileOutTemp
          //read in existing xml FileOut
          //as read in each line write it out to temp file
          //update variables as read through file
          //when finished delete FileOut and copy FileOutTemp to FileOut


          using (XmlTextWriter x = new XmlTextWriter(FileOutTemp, System.Text.Encoding.UTF8))
          {

            XmlNodeReader r;
            XmlDocument document;
            //bool rtn = true;

            if (System.IO.File.Exists(FileOut))
            {
              try
              {
                document = new XmlDocument();
                document.Load(FileOut);
                r = new XmlNodeReader(document);
              }
              catch (Exception ex)
              {
                throw (ex);
              }
              while (r.Read())
              {
                switch (r.NodeType)
                {

                  case XmlNodeType.XmlDeclaration:
                    x.Formatting = Formatting.Indented;
                    x.WriteStartDocument();
                    break;
                  case XmlNodeType.Element:


                    if (r.Name == "TransForms")
                    {
                      x.WriteStartElement("TransForms");
                      //don't call x.WriteEndElement(); because there are sub-elements within Transforms
                      break;
                    }
                    if (r.Name == "EITransform")
                    {
                      x.WriteStartElement("EITransform");
                      //don't call x.WriteEndElement(); because there are sub-elements within EITransform
                      break;
                    }
                    if (r.Name == "variables")
                    {
                      x.WriteStartElement("variables");
                      //don't call x.WriteEndElement(); because there are sub-elements within variables
                      break;
                    }


                    //robot
                    if (RobotFile)
                    {
                      if (r.Name == "robot")
                      {
                        x.WriteStartElement("robot");
                        x.WriteEndElement();
                        break;
                      }
                    }


                    //<machineaxes multiply="1000000.0">
                    if (r.Name == "machineaxes")
                    {
                      x.WriteStartElement("machineaxes");
                      x.WriteAttributeString("multiply", myCalculator._GetVariable("MachineAxisScaler").ToString("F1"));
                      //don't call x.WriteEndElement(); because there are sub-elements within machineaxes

                      break;
                    }

                    //<static multiply="10000.0">
                    if (r.Name == "static")
                    {
                      x.WriteStartElement("static");
                      x.WriteAttributeString("multiply", myCalculator._GetVariable("StaticVariableScaler").ToString("F1"));
                      //don't call x.WriteEndElement(); because there are sub-elements within static

                      x.WriteStartElement("_Debug");
                      x.WriteAttributeString("value", (Math.Round(this.myCalculator._GetVariable("_Debug")) == 1).ToString());
                      x.WriteEndElement();

                      break;
                    }

                    //static variables
                    //WriteStaticVariables(xw);
                    string[] vars = GetAttributeNames();
                    foreach (string var in vars)
                    {
                      if (var.StartsWith("BaseShift"))
                        continue;


                      if (r.Name == var)
                      {
                        x.WriteStartElement(var);
                        x.WriteAttributeString("value", this.myCalculator._GetVariable(r.Name).ToString("F4"));
                        x.WriteEndElement();
                      }
                    }

                    //static variables, BaseShift
                    if (r.Name.StartsWith("BaseShift"))
                    {
                      //continue;
                      x.WriteStartElement("BaseShift");


                      string value;
                      value = this.myCalculator._GetVariable("BaseShift_X").ToString("F4");
                      x.WriteAttributeString("X", value);
                      value = this.myCalculator._GetVariable("BaseShift_Y").ToString("F4");
                      x.WriteAttributeString("Y", value);
                      value = this.myCalculator._GetVariable("BaseShift_Z").ToString("F4");
                      x.WriteAttributeString("Z", value);
                      value = this.myCalculator._GetVariable("BaseShift_rX").ToString("F4");
                      x.WriteAttributeString("rX", value);
                      value = this.myCalculator._GetVariable("BaseShift_rY").ToString("F4");
                      x.WriteAttributeString("rY", value);
                      value = this.myCalculator._GetVariable("BaseShift_rZ").ToString("F4");
                      x.WriteAttributeString("rZ", value);

                      x.WriteEndElement();

                    }

                    //static variables, otherStaticVars
                    for (int ii = 0; ii < this.myOtherStaticVarNames.Count; ii++)
                    {
                      if (r.Name == this.myOtherStaticVarNames[ii].ToString())
                      {
                        x.WriteStartElement(this.myOtherStaticVarNames[ii].ToString());
                        x.WriteAttributeString("value", this.myCalculator._GetVariable((string)this.myOtherStaticVarNames[ii]).ToString("F4"));
                        if (this.myOtherStaticVarNames[ii].ToString().IndexOf("_") != 0)
                          x.WriteAttributeString("Hide", "true");
                        x.WriteEndElement();
                      }
                    }




                    //axis names
                    //WriteMachineAxes(x);
                    vars = GetAxisNames();
                    foreach (string var in vars)
                    {
                      if (r.Name == var)
                      {
                        x.WriteStartElement(var);
                        x.WriteAttributeString("rotary", ((cAxis)myAxisHash[var]).Rotary.ToString());
                        x.WriteAttributeString("scalefactor", ((cAxis)myAxisHash[var]).ScaleFactor.ToString("F6"));

                        x.WriteAttributeString("offset", ((cAxis)myAxisHash[var]).Offset.ToString("F6"));

                        if (((cAxis)myAxisHash[var]).Rollover == true)
                          x.WriteAttributeString("rollover", "true");
                        if (((cAxis)myAxisHash[var]).Hidden)
                          x.WriteAttributeString("Hide", "true");
                        x.WriteEndElement();
                      }
                    }

                    if (r.Name == "comptables")
                    {
                      // myCompTable.ToXmlFile(x, myCalculator);
                      x.WriteStartElement("comptables");
                      x.WriteAttributeString("linear_multiply", myCalculator._GetVariable("LinearCompScaler").ToString("F1"));
                      x.WriteAttributeString("angle_multiply", myCalculator._GetVariable("AngleCompScaler").ToString("F1"));
                      break;
                    }

                    if (r.Name == "table")
                    {
                      x.WriteStartElement("table");
                      x.WriteAttributeString("axis", r.GetAttribute("axis"));
                      //use CompTableAxis and station_count to keep track of which station is the next to write out
                      //do this because it would take a lot of work to match up all the axes values in the comp tables to the line read in
                      CompTableAxis = r.GetAttribute("axis");
                      station_count = 0;
                      break;
                    }

                    //comp table stations:
                    if (r.Name == "station")
                    {
                      //ToXmlFile2 relies on CompTableAxis and station_count
                      //updates station_count
                      station_count = myCompTable.ToXmlFile2(x, myCalculator, CompTableAxis, station_count);
                    }


                    //transform:
                    if (r.Name == "DH")
                    {
                      string DH_number = r.GetAttribute("n");
                      int index = Convert.ToInt32(DH_number);
                      if (myDHs[index] is LinearAlgebra.DHMatrix)
                      {
                        Electroimpact.LinearAlgebra.DHMatrix dh = (Electroimpact.LinearAlgebra.DHMatrix)myDHs[index];
                        x.WriteStartElement("DH");
                        {
                          x.WriteAttributeString("n", (index).ToString());
                          x.WriteAttributeString("theta", dh.myStrings.Theta);
                          x.WriteAttributeString("d", dh.myStrings.d);
                          x.WriteAttributeString("a", dh.myStrings.a);
                          x.WriteAttributeString("alpha", dh.myStrings.Alpha);
                        }
                        x.WriteEndElement();
                      }
                      else if (myDHs[index] is LinearAlgebra.string6DOF)
                      {
                        x.WriteStartElement("sixdof");
                        {
                          x.WriteAttributeString("n", (index).ToString());
                          x.WriteAttributeString("X", ((LinearAlgebra.string6DOF)myDHs[index]).myStrings.X);
                          x.WriteAttributeString("Y", ((LinearAlgebra.string6DOF)myDHs[index]).myStrings.Y);
                          x.WriteAttributeString("Z", ((LinearAlgebra.string6DOF)myDHs[index]).myStrings.Z);
                          x.WriteAttributeString("rX", ((LinearAlgebra.string6DOF)myDHs[index]).myStrings.rX);
                          x.WriteAttributeString("rY", ((LinearAlgebra.string6DOF)myDHs[index]).myStrings.rY);
                          x.WriteAttributeString("rZ", ((LinearAlgebra.string6DOF)myDHs[index]).myStrings.rZ);
                        }
                        x.WriteEndElement();
                      }
                    }



                    break;
                  case XmlNodeType.Attribute:
                    x.WriteAttributeString(r.Name, r.Value);
                    //x.WriteAttributeString(r.Name,r.ReadString);
                    break;
                  case XmlNodeType.Comment:
                    x.WriteComment(r.Value);
                    break;
                  case XmlNodeType.EndElement:
                    try
                    {
                      x.WriteEndElement();
                    }
                    catch (Exception ex)
                    {
                      r.Close();
                      x.Flush();
                      x.Close();
                      break;
                    }
                    break;
                  default:
                    break;
                }
              }
              r.Close();
            } //end of add to existing file

            //x.WriteEndElement(); // end solver_file

            x.Flush();
            x.Close();



            //delete original file
            System.IO.File.Delete(FileOut);

            //copy temp file to final file
            System.IO.File.Copy(FileOutTemp, FileOut);

            //delete temp file
            System.IO.File.Delete(FileOutTemp);
          }
          #endregion

        }
        else
        {

          #region NoComments
          using (XmlTextWriter xw = new XmlTextWriter(FileOut, System.Text.Encoding.UTF8))
          {
            xw.Formatting = Formatting.Indented;
            xw.WriteStartDocument();
            {
              xw.WriteStartElement("EITransform");
              {
                xw.WriteAttributeString("filename", FileOut);
                xw.WriteAttributeString("time", System.DateTime.Now.ToString());
                xw.WriteAttributeString("solver_version", version);

                xw.WriteStartElement("variables");
                {
                  if (RobotFile)
                  {
                    xw.WriteStartElement("robot");
                    xw.WriteEndElement();
                  }
                  xw.WriteStartElement("static");
                  {
                    xw.WriteAttributeString("multiply", myCalculator._GetVariable("StaticVariableScaler").ToString("F1"));
                    WriteStaticVariables(xw);
                  }
                  xw.WriteEndElement();

                  xw.WriteStartElement("machineaxes");
                  {
                    xw.WriteAttributeString("multiply", myCalculator._GetVariable("MachineAxisScaler").ToString("F1"));
                    WriteMachineAxes(xw);
                  }
                  xw.WriteEndElement();

                  myCompTable.ToXmlFile(xw, myCalculator);

                }
                xw.WriteEndElement();

                xw.WriteStartElement("TransForms");
                {
                  {
                    //the first DH is implied
                    for (int ii = 1; ii < myDHs.Count; ii++)
                    {
                      //If you made it here.  Todd hates this damn switch.  Good problem to solve.
                      if (myDHs[ii] is LinearAlgebra.DHMatrix)
                      {
                        Electroimpact.LinearAlgebra.DHMatrix dh = (Electroimpact.LinearAlgebra.DHMatrix)myDHs[ii];
                        xw.WriteStartElement("DH");
                        {
                          xw.WriteAttributeString("n", (ii).ToString());
                          xw.WriteAttributeString("theta", dh.myStrings.Theta);
                          xw.WriteAttributeString("d", dh.myStrings.d);
                          xw.WriteAttributeString("a", dh.myStrings.a);
                          xw.WriteAttributeString("alpha", dh.myStrings.Alpha);
                        }
                        xw.WriteEndElement();
                      }
                      else if (myDHs[ii] is LinearAlgebra.string6DOF)
                      {
                        xw.WriteStartElement("sixdof");
                        {
                          xw.WriteAttributeString("n", (ii).ToString());
                          xw.WriteAttributeString("X", ((LinearAlgebra.string6DOF)myDHs[ii]).myStrings.X);
                          xw.WriteAttributeString("Y", ((LinearAlgebra.string6DOF)myDHs[ii]).myStrings.Y);
                          xw.WriteAttributeString("Z", ((LinearAlgebra.string6DOF)myDHs[ii]).myStrings.Z);
                          xw.WriteAttributeString("rX", ((LinearAlgebra.string6DOF)myDHs[ii]).myStrings.rX);
                          xw.WriteAttributeString("rY", ((LinearAlgebra.string6DOF)myDHs[ii]).myStrings.rY);
                          xw.WriteAttributeString("rZ", ((LinearAlgebra.string6DOF)myDHs[ii]).myStrings.rZ);
                        }
                        xw.WriteEndElement();
                      }
                    }
                  }
                }
                xw.WriteEndElement();
              }
              xw.WriteEndElement();
            }
            xw.WriteEndDocument();
            xw.Flush();
            xw.Close();
          }
          #endregion
        }
      }
      /// <summary>
      /// Gets all static variable names.  For use with the solver.  During solving, these "static" variables are manipulated
      /// to match tracker data.
      /// </summary>
      /// <returns>a string array of static variable names</returns>
      public string[] GetAttributeNames()
      {

        int nAttribs = this.myStaticVarNames.Count;
        foreach (string name in this.myOtherStaticVarNames)
        {

        }
        string[] rtn = new string[this.myStaticVarNames.Count];
        for (int ii = 0; ii < rtn.Length; ii++)
        {
          rtn[ii] = (string)myStaticVarNames[ii];
        }
        return rtn;
      }
      /// <summary>
      /// Returns a string array of this machine's axis names.
      /// </summary>
      /// <returns>string array</returns>
      private string[] GetAxisNames()
      {
        string[] ret = new string[this.myMachineAxisNames.Count];
        for (int ii = 0; ii < myMachineAxisNames.Count; ii++)
        {
          ret[ii] = (string)myMachineAxisNames[ii];
        }
        return ret;
      }
      /// <summary>
      /// Returns a string array of this machine's axis names.
      /// </summary>
      /// <returns>string array</returns>
      public string[] GetAxisNames(bool IncludeHidden = true)
      {
        if (IncludeHidden)
          return GetAxisNames();
        System.Collections.Generic.List<string> names = new List<string>();
        for (int ii = 0; ii < myMachineAxisNames.Count; ii++)
        {
          if (!((cAxis)myAxisHash[myMachineAxisNames[ii]]).Hidden)
            names.Add(myMachineAxisNames[ii].ToString());
        }
        string[] ret = names.ToArray();
        return ret;
      }
      /// <summary>
      /// Indicates if an axis and the desired attribute exists
      /// </summary>
      /// <param name="AxName">string One of the axis names in the machine.</param>
      /// <param name="AttributeName">string An attribute name</param>
      /// <returns>bool indicating if the axisname, attribute name combo exists in this machine.</returns>
      public bool ContainsAxisAttribute(string AxName, string AttributeName)
      {
        if (this.myAxisHash.ContainsKey(AxName))
        {
          cAxis ax = (cAxis)myAxisHash[AxName];
          switch (AttributeName)
          {
            case "offset":
              return true;
            case "scalefactor":
              return true;
            default:
              return false;
          }
        }
        else
          return false;
      }

      public string[] GetCompensationTableAxisNames
      {
        get
        {
          return myCompTable.GetCompensationTableAxisNames;
        }
      }

      public int GetCompensationTableStationCount(string AxisName)
      {
        return myCompTable.GetCompensationTableStationCount(AxisName);
      }

      public int GetCompensationTableStationCount(string AxisName, double min, double max)
      {
        System.Collections.Generic.List<Electroimpact.Machine.cCompTable.cCompStation> t = myCompTable.GetCompensationTable(AxisName);
        int lower = -1;
        int upper = t.Count;
        for (int ii = 0; ii < t.Count; ii++)
        {
          if (t[ii].loc >= min && lower == -1)
            lower = ii;
          if (t[ii].loc <= max)
            upper = ii;
          if (t[ii].loc > max)
            break;
        }
        return upper - lower + 1;
      }

      public List<Electroimpact.Machine.cCompTable.cCompStation> GetCompensationTable(string AxisName)
      {
        return (List<Electroimpact.Machine.cCompTable.cCompStation>)myCompTable.GetCompensationTable(AxisName);
      }

      public List<Electroimpact.Machine.cCompTable.cCompStation> GetCompensationTable(string AxisName, double min, double max, out int lower, out int upper)
      {
        List<Electroimpact.Machine.cCompTable.cCompStation> t = GetCompensationTable(AxisName);
        lower = -1;
        upper = t.Count;
        for (int ii = 0; ii < t.Count; ii++)
        {
          if (t[ii].loc >= min && lower == -1)
            lower = ii;
          if (t[ii].loc <= max)
            upper = ii;
          if (t[ii].loc > max)
            break;
        }
        return t;
      }

      public void PutCompensationTable(string AxisName, List<Electroimpact.Machine.cCompTable.cCompStation> Table)
      {
        cAxis ax = (cAxis)myAxisHash[AxisName];
        if (ax.Rotary && ax.Rollover)
        {
          cCompTable.cCompStation cs0 = null;// = new cCompTable.cCompStation();
          for (int ii = 0; ii < Table.Count; ii++)
          {
            cCompTable.cCompStation cs = Table[ii];
            if (Math.Abs(cs.loc) < 1.0e-6)
            {
              cs0 = cs;
            }
            if (cs0 != null && Math.Abs(cs.loc - 360.0) < 1.0e-6)
              cs.Copy(cs0, 360.0);
          }
        }
        myCompTable.PutCompensationTable(AxisName, Table);
      }

      public bool CompTableOverride
      {
        get { return this.comptableoverride; }
        set { this.comptableoverride = value; }
      }

      /// <summary>
      /// Puts a number in the ...
      /// </summary>
      /// <param name="CompTableVaiableName"></param>
      /// <param name="value"></param>
      public void OverrideCompTableValues(string CompTableVaiableName, double value)
      {
        myCalculator._AssignVariable(CompTableVaiableName, value);
      }

      public string[] CompTableVariableNames
      {
        get { return myCompTable.CompTableVariableNames; }
      }

      public void StoreToCNC(string cnc_address)
      {
        //Electroimpact.FANUC.Err_Code err;
        //Electroimpact.FANUC.OpenCNC c = new Electroimpact.FANUC.OpenCNC(cnc_address, out err);

      }

      public string[] ToString()
      {
        System.Collections.Generic.List<string> positions = new List<string>();
        positions.Add(" X: " + X.ToString("F4").PadLeft(10));
        positions.Add(" Y: " + Y.ToString("F4").PadLeft(10));
        positions.Add(" Z: " + Z.ToString("F4").PadLeft(10));
        positions.Add(" A: " + A.ToString("F4").PadLeft(10));
        positions.Add(" B: " + B.ToString("F4").PadLeft(10));
        positions.Add(" C: " + C[0].ToString("F4").PadLeft(10));
        positions.Add("rX: " + rZrYrX[0, 0].ToString("F4").PadLeft(10));
        positions.Add("rY: " + rZrYrX[1, 0].ToString("F4").PadLeft(10));
        positions.Add("rZ: " + rZrYrX[2, 0].ToString("F4").PadLeft(10));
        return positions.ToArray();
      }

      #endregion

      public void GenerateCcode(string output_filename)
      {


        throw new NotImplementedException();
      }
    }
    #endregion

    #region CompTable
    public class cCompTable
    {
      #region Members and Constructor

      #region cCompStation
      public class cCompStation
      {
        public double loc;
        public double X = 0;
        public double Y = 0;
        public double Z = 0;
        public double rX = 0;
        public double rY = 0;
        public double rZ = 0;

        public cCompStation()
        { }


        public cCompStation(cCompStation csin)
        {
          loc = csin.loc;
          X = csin.X;
          Y = csin.Y;
          Z = csin.Z;
          rX = csin.rX;
          rY = csin.rY;
          rZ = csin.rZ;
        }

        public void Copy(cCompStation csin, double location)
        {
          loc = location;
          X = csin.X;
          Y = csin.Y;
          Z = csin.Z;
          rX = csin.rX;
          rY = csin.rY;
          rZ = csin.rZ;
        }

        public override string ToString()
        {
          string ret = loc.ToString("F4") + ", " +
                       X.ToString("F4") + ", " +
                       Y.ToString("F4") + ", " +
                       Z.ToString("F4") + ", " +
                       rX.ToString("F4") + ", " +
                       rY.ToString("F4") + ", " +
                       rZ.ToString("F4");
          return ret;
        }
      }
      #endregion

      //private System.Collections.Generic.KeyValuePair<string, System.Collections.Generic.List> myCompTables = new KeyValuePair<string, cCompStation>();
      private Hashtable myCompTables = new Hashtable();
      private System.Collections.ArrayList myAxisKeys = new ArrayList();

      public cCompTable()
      {
      }

      #endregion

      #region Publics

      public void Clear()
      {
        myCompTables.Clear();
        myAxisKeys.Clear();
      }

      public void NullTableValues()
      {
        for (int ii = 0; ii < myAxisKeys.Count; ii++)
        {
          List<cCompStation> l = (List<cCompStation>)myCompTables[myAxisKeys[ii]];
          for (int jj = 0; jj < l.Count; jj++)
          {
            l[jj].X = 0;
            l[jj].Y = 0;
            l[jj].Z = 0;
            l[jj].rX = 0;
            l[jj].rY = 0;
            l[jj].rZ = 0;
          }
        }
      }

      /// <summary>
      /// Returns an ArrayList of strings comprised of the variable names added.
      /// </summary>
      /// <param name="r"></param>
      /// <returns></returns>
      public ArrayList AddTable(XmlNodeReader r)
      {
        string axis = r.GetAttribute("axis");
        List<cCompStation> Table = new List<cCompStation>();
        System.Collections.ArrayList VariableNames = new ArrayList();
        while (r.Read())
        {
          if (r.NodeType == XmlNodeType.EndElement)
          {
            if (r.Name != "table")
            {
              Exception ex = new Exception("XML Format Error!");
              throw (ex);
            }
            break;
          }
          if (r.NodeType == XmlNodeType.Element)
          {
            if (r.Name == "station")
            {
              cCompStation s = new cCompStation();
              s.loc = double.TryParse(r.GetAttribute("loc"), out s.loc) ? s.loc : 0;
              s.X = double.TryParse(r.GetAttribute("X"), out s.X) ? s.X : 0;
              s.Y = double.TryParse(r.GetAttribute("Y"), out s.Y) ? s.Y : 0;
              s.Z = double.TryParse(r.GetAttribute("Z"), out s.Z) ? s.Z : 0;
              s.rX = double.TryParse(r.GetAttribute("rX"), out s.rX) ? s.rX : 0;
              s.rY = double.TryParse(r.GetAttribute("rY"), out s.rY) ? s.rY : 0;
              s.rZ = double.TryParse(r.GetAttribute("rZ"), out s.rZ) ? s.rZ : 0;
              Table.Add(s);
            }
          }
        }
        string key;

        if (myAxisKeys.Contains(axis))
          myCompTables[axis] = Table;
        else
        {
          myAxisKeys.Add(axis);
          myCompTables.Add(axis, Table);
          key = "l" + axis + "x";
          VariableNames.Add(key);

          key = "l" + axis + "y";
          VariableNames.Add(key);

          key = "l" + axis + "z";
          VariableNames.Add(key);

          key = "r" + axis + "x";
          VariableNames.Add(key);

          key = "r" + axis + "y";
          VariableNames.Add(key);

          key = "r" + axis + "z";
          VariableNames.Add(key);
        }

        return VariableNames;
      }

      public bool ContainsAxis(string axisname)
      {
        return myCompTables.ContainsKey(axisname);
      }

      /// <summary>
      /// Used to return the linear interpolated values of X,Y,Z,Roll,Pitch,Yaw of the comptable for the specified axis at the specified location.
      /// </summary>
      /// <param name="axisname">string...the axis name</param>
      /// <param name="axisposition">double the axis' position</param>
      /// <returns>Hashtable with key of l(axisname)x (y and z) and r(axisname)r (p and y)</returns>
      public System.Collections.Hashtable GetCompVariables(string axisname, double axisposition)
      {
        if (!myCompTables.ContainsKey(axisname))
        {
          throw new System.Exception("This comptable does not contain a table for the " + axisname + " axis!");
        }
        cCompStation cs = CreateCompStation((List<cCompStation>)myCompTables[axisname], axisposition);
        System.Collections.Hashtable rtn = new Hashtable();

        string key;
        key = "l" + axisname + "x";
        rtn.Add(key, cs.X);

        key = "l" + axisname + "y";
        rtn.Add(key, cs.Y);

        key = "l" + axisname + "z";
        rtn.Add(key, cs.Z);

        key = "r" + axisname + "x";
        rtn.Add(key, cs.rX);

        key = "r" + axisname + "y";
        rtn.Add(key, cs.rY);

        key = "r" + axisname + "z";
        rtn.Add(key, cs.rZ);
        return rtn;
      }

      public int ToXmlFile2(XmlTextWriter xw, Electroimpact.StringCalc.cStringCalc calc, string AxisIn, int ScountIn)
      {
        for (int ii = 0; ii < myAxisKeys.Count; ii++)
        {
          if (AxisIn == (string)myAxisKeys[ii])
          {
            List<Electroimpact.Machine.cCompTable.cCompStation> Table = (List<Electroimpact.Machine.cCompTable.cCompStation>)myCompTables[(string)myAxisKeys[ii]];
            int jj = ScountIn;
            cCompStation cs = (cCompStation)Table[jj];
            xw.WriteStartElement("station");
            {
              xw.WriteAttributeString("loc", cs.loc.ToString("F4"));
              xw.WriteAttributeString("X", cs.X.ToString("F4"));
              xw.WriteAttributeString("Y", cs.Y.ToString("F4"));
              xw.WriteAttributeString("Z", cs.Z.ToString("F4"));
              xw.WriteAttributeString("rX", cs.rX.ToString("F10"));
              xw.WriteAttributeString("rY", cs.rY.ToString("F10"));
              xw.WriteAttributeString("rZ", cs.rZ.ToString("F10"));
            }
            xw.WriteEndElement();
          }
        }
        int temp = ScountIn + 1; //increment temp var because can't change ScountIn
        return temp;
      }


      public void ToXmlFile(XmlTextWriter xw, Electroimpact.StringCalc.cStringCalc calc)
      {
        if (myCompTables.Count > 0)
        {
          //string[] axisnames = new string[myCompTables.Count];
          //myCompTables.Keys.CopyTo(axisnames, 0);
          xw.WriteStartElement("comptables");
          {
            xw.WriteAttributeString("linear_multiply", calc._GetVariable("LinearCompScaler").ToString("F1"));
            xw.WriteAttributeString("angle_multiply", calc._GetVariable("AngleCompScaler").ToString("F1"));

            for (int ii = 0; ii < myAxisKeys.Count; ii++)
            {
              xw.WriteStartElement("table");
              {
                xw.WriteAttributeString("axis", (string)myAxisKeys[ii]);
                List<Electroimpact.Machine.cCompTable.cCompStation> Table = (List<Electroimpact.Machine.cCompTable.cCompStation>)myCompTables[(string)myAxisKeys[ii]];
                for (int jj = 0; jj < Table.Count; jj++)
                {
                  xw.WriteStartElement("station");
                  {
                    cCompStation cs = (cCompStation)Table[jj];
                    xw.WriteAttributeString("loc", cs.loc.ToString("F4"));
                    xw.WriteAttributeString("X", cs.X.ToString("F4"));
                    xw.WriteAttributeString("Y", cs.Y.ToString("F4"));
                    xw.WriteAttributeString("Z", cs.Z.ToString("F4"));
                    xw.WriteAttributeString("rX", cs.rX.ToString("F10"));
                    xw.WriteAttributeString("rY", cs.rY.ToString("F10"));
                    xw.WriteAttributeString("rZ", cs.rZ.ToString("F10"));
                  }
                  xw.WriteEndElement();
                }
              }
              xw.WriteEndElement();
            }
          }
          xw.WriteEndElement();
        }
      }

      public string[] GetCompensationTableAxisNames
      {
        get
        {
          string[] names = new string[myAxisKeys.Count];
          for (int ii = 0; ii < myAxisKeys.Count; ii++)
          {
            names[ii] = (string)myAxisKeys[ii];
          }
          return names;
        }
      }

      public string[] CompTableVariableNames
      {
        get
        {
          string[] axes = GetCompensationTableAxisNames;
          string[] r = new string[axes.Length * 6];
          for (int ii = 0; ii < axes.Length; ii++)
          {
            r[ii * 6 + 0] = "l" + axes[ii] + "x";
            r[ii * 6 + 1] = "l" + axes[ii] + "y";
            r[ii * 6 + 2] = "l" + axes[ii] + "z";
            r[ii * 6 + 3] = "r" + axes[ii] + "x";
            r[ii * 6 + 4] = "r" + axes[ii] + "y";
            r[ii * 6 + 5] = "r" + axes[ii] + "z";
          }
          return r;
        }
      }

      /// <summary>
      /// Returns the number of cCompStations in the comptable for this AxisName.
      /// </summary>
      /// <param name="AxisName"></param>
      /// <returns></returns>
      public int GetCompensationTableStationCount(string AxisName)
      {
        if (ContainsAxis(AxisName))
        {
          return ((List<cCompStation>)myCompTables[AxisName]).Count;
        }
        else
          throw new System.Exception("AxisName: " + AxisName + " is in valid.");
      }

      /// <summary>
      /// Used in Solver only!
      /// </summary>
      /// <param name="AxisName">string The axis name</param>
      /// <returns>System.Collections.ArrayList of cCompStations</returns>
      public List<Electroimpact.Machine.cCompTable.cCompStation> GetCompensationTable(string AxisName)
      {
        if (ContainsAxis(AxisName))
        {
          return (List<Electroimpact.Machine.cCompTable.cCompStation>)myCompTables[AxisName];
        }
        else
          throw new Exception("error in: " + System.Reflection.MethodInfo.GetCurrentMethod().ToString());
      }

      /// <summary>
      /// Only a solver funtion.
      /// </summary>
      /// <param name="AxisName">string</param>
      /// <param name="table">a System.Collections.ArrayList of cCompStation (s). used to modify the axis comp from solver.</param>
      public void PutCompensationTable(string AxisName, List<Electroimpact.Machine.cCompTable.cCompStation> table)
      {
        if (ContainsAxis(AxisName))
        {
          this.myCompTables[AxisName] = table;
        }
        else
          throw new Exception("error in: " + System.Reflection.MethodInfo.GetCurrentMethod().ToString());
      }

      #endregion

      #region Privates

      private cCompStation CreateCompStation(List<cCompStation> Stations, double AxVal)
      {
        int test = Stations.Count / 2;
        int lastBigger = Stations.Count;
        while (true)
        {
          cCompStation cs = Stations[test];
          cCompStation csp = Stations[test - 1];
          if (AxVal <= cs.loc)
          {
            if (AxVal > csp.loc)
            {//You got it!
              break;
            }
            else
            {
              lastBigger = test;
              test /= 2;
              if (test == 0)
                break; //Failed to the negative
            }
          }
          else
          {
            if (test == Stations.Count - 1)
            {
              return (cCompStation)Stations[Stations.Count - 1]; //Off the table to the positive.
            }
            int dist = (lastBigger - test) / 2;
            dist = dist > 1 ? dist : 1;
            test += dist;
          }
        }
        if (test > 0)
        {
          cCompStation cs = (cCompStation)Stations[test];
          cCompStation csp = (cCompStation)Stations[test - 1];
          cCompStation rtn = new cCompStation();

          rtn.X = LinearInterpolate(AxVal, cs.loc, csp.loc, cs.X, csp.X);
          rtn.Y = LinearInterpolate(AxVal, cs.loc, csp.loc, cs.Y, csp.Y);
          rtn.Z = LinearInterpolate(AxVal, cs.loc, csp.loc, cs.Z, csp.Z);
          rtn.rX = LinearInterpolate(AxVal, cs.loc, csp.loc, cs.rX, csp.rX);
          rtn.rY = LinearInterpolate(AxVal, cs.loc, csp.loc, cs.rY, csp.rY);
          rtn.rZ = LinearInterpolate(AxVal, cs.loc, csp.loc, cs.rZ, csp.rZ);
          return rtn;
        }
        else
          return (cCompStation)Stations[0]; //Off the table to the negative.
      }

      private double LinearInterpolate(double loc, double X1, double X2, double Y1, double Y2)
      {
        double m = (Y2 - Y1) / (X2 - X1);
        double b = Y1 - m * X1;
        return m * loc + b;
      }

      #endregion

    }
    #endregion

    #region Robot Class
    /// <summary>
    /// Things to do:
    /// 
    /// 1) PayloadCalc needs to be passed parameters from setup file.
    /// </summary>
    public class cRobot : cMachine//, IMachine
    {
      #region Members and Constructor
      PayloadCalc myPayloadCalc;
      public cRobot(string FileName, bool UseDegrees)
      {
        this.FileName = FileName;
        this.HookedUp = base.HookMeUp(FileName);
        this.myCalculator.Degrees = UseDegrees;

        #region setup myPayloadCalc
        myPayloadCalc = new PayloadCalc();
        if (myCalculator._ContainsVariable("_m2"))
          myPayloadCalc.m2 = myCalculator._GetVariable("_m2");
        if (myCalculator._ContainsVariable("_m3"))
          myPayloadCalc.m3 = myCalculator._GetVariable("_m3");
        if (myCalculator._ContainsVariable("_mH"))
          myPayloadCalc.mH = myCalculator._GetVariable("_mH");
        if (myCalculator._ContainsVariable("_mP"))
          myPayloadCalc.mP = myCalculator._GetVariable("_mP");
        if (myCalculator._ContainsVariable("_ls2"))
          myPayloadCalc.ls2 = myCalculator._GetVariable("_ls2");
        if (myCalculator._ContainsVariable("_ls3"))
          myPayloadCalc.ls3 = myCalculator._GetVariable("_ls3");
        if (myCalculator._ContainsVariable("_lsH"))
          myPayloadCalc.lsH = myCalculator._GetVariable("_lsH");
        if (myCalculator._ContainsVariable("_lsPx"))
          myPayloadCalc.lsPx = myCalculator._GetVariable("_lsPx");
        if (myCalculator._ContainsVariable("_lsPy"))
          myPayloadCalc.lsPy = myCalculator._GetVariable("_lsPy");
        if (myCalculator._ContainsVariable("_lsPz"))
          myPayloadCalc.lsPz = myCalculator._GetVariable("_lsPz");
        if (myCalculator._ContainsVariable("J1z"))
          myPayloadCalc.J1z = myCalculator._GetVariable("J1z");
        if (myCalculator._ContainsVariable("J4z"))
          myPayloadCalc.J4z = myCalculator._GetVariable("J4z");
        if (myCalculator._ContainsVariable("J5z"))
          myPayloadCalc.J5z = myCalculator._GetVariable("J5z");
        if (myCalculator._ContainsVariable("J1y"))
          myPayloadCalc.J1y = myCalculator._GetVariable("J1y");
        if (myCalculator._ContainsVariable("J2z"))
          myPayloadCalc.J2z = myCalculator._GetVariable("J2z");
        if (myCalculator._ContainsVariable("J3y"))
          myPayloadCalc.J3y = myCalculator._GetVariable("J3y");
        #endregion
      }

      protected override void UpdateTransforms()
      {
        if (base.Rebuild)
        {
          PayloadCalc p = myPayloadCalc;
          double[] jointPositions = new double[6];
          jointPositions[0] = base.GetAxisPostion("J1");
          jointPositions[1] = base.GetAxisPostion("J2");
          jointPositions[2] = base.GetAxisPostion("J3");
          jointPositions[3] = base.GetAxisPostion("J4");
          jointPositions[4] = base.GetAxisPostion("J5");
          jointPositions[5] = base.GetAxisPostion("J6");
          double[] moments = p.CalcMoments(jointPositions);
          for (int ii = 0; ii < moments.Length; ii++)
            moments[ii] /= 1.0e6;
          //Console.WriteLine(moments[0].ToString("F3"));
          //Console.WriteLine(moments[1].ToString("F3"));
          //Console.WriteLine(moments[2].ToString("F3"));
          //Console.WriteLine(moments[3].ToString("F3"));
          //Console.WriteLine(moments[4].ToString("F3"));
          //Console.WriteLine(moments[5].ToString("F3"));
          base.WriteAxisPosition("tJ1", moments[0]);
          base.WriteAxisPosition("tJ2", moments[1]);
          base.WriteAxisPosition("tJ3", moments[2]);
          base.WriteAxisPosition("tJ4", moments[3]);
          base.WriteAxisPosition("tJ5", moments[4]);
          base.WriteAxisPosition("tJ6", moments[5]);
          base.UpdateTransforms();
        }
      }

      #endregion

      #region Russ's Moment Software
      public class PayloadCalc
      {
        public double J1z { set { x1 = value; } }
        private double x1 = 500;//J1z												//KR360/1 link parameters in mm
        public double J4z { set { x2 = value; } }
        private double x2 = 1025;//J4z
        public double J5z { set { x3 = value; } }
        private double x3 = 290;//J5z
        public double J1y { set { y1 = value; } }
        private double y1 = 1045;//J1y
        public double J2z { set { y2 = value; } }
        private double y2 = 1300;//J2z
        public double J3y { set { y3 = value; } }
        private double y3 = 55;//J3y
        public double m2 = 0;	//434kg config file
        public double m3 = 0; //570kg config file
        public double mH = 0;	//194kg config file
        public double mP = 352; // 273;	//124kg Joe's head Pounds config file
        public double ls2 = 500;  //link 2 cg distance config file
        public double ls3 = -326;  //link 3 cg distance config file
        public double lsH = 225;  //link wrist cg distance config file
        public double lsPx = 126;  //end-effector cg config file
        public double lsPy = -14;  //end-effector cg config file
        public double lsPz = 155;  //end-effector cg config file

        #region Some Classes From Russ's
        public class Matrix
        {
          private double[,] M = new double[4, 4];
          public Matrix()
          {

          }
          public Matrix(double i1, double j1, double k1, double r1,
                        double i2, double j2, double k2, double r2,
                        double i3, double j3, double k3, double r3,
                        double i4, double j4, double k4, double r4)
          {
            M[0, 0] = i1;
            M[0, 1] = j1;
            M[0, 2] = k1;
            M[0, 3] = r1;

            M[1, 0] = i2;
            M[1, 1] = j2;
            M[1, 2] = k2;
            M[1, 3] = r2;

            M[2, 0] = i3;
            M[2, 1] = j3;
            M[2, 2] = k3;
            M[2, 3] = r3;

            M[3, 0] = i4;
            M[3, 1] = j4;
            M[3, 2] = k4;
            M[3, 3] = r4;
          }
          public Matrix(double[,] VarArray)
          {
            for (int ii = 0; ii < 4; ii++)
            {
              for (int jj = 0; jj < 4; jj++)
              {
                M[ii, jj] = VarArray[ii, jj];
              }
            }
          }

          public double this[int row, int col]
          {
            get
            {
              return M[row, col];
            }
            set
            {
              M[row, col] = value;
            }
          }
          public static Matrix operator *(Matrix a, Matrix b)
          {
            Matrix S = new Matrix();
            double[,] c = a.M;
            double[,] d = b.M;

            S[0, 0] = d[0, 0] * c[0, 0]
                    + d[1, 0] * c[0, 1]
                    + d[2, 0] * c[0, 2]
                    + d[3, 0] * c[0, 3];

            S[1, 0] = d[0, 0] * c[1, 0]
                    + d[1, 0] * c[1, 1]
                    + d[2, 0] * c[1, 2]
                    + d[3, 0] * c[1, 3];

            S[2, 0] = d[0, 0] * c[2, 0]
                    + d[1, 0] * c[2, 1]
                    + d[2, 0] * c[2, 2]
                    + d[3, 0] * c[2, 3];

            S[3, 0] = d[0, 0] * c[3, 0]
                    + d[1, 0] * c[3, 1]
                    + d[2, 0] * c[3, 2]
                    + d[3, 0] * c[3, 3];

            S[0, 1] = d[0, 1] * c[0, 0]
                    + d[1, 1] * c[0, 1]
                    + d[2, 1] * c[0, 2]
                    + d[3, 1] * c[0, 3];

            S[1, 1] = d[0, 1] * c[1, 0]
                    + d[1, 1] * c[1, 1]
                    + d[2, 1] * c[1, 2]
                    + d[3, 1] * c[1, 3];

            S[2, 1] = d[0, 1] * c[2, 0]
                    + d[1, 1] * c[2, 1]
                    + d[2, 1] * c[2, 2]
                    + d[3, 1] * c[2, 3];

            S[3, 1] = d[0, 1] * c[3, 0]
                    + d[1, 1] * c[3, 1]
                    + d[2, 1] * c[3, 2]
                    + d[3, 1] * c[3, 3];

            S[0, 2] = d[0, 2] * c[0, 0]
                    + d[1, 2] * c[0, 1]
                    + d[2, 2] * c[0, 2]
                    + d[3, 2] * c[0, 3];

            S[1, 2] = d[0, 2] * c[1, 0]
                    + d[1, 2] * c[1, 1]
                    + d[2, 2] * c[1, 2]
                    + d[3, 2] * c[1, 3];

            S[2, 2] = d[0, 2] * c[2, 0]
                    + d[1, 2] * c[2, 1]
                    + d[2, 2] * c[2, 2]
                    + d[3, 2] * c[2, 3];

            S[3, 2] = d[0, 2] * c[3, 0]
                    + d[1, 2] * c[3, 1]
                    + d[2, 2] * c[3, 2]
                    + d[3, 2] * c[3, 3];

            S[0, 3] = d[0, 3] * c[0, 0]
                    + d[1, 3] * c[0, 1]
                    + d[2, 3] * c[0, 2]
                    + d[3, 3] * c[0, 3];

            S[1, 3] = d[0, 3] * c[1, 0]
                    + d[1, 3] * c[1, 1]
                    + d[2, 3] * c[1, 2]
                    + d[3, 3] * c[1, 3];

            S[2, 3] = d[0, 3] * c[2, 0]
                    + d[1, 3] * c[2, 1]
                    + d[2, 3] * c[2, 2]
                    + d[3, 3] * c[2, 3];

            S[3, 3] = d[0, 3] * c[3, 0]
                    + d[1, 3] * c[3, 1]
                    + d[2, 3] * c[3, 2]
                    + d[3, 3] * c[3, 3];
            return S;
          }
          public static Matrix operator *(Matrix a, double b)
          {
            a[0, 0] = a[0, 0] * b;
            a[0, 1] = a[0, 1] * b;
            a[0, 2] = a[0, 2] * b;
            a[0, 3] = a[0, 3] * b;

            a[1, 0] = a[1, 0] * b;
            a[1, 1] = a[1, 1] * b;
            a[1, 2] = a[1, 2] * b;
            a[1, 3] = a[1, 3] * b;

            a[2, 0] = a[2, 0] * b;
            a[2, 1] = a[2, 1] * b;
            a[2, 2] = a[2, 2] * b;
            a[2, 3] = a[2, 3] * b;

            a[3, 0] = a[3, 0] * b;
            a[3, 1] = a[3, 1] * b;
            a[3, 2] = a[3, 2] * b;
            a[3, 3] = a[3, 3] * b;
            return a;
          }
          public Vector i()
          {
            Vector p = new Vector(this.M[0, 0], this.M[1, 0], this.M[2, 0]);
            return p;
          }
          public Vector j()
          {
            Vector p = new Vector(M[0, 1], M[1, 1], M[2, 1]);
            return p;
          }
          public Vector k()
          {
            Vector p = new Vector(M[0, 2], M[1, 2], M[2, 2]);
            return p;
          }
          public Vector r()
          {
            Vector p = new Vector(M[0, 3], M[1, 3], M[2, 3]);
            return p;
          }
          public static Matrix Inv(Matrix a)
          {
            Matrix b = new Matrix();
            Vector ar = a.r();
            Vector ai = a.i();
            Vector aj = a.j();
            Vector ak = a.k();
            double d = Vector.dot(ar, ai);
            double e = Vector.dot(ar, aj);
            double f = Vector.dot(ar, ak);
            b[0, 0] = a[0, 0];
            b[0, 1] = a[1, 0];
            b[0, 2] = a[2, 0];
            b[0, 3] = -d;
            b[1, 0] = a[0, 1];
            b[1, 1] = a[1, 1];
            b[1, 2] = a[2, 1];
            b[1, 3] = -e;
            b[2, 0] = a[0, 2];
            b[2, 1] = a[1, 2];
            b[2, 2] = a[2, 2];
            b[2, 3] = -f;
            b[3, 0] = 0;
            b[3, 1] = 0;
            b[3, 2] = 0;
            b[3, 3] = 1;
            return b;
          }
          public Pose XYZABC()
          {
            Pose J = new Pose();
            double degrees = 180 / Math.PI;

            J[0] = this.M[0, 3];
            J[1] = this.M[1, 3];
            J[2] = this.M[2, 3];
            if (Math.Abs(this.M[0, 0]) < Math.Pow(2, -52) && Math.Abs(this.M[1, 0]) < Math.Pow(2, -52))
            {
              J[3] = 0;
              J[4] = Math.Atan2(-this.M[2, 0], this.M[0, 0]) * degrees;
              J[5] = Math.Atan2(-this.M[1, 2], this.M[1, 1]) * degrees;
            }
            else
            {
              J[3] = Math.Atan2(this.M[1, 0], this.M[0, 0]) * degrees;
              double sp = Math.Sin(J[3] / degrees);
              double cp = Math.Cos(J[3] / degrees);
              J[4] = Math.Atan2(-this.M[2, 0], cp * this.M[0, 0] + sp * this.M[1, 0]) * degrees;
              J[5] = Math.Atan2(sp * this.M[0, 2] - cp * this.M[1, 2], cp * this.M[1, 1] - sp * this.M[0, 1]) * degrees;
            }
            return J;
          }
        }
        public class Vector
        {
          private double[] V = new double[3];
          public Vector(double a, double b, double c)
          {
            V[0] = a;
            V[1] = b;
            V[2] = c;
          }
          public Vector()
          {

          }
          public double i()
          {
            double s = this.V[0];
            return s;
          }
          public double j()
          {
            double s = this.V[1];
            return s;
          }
          public double k()
          {
            double s = this.V[2];
            return s;
          }
          public static Vector operator *(double s, Vector a)
          {
            a[0] = a[0] * s;
            a[1] = a[1] * s;
            a[2] = a[2] * s;

            return a;
          }
          public static Vector operator /(Vector a, double s)
          {
            a[0] = a[0] / s;
            a[1] = a[1] / s;
            a[2] = a[2] / s;

            return a;
          }
          public static Vector operator +(Vector a, Vector b)
          {
            Vector Peanut = new Vector();
            Peanut[0] = a[0] + b[0];
            Peanut[1] = a[1] + b[1];
            Peanut[2] = a[2] + b[2];
            return Peanut;
          }
          public static Vector operator -(Vector a, Vector b)
          {
            Vector Peanut = new Vector();
            Peanut[0] = a[0] - b[0];
            Peanut[1] = a[1] - b[1];
            Peanut[2] = a[2] - b[2];
            return Peanut;
          }
          public static double dot(Vector a, Vector b)
          {
            double s;
            s = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
            return s;
          }
          public Vector Cross(Vector dog)
          {
            Vector Peanut = new Vector(
                          V[1] * dog[2] - dog[1] * V[2],
                          V[2] * dog[0] - dog[2] * V[0],
                          V[0] * dog[1] - dog[0] * V[1]);
            return Peanut;
          }
          public double this[int index]
          {
            get
            {
              return V[index];
            }
            set
            {
              V[index] = value;
            }
          }
        }
        public class Pose
        {
          private double[] P = new double[6];
          public Pose(double xx, double yy, double zz, double aa, double bb, double cc)
          {
            P[0] = xx;
            P[1] = yy;
            P[2] = zz;
            P[3] = aa;
            P[4] = bb;
            P[5] = cc;
          }
          public Pose()
          {

          }
          public double this[int index]
          {
            get
            {
              return P[index];
            }
            set
            {
              P[index] = value;
            }
          }
          public Matrix TR()
          {
            Matrix M = new Matrix();
            double Rz = this.P[3] * Math.PI / 180;
            double Ry = this.P[4] * Math.PI / 180;
            double Rx = this.P[5] * Math.PI / 180;
            double sinRz = Math.Sin(Rz);
            double sinRy = Math.Sin(Ry);
            double sinRx = Math.Sin(Rx);
            double cosRz = Math.Cos(Rz);
            double cosRy = Math.Cos(Ry);
            double cosRx = Math.Cos(Rx);
            double x = this.P[0];
            double y = this.P[1];
            double z = this.P[2];

            M[0, 0] = cosRz * cosRy;
            M[0, 1] = cosRz * sinRy * sinRx - sinRz * cosRx;
            M[0, 2] = cosRz * sinRy * cosRx + sinRz * sinRx;
            M[0, 3] = x;
            M[1, 0] = sinRz * cosRy;
            M[1, 1] = sinRz * sinRy * sinRx + cosRz * cosRx;
            M[1, 2] = sinRz * sinRy * cosRx - cosRz * sinRx;
            M[1, 3] = y;
            M[2, 0] = -sinRy;
            M[2, 1] = cosRy * sinRx;
            M[2, 2] = cosRy * cosRx;
            M[2, 3] = z;
            M[3, 0] = 0;
            M[3, 1] = 0;
            M[3, 2] = 0;
            M[3, 3] = 1;

            return M;
          }

        }
        public struct SJoint
        {
          public double A1;
          public double A2;
          public double A3;
          public double A4;
          public double A5;
          public double A6;
        }

        #endregion

        private static SJoint jointAngles;
        private static double[] JointMoments = new double[6];

        public double[] CalcMoments(double[] JointAngles)
        {

          jointAngles.A1 = JointAngles[0];
          jointAngles.A2 = JointAngles[1];
          jointAngles.A3 = JointAngles[2];
          jointAngles.A4 = JointAngles[3];
          jointAngles.A5 = JointAngles[4];
          jointAngles.A6 = JointAngles[5];
          DoTheMath();
          return JointMoments;
        }

        static double rad = Math.PI / 180.0;
        static double[,] DataWC2Robot = { { 0, -1, 0, 0 }, { 0, 0, 1, 0 }, { -1, 0, 0, 0 }, { 0, 0, 0, 1 } };
        static double[,] DataTFacePlate = { { -1, 0, 0, 0 }, { 0, 1, 0, 0 }, { 0, 0, -1, 0 }, { 0, 0, 0, 1 } };
        static double[,] DataTTDXYZ = new double[4, 4]; //{ { 1, 0, 0, toolDef.X }, { 0, 1, 0, toolDef.Y }, { 0, 0, 1, toolDef.Z }, { 0, 0, 0, 1 } };
        static double[,] DataTTDRotZ = new double[4, 4]; //{ { Math.Cos(toolDef.Roll), -Math.Sin(toolDef.Roll), 0, 0 }, { Math.Sin(toolDef.Roll), Math.Cos(toolDef.Roll), 0, 0 }, { 0, 0, 1, 0 }, { 0, 0, 0, 1 } };
        static double[,] DataTTDRotY = new double[4, 4]; //{ { Math.Cos(toolDef.Pitch), 0, Math.Sin(toolDef.Pitch), 0 }, { 0, 1, 0, 0 }, { -Math.Sin(toolDef.Pitch), 0, Math.Cos(toolDef.Pitch), 0 }, { 0, 0, 0, 1 } };
        static double[,] DataTTDRotX = new double[4, 4]; //{ { 1, 0, 0, 0 }, { 0, Math.Cos(toolDef.Yaw), -Math.Sin(toolDef.Yaw), 0 }, { 0, Math.Sin(toolDef.Yaw), Math.Cos(toolDef.Yaw), 0 }, { 0, 0, 0, 1 } };

        private void DoTheMath()
        {
          #region Variables

          double mA2 = mP + mH + m3 + m2;
          double mA3 = mP + mH + m3;
          double mA4 = mP;
          double mA5 = mP;
          double mA6 = mP;

          #endregion

          double A1 = jointAngles.A1 * rad;				//Joint angle expressed in rad
          double A2 = jointAngles.A2 * rad;
          double A3 = jointAngles.A3 * rad;
          double A4 = jointAngles.A4 * rad;
          double A5 = jointAngles.A5 * rad;
          double A6 = jointAngles.A6 * rad;

          #region Build Matricies

          //refer to Paul, "Robot Manipulators", 1981 for roll/pitch/yaw transform method
          //For derivation of robot transformations, contact Russ DeVlieg @ Electroimpact and/or consult the project's technical file

          double[,] DataTb1 = { { Math.Cos(A1), -Math.Sin(A1), 0, 0 }, { -Math.Sin(A1), -Math.Cos(A1), 0, 0 }, { 0, 0, -1, 0 }, { 0, 0, 0, 1 } };
          double[,] DataT12 = { { Math.Cos(A2), -Math.Sin(A2), 0, x1 }, { 0, 0, -1, 0 }, { Math.Sin(A2), Math.Cos(A2), 0, -y1 }, { 0, 0, 0, 1 } };
          double[,] DataT23 = { { Math.Cos(A3), -Math.Sin(A3), 0, y2 }, { Math.Sin(A3), Math.Cos(A3), 0, 0 }, { 0, 0, 1, 0 }, { 0, 0, 0, 1 } };
          double[,] DataT34 = { { 0, 0, -1, 0 }, { -Math.Cos(A4), Math.Sin(A4), 0, -y3 }, { Math.Sin(A4), Math.Cos(A4), 0, 0 }, { 0, 0, 0, 1 } };
          double[,] DataT45 = { { -Math.Sin(A5), -Math.Cos(A5), 0, 0 }, { 0, 0, 1, 0 }, { -Math.Cos(A5), Math.Sin(A5), 0, -x2 }, { 0, 0, 0, 1 } };
          double[,] DataT56 = { { 0, 0, -1, x3 }, { -Math.Cos(A6), Math.Sin(A6), 0, 0 }, { Math.Sin(A6), Math.Cos(A6), 0, 0 }, { 0, 0, 0, 1 } };

          double[,] DataT2m2 = { { 1, 0, 0, ls2 }, { 0, 1, 0, 0 }, { 0, 0, 1, 0 }, { 0, 0, 0, 1 } };
          double[,] DataT4m3 = { { 1, 0, 0, 0 }, { 0, 1, 0, 0 }, { 0, 0, 1, -ls3 }, { 0, 0, 0, 1 } };
          double[,] DataT5mH = { { 1, 0, 0, -lsH }, { 0, 1, 0, 0 }, { 0, 0, 1, 0 }, { 0, 0, 0, 1 } };
          double[,] DataTFacePlatemP = { { 1, 0, 0, lsPx }, { 0, 1, 0, lsPy }, { 0, 0, 1, lsPz }, { 0, 0, 0, 1 } };

          Matrix WC2Robot = new Matrix(DataWC2Robot);
          Matrix Tb1 = new Matrix(DataTb1);
          Matrix T12 = new Matrix(DataT12);
          Matrix T23 = new Matrix(DataT23);
          Matrix T34 = new Matrix(DataT34);
          Matrix T45 = new Matrix(DataT45);
          Matrix T56 = new Matrix(DataT56);
          Matrix TFacePlate = new Matrix(DataTFacePlate);

          Matrix T2m2 = new Matrix(DataT2m2);
          Matrix T4m3 = new Matrix(DataT4m3);
          Matrix T5mH = new Matrix(DataT5mH);
          Matrix TFacePlatemP = new Matrix(DataTFacePlatemP);

          #endregion

          #region Payload-induced Deflection Calculations
          //Transformations from Joint(n) to TCP
          //Example: T2T is the transform from joint 2 to the TCP of the end effector
          Matrix A = WC2Robot * Tb1;
          Matrix B = A * T12;
          Matrix C = B * T23;
          Matrix D = C * T34;
          Matrix E = D * T45;
          Matrix F = E * T56;

          Matrix invA = Matrix.Inv(A);
          Matrix invB = Matrix.Inv(B);
          Matrix invC = Matrix.Inv(C);
          Matrix invD = Matrix.Inv(D);
          Matrix invE = Matrix.Inv(E);
          Matrix invF = Matrix.Inv(F);

          //Transformations and vectors from joints to individual masses
          Matrix TmPA6 = TFacePlate * TFacePlatemP;
          Vector rmPA6 = TmPA6.r();
          Matrix TmPA5 = T56 * TmPA6;
          Vector rmPA5 = TmPA5.r();
          Matrix TmPA4 = T45 * TmPA5;
          Vector rmPA4 = TmPA4.r();
          Matrix TmPA3 = T34 * TmPA4;
          Vector rmPA3 = TmPA3.r();
          Matrix TmPA2 = T23 * TmPA3;
          Vector rmPA2 = TmPA2.r();
          Matrix TmHA3 = T34 * T45 * T5mH;
          Vector rmHA3 = TmHA3.r();
          Matrix TmHA2 = T23 * TmHA3;
          Vector rmHA2 = TmHA2.r();
          Matrix Tm3A3 = T34 * T4m3;
          Vector rm3A3 = Tm3A3.r();
          Matrix Tm3A2 = T23 * Tm3A3;
          Vector rm3A2 = Tm3A2.r();
          Matrix Tm2A2 = T2m2;
          Vector rm2A2 = Tm2A2.r();

          Vector rA6 = (mP * rmPA6) / mA6;
          Vector rA5 = (mP * rmPA5) / mA5;
          Vector rA4 = (mP * rmPA4) / mA4;
          Vector rA3 = ((mP * rmPA3) + (mH * rmHA3) + (m3 * rm3A3)) / mA3;
          Vector rA2 = ((mP * rmPA2) + (mH * rmHA2) + (m3 * rm3A2) + (m2 * rm2A2)) / mA2;

          Vector MomentA1 = new Vector(0, 0, 0);	//R X F(unit vector of Z axis) tooltip relative to joint[x]
          Vector MomentA2 = -mA2 * (rA2.Cross(invB.j()));	//Result is the torque vector about each joint unduced by a unit clamp force
          Vector MomentA3 = -mA3 * (rA3.Cross(invC.j()));	//Clamp force is assumed to be in the direction of the Z axis of the end effector
          Vector MomentA4 = -mA4 * (rA4.Cross(invD.j()));
          Vector MomentA5 = -mA5 * (rA5.Cross(invE.j()));
          Vector MomentA6 = -mA6 * (rA6.Cross(invF.j()));

          JointMoments[0] = MomentA1.k();
          JointMoments[1] = MomentA2.k();
          JointMoments[2] = MomentA3.k();
          JointMoments[3] = MomentA4.k();
          JointMoments[4] = MomentA5.k();
          JointMoments[5] = MomentA6.k();
          #endregion

          #region Debug
          //System.Console.WriteLine(JointMoments[0]);
          //System.Console.WriteLine(JointMoments[1]);
          //System.Console.WriteLine(JointMoments[2]);
          //System.Console.WriteLine(JointMoments[3]);
          //System.Console.WriteLine(JointMoments[4]);
          //System.Console.WriteLine(JointMoments[5]);
          #endregion

        }
      }
      #endregion
    }
    #endregion
  }
  //namespace FileParser
  //{
  //  public class cFileParse : IDisposable
  //  {
  //    public Electroimpact.StringCalc.cStringCalc _StringCalc = new Electroimpact.StringCalc.cStringCalc();

  //    #region Publics
  //    public cFileParse()
  //    {
  //    }
  //    ~cFileParse()
  //    {
  //      Dispose(false);
  //    }
  //    public bool GetArgument(string szLine, string varname, out double dPos)
  //    {
  //      string sArg;
  //      int StartArg;
  //      int Length;
  //      return this.GetArgument(szLine, varname, out dPos, out sArg, out StartArg, out Length);
  //    }
  //    public bool GetArgument(string szLine, string varname, out double dPos, out string sArg, out int StartArg, out int Length)
  //    {
  //      csString iString = new csString();
  //      bool bTestLeft = false;
  //      bool bTestRight = false;
  //      int nStartArg = 0;
  //      int nEndArg = 0;
  //      dPos = 0;
  //      sArg = "0";
  //      StartArg = -1;
  //      Length = -1;
  //      while (bTestLeft == false || bTestRight == false)
  //      {
  //        nStartArg = szLine.IndexOf(varname, nStartArg, szLine.Length - nStartArg);
  //        if (szLine.Length <= nStartArg + varname.Length) nStartArg = -1;
  //        if (nStartArg == -1)
  //          break;
  //        if (nStartArg == 0)
  //          bTestLeft = true;
  //        else
  //          bTestLeft = this.CheckLeftSide(szLine.Substring(nStartArg - 1, 1));

  //        bTestRight = this.CheckRightSide(szLine.Substring(nStartArg + varname.Length, 1));

  //        nStartArg += varname.Length;
  //      }

  //      if (bTestRight && bTestLeft)
  //      {
  //        StartArg = nStartArg;
  //        iString.String = szLine.Substring(nStartArg, szLine.Length - nStartArg);
  //        nEndArg = nStartArg = 0;
  //        string ch = iString.GetLeft(1);
  //        bool gotEq = false;
  //        if (ch == "")
  //          return false;
  //        if (this.CheckForOpenParen(ch))
  //        {
  //          int LP = 1;
  //          while ((ch = iString.GetLeft(1)) != "")
  //          {
  //            nEndArg++;
  //            switch (ch)
  //            {
  //              case "{":
  //              case "(":
  //              case "[":
  //                LP++;
  //                break;
  //              case "}":
  //              case ")":
  //              case "]":
  //                LP--;
  //                break;
  //              default:
  //                break;
  //            }
  //            gotEq = LP == 0;
  //            if (gotEq)
  //              break;
  //          }
  //          if (gotEq)
  //          {
  //            string eq = iString.String.Substring(nStartArg, nEndArg + 1);
  //            dPos = this._StringCalc.SimpleCalc(eq);
  //            return dPos != double.NaN;
  //          }
  //        }
  //        else
  //        {
  //          iString.Rewind();
  //          string rest = iString.String;
  //          string szCompare = "0123456789+-."; //The '+' is questionable.

  //          for (int ii = 0; ii < rest.Length; ii++)
  //          {
  //            if (!iString.InString(rest.Substring(ii, 1), szCompare))
  //              break;
  //            nEndArg++;
  //          }
  //          string szResult = rest.Substring(0, nEndArg);
  //          dPos = iString.ToDouble(szResult);
  //          Length = nEndArg;
  //          sArg = szResult;
  //          return true;
  //        }
  //      }
  //      return false;
  //    }
  //    /// <summary>
  //    /// Replaces whatever argument held in varname with newval.  Returns true if it found the argument and changed it. 
  //    /// </summary>
  //    /// <param name="szLine"></param>
  //    /// <param name="varname"></param>
  //    /// <param name="newval"></param>
  //    /// <param name="NewLine"></param>
  //    /// <returns></returns>
  //    public bool ReplaceArgument(string szLine, string varname, double newval, out string NewLine)
  //    {
  //      double dhold;
  //      string sz;
  //      int Start, Length;
  //      NewLine = szLine;
  //      if (GetArgument(szLine, varname, out dhold, out sz, out Start, out Length))
  //      {
  //        string left = szLine.Substring(0, Start);
  //        string right = szLine.Substring(Start + Length);
  //        NewLine = left + newval.ToString("F3") + right;
  //        return true;
  //      }
  //      return false;
  //    }
  //    public bool ReplaceArgument(string szLine, string varname, string newval, out string NewLine)
  //    {
  //      double dhold;
  //      string sz;
  //      int Start, Length;
  //      NewLine = szLine;
  //      if (GetArgument(szLine, varname, out dhold, out sz, out Start, out Length))
  //      {
  //        string left = szLine.Substring(0, Start);
  //        string right = szLine.Substring(Start + Length);
  //        NewLine = left + newval + right;
  //        return true;
  //      }
  //      return false;
  //    }
  //    #endregion

  //    #region Privates
  //    private bool CheckLeftSide(string szIn)
  //    {
  //      //Watch out...I added a weird case "("
  //      string[] szOK = { "(", ".", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", " ", "]", "}", ")", "-", "+" };

  //      for (int ii = 0; ii < szOK.Length; ii++)
  //      {
  //        if (szIn == szOK[ii])
  //        {
  //          return true;
  //        }
  //      }
  //      return false;
  //    }

  //    private bool CheckRightSide(string szIn)
  //    {
  //      string[] szOK = { ".", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", " ", "[", "{", "(", "-", "+" };

  //      for (int ii = 0; ii < szOK.Length; ii++)
  //      {
  //        if (szIn == szOK[ii])
  //        {
  //          return true;
  //        }
  //      }
  //      return false;
  //    }

  //    private bool CheckForOpenParen(string szIn)
  //    {
  //      string[] szOK = { "(", "{", "[" };

  //      for (int ii = 0; ii < szOK.Length; ii++)
  //      {
  //        if (szIn == szOK[ii])
  //        {
  //          return true;
  //        }
  //      }
  //      return false;
  //    }

  //    #endregion

  //    #region IDisposable Members

  //    void IDisposable.Dispose()
  //    {
  //      Dispose(true);
  //      System.GC.SuppressFinalize(this);
  //    }

  //    void Dispose(bool explicitCall)
  //    {
  //      if (explicitCall)
  //      {
  //        if (_StringCalc != null)
  //          _StringCalc = null;
  //      }
  //    }

  //    #endregion
  //  }
  //}
}

