/*
 * Created by SharpDevelop.
 * User: scott_stapleton
 * Date: 4/18/2019
 * Time: 6:13 PM
 * 
 * To change this template use Tools | Options | Coding | Edit Standard Headers.
 */
using System;
using RandomMath;

namespace FiniteElementSimple.Elements
{
	/// <summary>
	/// Description of MatrixCylindricalElement.
	/// </summary>
	public class Cylindrical_Node3Element2D : TriangularElement2D
	{
		#region private Properties
		double[] center;

		#endregion

		#region Public Properties

		#endregion


		/// <summary>
		/// element is a triangle, with Node 1 at 0,0, Node 2 at 1,0, and Node 3 at 0,1 (coordinates are in the form of: (r, theta)).  0 and 1 should be the same r.
		/// </summary>
		/// <param name="elementMaterial"></param>
		/// <param name="localToGlobalConnectivity_1index"></param>
		/// <param name="thickness"></param>
		/// <param name="center"></param>
		/// <param name="nodalLocations">This should have nearly the same radius as corner 2</param>
		public Cylindrical_Node3Element2D(Material elementMaterial, int[] localToGlobalConnectivity_1index,
									 double thickness, double[] center, double[][] nodalLocations)
			: base(elementMaterial, localToGlobalConnectivity_1index, thickness, nodalLocations, 3) //TODO: check if this works for all triangular elements
		{
			this.center = center;
			//Since the nodal coordinates have been given to the base, one can assume that the local nodal locations are in the global coordinate system right now
			//So we need to put them into the local cylindrical system

			//Vectors from the center to the corners
			double[] v1 = VectorMath.Subtract(nodalLocations[localToGlobalConnectivity_1index[0]], center);
			double r1 = VectorMath.Norm(v1);
			double[] v2 = VectorMath.Subtract(nodalLocations[localToGlobalConnectivity_1index[1]], center);
			double r2 = VectorMath.Norm(v2);
			double[] v3 = VectorMath.Subtract(nodalLocations[localToGlobalConnectivity_1index[2]], center);
			double r3 = VectorMath.Norm(v3);


			//loacl cylindrical coordinate system

			//convert the nodal locations into r, in the form of r1, theta1, ....
			this.nodalLocations = new double[6];
			//Node 0: Assume that r=0 is at x=1, y=0
			this.nodalLocations[0] = r1;
			this.nodalLocations[1] = CalculateRAndTheta(r1, v1);
			//Node 1: on matrix, lower corner
			this.nodalLocations[2] = r2;
			this.nodalLocations[3] = CalculateRAndTheta(r2, v2);
			//Node 2: on matrix, upper corner
			this.nodalLocations[4] = r3;
			this.nodalLocations[5] = CalculateRAndTheta(r3, v3);

		}

		public static double CalculateRAndTheta(double r, double[] v)
        {
			//If r=0 (node is at the origin), make the angle 0
            if (r.Equals(0.0))
            {
				//return 0.0;
				return Math.PI / 2d;
            }
			//This expression is to get the angle from the x axis
			return Math.Acos(v[0] / r);
        }

		public override double[,] ShapeFunction(double xi, double eta, double zeta)
		{
			double r = xi;
			double theta = eta;
			//double N1, N2, N3, N4, N5;

			double N1 = 1.0 - r - theta;
			double N2 = r;
			double N3 = theta;
			return new double[,]{{N1, 0.0, N2, 0.0, N3, 0.0},
				{0.0, N1, 0.0, N2, 0.0, N3}};
		}

		public override double[,] DNdxi(double xi, double eta, double zeta)
		{
			double dN1dxi = -1.0;
			double dN1deta = -1.0;

			double dN2dxi = -1.0;
			double dN2deta = 0.0;

			double dN3dxi = 0.0;
			double dN3deta = -1.0;


			return new double[,]{{dN1dxi, 0.0, dN2dxi, 0.0, dN3dxi, 0.0},
				{dN1deta, 0.0, dN2deta, 0.0, dN3deta, 0.0},
				{0.0, dN1dxi, 0.0, dN2dxi, 0.0, dN3dxi},
			{0.0, dN1deta, 0.0, dN2deta, 0.0, dN3deta}};
		}

		public override double[,] B(double xi, double eta, double zeta)
		{
			//Change B here for the Cylindrical Coordinate System
			double r = xi;

			double[,] N = ShapeFunction(xi, eta, zeta);
			double[,] dNhat = DNdxi(xi, eta, zeta);
			double[,] Jinv = MatrixMath.InvertMatrix(J(xi, eta, zeta));

			double[,] Jinvhat = new double[4, 4];
			MatrixMath.CopyToMatrix(ref Jinvhat, Jinv, 0, 0);
			MatrixMath.CopyToMatrix(ref Jinvhat, Jinv, 2, 2);

			double[,] A = AMatrix(r);
			double[,] A0 = A0Matrix(r);

			double[,] AJinvhat = MatrixMath.Multiply(A, Jinvhat);
			double[,] AJinvhatdNhat = MatrixMath.Multiply(AJinvhat, dNhat);

			double[,] A0N = MatrixMath.Multiply(A0, N);

			return MatrixMath.Add(AJinvhatdNhat, A0N);
		}

		public override double[] GlobalXPosition(double xi, double eta, double zeta)
		{
			double[,] N = ShapeFunction(xi, eta, zeta);
			double[] x = MatrixMath.Multiply(N, nodalLocations);
			return ConvertCylindricalToCartesian(x[0], x[1], 0, center);
		}

		public static double[] ConvertCylindricalToCartesian(double r, double theta, double z, double[] origin)
		{
			double x = r * Math.Cos(theta) + origin[0];
			double y = r * Math.Sin(theta) + origin[1];
			return new double[] { x, y, z };
		}

		public static double[] ConvertCartesianToCylindrical(double x, double y, double z)
		{
			double r = Math.Sqrt(x * x + y * y + z * z);
			double theta = Math.Atan2(x, y);
			return new double[] { r, theta, z };
		}

		public static double[,] AMatrix(double r)
		{
			return new double[,]{{1.0, 0.0, 0.0, 0.0},
			{0.0, 0.0, 0.0, 1.0 / r},
			{0.0, 1.0 / (2.0 * r),  1.0 / (2.0 * r), 0.0}};
		}

		public static double[,] A0Matrix(double r)
		{
			return new double[,]{{0.0, 0.0},
			{1.0 / r, 0.0},
			{0.0, -1.0 / (2.0 * r)}};
		}
	}
}



