
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
	/// Description of Quadratic8NodeElement.
	/// </summary>
	public class Cylindrical_Node4Element2D : Node4Element2D
	{
		#region private Properties
		private double[] center;

		#endregion
		/// <summary>
		/// 
		/// </summary>
		/// <param name="elementMaterial"></param>
		/// <param name="localToGlobalConnectivity_1index"></param>
		/// <param name="thickness"></param>
		/// <param name="center"></param>
		/// <param name="nodalLocations">Nodal locations need to be in pairs, with the first two being the lower r, the next two the higher r. Nodes numbered counter-clockwise</param>
		public Cylindrical_Node4Element2D(Material elementMaterial, int[] localToGlobalConnectivity_1index,
									 double thickness, double[] center, double[][] nodalLocations)
			: base(elementMaterial, localToGlobalConnectivity_1index, thickness, nodalLocations)
		{
			this.center = center;

			//Since the nodal coordinates have been given to the base, one can assume that the local nodal locations are in the global coordinate system right now
			//So we need to put them into the local

			//Vectors from the center to the corners
			double[] v1 = VectorMath.Subtract(nodalLocations[localToGlobalConnectivity_1index[0]], center);
			double r1 = VectorMath.Norm(v1);
			double[] v2 = VectorMath.Subtract(nodalLocations[localToGlobalConnectivity_1index[1]], center);
			double r2 = VectorMath.Norm(v2);
			double[] v3 = VectorMath.Subtract(nodalLocations[localToGlobalConnectivity_1index[2]], center);
			double r3 = VectorMath.Norm(v3);
			double[] v4 = VectorMath.Subtract(nodalLocations[localToGlobalConnectivity_1index[3]], center);
			double r4 = VectorMath.Norm(v4);


			//loacl cylindrical coordinate system

			//convert the nodal locations into r, in the form of r1, theta1, ....
			this.nodalLocations = new double[8];
			//Node 0: Assume that r=0 is at x=1, y=0
			this.nodalLocations[0] = r1;
			this.nodalLocations[1] = Cylindrical_Node3Element2D.CalculateRAndTheta(r1, v1);
			//Node 1: on matrix, lower corner
			this.nodalLocations[2] = r2;
			this.nodalLocations[3] = Cylindrical_Node3Element2D.CalculateRAndTheta(r2, v2);
			//Node 2: on matrix, upper corner
			this.nodalLocations[4] = r3;
			this.nodalLocations[5] = Cylindrical_Node3Element2D.CalculateRAndTheta(r3, v3);

			this.nodalLocations[6] = r4;
			this.nodalLocations[7] = Cylindrical_Node3Element2D.CalculateRAndTheta(r4, v4);
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

			double[,] A = Cylindrical_Node3Element2D.AMatrix(r);
			double[,] A0 = Cylindrical_Node3Element2D.A0Matrix(r);

			double[,] AJinvhat = MatrixMath.Multiply(A, Jinvhat);
			double[,] AJinvhatdNhat = MatrixMath.Multiply(AJinvhat, dNhat);

			double[,] A0N = MatrixMath.Multiply(A0, N);

			return MatrixMath.Add(AJinvhatdNhat, A0N);
		}

		public override double Det_Of_J(double xi, double eta, double zeta)
		{
			//Slide the xi in there because the indegral of the area is r*dr*dtheta
			double[,] myJ = J(xi, eta, zeta);
			return MatrixMath.Determinant(myJ) * xi;
		}

		public override double[] GlobalXPosition(double xi, double eta, double zeta)
		{
			double[,] N = ShapeFunction(xi, eta, zeta);
			double[] x = MatrixMath.Multiply(N, nodalLocations);
			return Cylindrical_Node3Element2D.ConvertCylindricalToCartesian(x[0], x[1], 0, center);
		}
	}
}

