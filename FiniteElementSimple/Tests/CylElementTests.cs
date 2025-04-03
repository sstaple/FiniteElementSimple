using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Drawing;
using FiniteElementSimple.Elements;

namespace FiniteElementSimple.Tests
{
    class CylElementTests
    {

        public static void MapElements()
        {

            LinElastic2DPlaneStress myMaterial = new LinElastic2DPlaneStress(70000, 0.33);
            double thickness = 1.0;

			double[] fiberCenter = new double[] { 1d, 1d };

            int[][] ConnectivityMatrix = {new int[]{0,1,2,3}, new int[] { 1, 4, 5, 2 }, new int[] { 4, 6, 5 } };
            double[][] NodalLocations = new double[][] { new double[] { 1, 1 }, new double[] { 2, 1 }, new double[] { 1, 2 }, new double[] { 1, 1 }, 
				new double[] { 3, 1 }, new double[] { 1, 3}, new double[] { 3, 3 }};

            List<BC> lLoads = new List<BC>();
            List<BC> lBCs = new List<BC>();
            List<Element> lElements = new List<Element>();

            lElements.Add(new CylindricalAtOrigin_Node4Element2D(myMaterial, ConnectivityMatrix[0], thickness, fiberCenter,	NodalLocations));
			lElements.Add(new Cylindrical_Node4Element2D(myMaterial, ConnectivityMatrix[1], thickness, fiberCenter, NodalLocations));
			lElements.Add(new Cylindrical_Node3Element2D(myMaterial, ConnectivityMatrix[2], thickness, fiberCenter, NodalLocations));


			Assembly myAssembly = new Assembly(lElements, lLoads, lBCs, 2);

			myAssembly.PlotOutline(50);

			CreateContourPlotSquare(myAssembly, 10, 10, 4, 3, false);
        }
		public static void MapTwoTriangleElements()
		{

			LinElastic2DPlaneStress myMaterial = new LinElastic2DPlaneStress(70000, 0.33);
			double thickness = 1.0;

			double[] fiberCenter = new double[] { 1d, 1d };

			int[][] ConnectivityMatrix = { new int[] { 1, 2, 0 }, new int[] { 2, 1, 3 } };
			double[][] NodalLocations = new double[][] { new double[] { 1, 1 }, new double[] { 2, 1 }, new double[] { 1, 2 }, new double[] { 3, 2 } };

			List<BC> lLoads = new List<BC>();
			List<BC> lBCs = new List<BC>();
			List<Element> lElements = new List<Element>();

			lElements.Add(new Cylindrical_Node3Element2D(myMaterial, ConnectivityMatrix[0], thickness, fiberCenter,
				NodalLocations));

			lElements.Add(new Cylindrical_Node3Element2D(myMaterial, ConnectivityMatrix[1], thickness, fiberCenter, NodalLocations));

			Assembly myAssembly = new Assembly(lElements, lLoads, lBCs, 2);

			CreateContourPlotTriangle(myAssembly, 25, 25, 3, 1, false);
		}

		public static void CreateContourPlotSquare(Assembly myAssembly, int nDivPerElX, int nDivPerElY, int u0_e1_s2_j3_n4, int direction, bool plotDeformedState)
		{

			//Initialize lists to hold data: x and y are coordinates, z is the quantity to be plotted (displacement, strain, stress)
			//Lists are arrays that can be resized dynamically
			List<double> lYData = new List<double>();
			List<double> lXData = new List<double>();
			List<double> lZData = new List<double>();
			string zTitle = "";

			//Compile the data

			//Loop through each element
			for (int i = 0; i < myAssembly.lElements.Count; i++)
			{

				//Loop through the number of divisions within the element in the x direction
				for (int j = 0; j < nDivPerElX + 1; j++)
				{

					//Set the xi value just by incrementing
					double tempXi = -1.0 + 2.0 / (nDivPerElX) * j;

					//Loop through the number of divisions within the element in the y direction
					for (int k = 0; k < nDivPerElY + 1; k++)
					{

						//Set the eta value just by incrementing
						double tempEta = -1.0 + 2.0 / (nDivPerElY) * k;

						//Get the global x and y value from the local xi and eta
						double[] tempX = myAssembly.lElements[i].GlobalXPosition(tempXi, tempEta, 0.0);
						double tempZ = 0;

						//Extrapolate the desired output (defined by the user in the function inputs
						switch (u0_e1_s2_j3_n4)
						{
							case 0:
								tempZ = myAssembly.lElements[i].Displacement(tempXi, tempEta, 0.0)[direction];
								zTitle = "Displacement, U_" + (direction + 1);
								break;
							case 1:
								tempZ = myAssembly.lElements[i].Strain(tempXi, tempEta, 0.0)[direction];
								zTitle = "Strain, E_" + (direction + 1);
								break;
							case 2:
								tempZ = myAssembly.lElements[i].Stress(tempXi, tempEta, 0.0)[direction];
								zTitle = "Stress, S_" + (direction + 1);
								break;
							case 3:
								tempZ = myAssembly.lElements[i].J(tempXi, tempEta, 0.0)[direction, direction];
								zTitle = "Jacobian, J_" + (direction + 1) + "," + (direction + 1);
								break;
							case 4:
								tempZ = myAssembly.lElements[i].ShapeFunction(tempXi, tempEta, 0.0)[0, direction * 2];
								zTitle = "Shape Function_" + (direction);
								break;
						}
						//If the points should be in the deformed shape, then add the displacements to the coordinates
						if (plotDeformedState)
						{
							double[] def = myAssembly.lElements[i].Displacement(tempXi, tempEta, 0.0);
							tempX[0] += def[0];
							tempX[1] += def[1];
						}
						//Add the x, y, and z values into their lists
						lXData.Add(tempX[0]);
						lYData.Add(tempX[1]);
						lZData.Add(tempZ);
					}
				}
			}
			//Then I make a color scheme and create the plot.  This part will be different for you: I was calling a function I 
			//Had already written
			Color[] colorScheme = new Color[] { Color.Blue, Color.Aqua, Color.LimeGreen, Color.Yellow, Color.Red };
			SinglePlot.SinglePlotForm myPlot = new SinglePlot.SinglePlotForm(zTitle, lXData.ToArray(), lYData.ToArray(),
																			 lZData.ToArray(), colorScheme);
			myPlot.Activate();
			myPlot.ShowDialog();
		}

		public static void CreateContourPlotTriangle(Assembly myAssembly, int nDivPerElX, int nDivPerElY, int u0_e1_s2_j3, int direction, bool plotDeformedState)
		{

			//Initialize lists to hold data: x and y are coordinates, z is the quantity to be plotted (displacement, strain, stress)
			//Lists are arrays that can be resized dynamically
			List<double> lYData = new List<double>();
			List<double> lXData = new List<double>();
			List<double> lZData = new List<double>();
			string zTitle = "";

			//Compile the data

			//Loop through each element
			for (int i = 0; i < myAssembly.lElements.Count; i++)
			{

				//Loop through the number of divisions within the element in the x direction
				for (int j = 0; j < nDivPerElX + 1; j++)
				{

					//Set the xi value just by incrementing
					double tempXi = 1.0 / (nDivPerElX) * j;

					//Loop through the number of divisions within the element in the y direction
					for (int k = 0; k < nDivPerElY + 1; k++)
					{

						//Set the eta value just by incrementing
						double tempEta = (1.0 - tempXi) / (nDivPerElY) * k;

						//Get the global x and y value from the local xi and eta
						double[] tempX = myAssembly.lElements[i].GlobalXPosition(tempXi, tempEta, 0.0);
						double tempZ = 0;

						//Extrapolate the desired output (defined by the user in the function inputs
						switch (u0_e1_s2_j3)
						{
							case 0:
								tempZ = myAssembly.lElements[i].Displacement(tempXi, tempEta, 0.0)[direction];
								zTitle = "Displacement, U_" + (direction + 1);
								break;
							case 1:
								tempZ = myAssembly.lElements[i].Strain(tempXi, tempEta, 0.0)[direction];
								zTitle = "Strain, E_" + (direction + 1);
								break;
							case 2:
								tempZ = myAssembly.lElements[i].Stress(tempXi, tempEta, 0.0)[direction];
								zTitle = "Stress, S_" + (direction + 1);
								break;
							case 3:
								tempZ = myAssembly.lElements[i].J(tempXi, tempEta, 0.0)[direction, direction];
								zTitle = "Jacobian, J_" + (direction + 1) + "," + (direction + 1);
								break;
						}
						//If the points should be in the deformed shape, then add the displacements to the coordinates
						if (plotDeformedState)
						{
							double[] def = myAssembly.lElements[i].Displacement(tempXi, tempEta, 0.0);
							tempX[0] += def[0];
							tempX[1] += def[1];
						}
						//Add the x, y, and z values into their lists
						lXData.Add(tempX[0]);
						lYData.Add(tempX[1]);
						lZData.Add(tempZ);
					}
				}
			}
			//Then I make a color scheme and create the plot.  This part will be different for you: I was calling a function I 
			//Had already written
			Color[] colorScheme = new Color[] { Color.Blue, Color.Aqua, Color.LimeGreen, Color.Yellow, Color.Red };
			SinglePlot.SinglePlotForm myPlot = new SinglePlot.SinglePlotForm(zTitle, lXData.ToArray(), lYData.ToArray(),
																			 lZData.ToArray(), colorScheme);
			myPlot.Activate();
			myPlot.ShowDialog();
		}
	}
}
