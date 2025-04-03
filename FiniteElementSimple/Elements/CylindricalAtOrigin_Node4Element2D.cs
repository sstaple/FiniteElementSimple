using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FiniteElementSimple.Elements
{
    class CylindricalAtOrigin_Node4Element2D : Cylindrical_Node4Element2D
    {
        /// <summary>
        /// This collapses nodes 0 and 3 to the same point.  Nodes 0 and 3 must be at the origin
        /// </summary>
        /// <param name="elementMaterial"></param>
        /// <param name="localToGlobalConnectivity_1index"></param>
        /// <param name="thickness"></param>
        /// <param name="center"></param>
        /// <param name="nodalLocations"></param>
        public CylindricalAtOrigin_Node4Element2D(Material elementMaterial, int[] localToGlobalConnectivity_1index,
                                     double thickness, double[] center, double[][] nodalLocations)
            : base(elementMaterial, localToGlobalConnectivity_1index, thickness, center, nodalLocations)
        {
            if (!this.nodalLocations[0].Equals(0.0) || !this.nodalLocations[6].Equals(0.0))
            {
                throw new Exception("For a 4-Noded Cylindrical Element at the origin, the 0 and 1 nodes must be at the origin");
            }

            //Now, make the angles of the 0 and 1 nodes the same as 3 and 2 nodes respectively
            this.nodalLocations[1] = this.nodalLocations[3];
            this.nodalLocations[7] = this.nodalLocations[5];
        }
    }
}
