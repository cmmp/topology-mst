package br.fapesp.topology;

import java.util.ArrayList;

import br.fapesp.myutils.MyUtils;

/**
 * Topological analysis based on the PhD thesis:
 * Robins, V. (2000). Computational Topology at Multiple Resolutions: 
 * Foundations and Applications to Fractals and Dynamics.
 * 
 * @author Cássio M. M. Pereira
 *
 */
public class TopologyAnalysis {
	
	/**
	 * Runs a topological analysis on the data, returning
	 * properties of its connectivity.
	 * 
	 * @param data a finite point-set cloud
	 * @return a matrix with the following columns: x, log(x), Ce, log(Ce), De, log(De), Ie, log(Ie), Coefs
	 * The first row of the coefficients column contains the gamma coefficient, and the
	 * second row contains the delta coefficient.
	 * 
	 */
	public double[][] analysis(double[][] data) {
		// create an exponential series because when its log is computed
		// we get linearly spaced points.
		double[] epsValues = MyUtils.genExpSeries(2, 20);
		
		int Nx = epsValues.length;
		
		// results has the following columns:
		// x, log(x), Ce, log(Ce), De, log(De), Ie, log(Ie), Coefs
		// the Coefs column has only two elements:
		// in the first row, the gamma coefficient and, in the second row,
		// the delta coefficient
		double[][] results = new double[Nx][9];
		
		// put the exponential series in the desired range:
		epsValues = MyUtils.rescaleToRange(epsValues, 1e-5, 1.0);
		
		double[][] D = MyUtils.getEuclideanMatrix(data);
		int[][] mst = MyUtils.fastPrim(data, D);
		
		int i = 0;
		ArrayList<Edge> cmst;
		
		for (double eps : epsValues) {
			
			// create a copy of the MST:
			cmst = new ArrayList<Edge>();
			for (int j = 0; j < mst.length; j++)
				if (D[mst[j][0]][mst[j][1]] <= eps)
					cmst.add(new Edge(mst[j][0], mst[j][1], D[mst[j][0]][mst[j][1]]));
			
			// compute C(e) by counting the number of eps-connected components, 
			// which is just one more than the number of edges with length 
			// greater than eps
			results[i][0] = cmst.size() + 1;
			results[i][1] = Math.log(results[i][0]);
			
			i++;
		}
		
		return results;
	}
	
	public static void main(String[] args) {
		
	}

}
