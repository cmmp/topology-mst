package br.fapesp.topology;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;

import org.apache.commons.math3.stat.regression.SimpleRegression;

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
	
	private static final Random RNG = new Random(1234);
	
	/**
	 * Runs a topological analysis on the data, returning
	 * properties of its connectivity.
	 * 
	 * @param data a finite point-set cloud
	 * @param mineps the minimum eps value for analysis
	 * @param maxeps the maximum eps value for analysis
	 * @param nsteps the number of elements from mineps to maxeps: 
	 * 		20 is a good number usually. Be advised the generated series will 
	 * 		present exponential growth, so it becomes linear when log is taken
	 * 		to compute coefficients.
	 * @return a matrix with the following columns: eps, Ce, De, Ie, Coefs
	 * The first row of the coefficients column contains the gamma coefficient, and the
	 * second row contains the delta coefficient.
	 * 
	 */
	public static double[][] analysis(double[][] data, double mineps, double maxeps, int nsteps) {
		// create an exponential series because when its log is computed
		// we get linearly spaced points.
		//double[] epsValues = MyUtils.genExpSeries(2, nsteps);
		double[] epsValues = MyUtils.linspace(mineps, maxeps, nsteps);
		
		int Nx = epsValues.length;
		
		// results has the following columns:
		// x, Ce, De, Ie, Coefs
		// the Coefs column has only two elements:
		// in the first row, the gamma coefficient and, in the second row,
		// the delta coefficient
		double[][] results = new double[Nx][5];
		
		// put the exponential series in the desired range:
		//epsValues = MyUtils.rescaleToRange(epsValues, mineps, maxeps);
		
		double[][] D = MyUtils.getEuclideanMatrix(data);
		int[][] mst = MyUtils.fastPrim(data, D);
		ArrayList<Edge> cmst = new ArrayList<Edge>();
		int[][] cmstar;
		
		int i = 0;
		int N = data.length;
		ArrayDeque<Integer> points = new ArrayDeque<Integer>(N);
		HashSet<Integer> epscomp = new HashSet<Integer>();
		ArrayList<Integer> epscompList;
		ArrayDeque<Integer> expandSet = new ArrayDeque<Integer>();
		int ce, ie;
		int u, v;
		double de;
		
		for (double eps : epsValues) {
			ce = 0; de = 0; ie = 0;
			
			results[i][0] = eps;
			
			// compute C(e) by counting the number of eps-connected components, 
			// which is just one more than the number of edges with length 
			// greater than eps
			
			for (int j = 0; j < mst.length; j++)
				if (D[mst[j][0]][mst[j][1]] <= eps)
					cmst.add(new Edge(mst[j][0],mst[j][1],D[mst[j][0]][mst[j][1]]));
				else
					ce++;
			
			cmstar = new int[cmst.size()][2];
			for (int j = 0; j < cmst.size(); j++) {
				cmstar[j][0] = cmst.get(j).u;
				cmstar[j][1] = cmst.get(j).v;
			}
			
			// create an adjacency list from the cmst:
			HashMap<Integer, ArrayList<Integer>> adjacency = new HashMap<Integer, ArrayList<Integer>>();
			ArrayList<Integer> neighbors;
			for (int j = 0; j < cmstar.length; j++) {
				if (!adjacency.containsKey(cmstar[j][0])) {
					neighbors = new ArrayList<Integer>();
					neighbors.add(cmstar[j][1]);
					adjacency.put(cmstar[j][0], neighbors);
				} else {
					neighbors = adjacency.get(cmstar[j][0]);
					neighbors.add(cmstar[j][1]);
					adjacency.remove(cmstar[j][0]);
					adjacency.put(cmstar[j][0], neighbors);
				}
				if (!adjacency.containsKey(cmstar[j][1])) {
					neighbors = new ArrayList<Integer>();
					neighbors.add(cmstar[j][0]);
					adjacency.put(cmstar[j][1], neighbors);
				} else {
					neighbors = adjacency.get(cmstar[j][1]);
					neighbors.add(cmstar[j][0]);
					adjacency.remove(cmstar[j][1]);
					adjacency.put(cmstar[j][1], neighbors);
				}
			}
			
			// compute D(e) by expanding each eps-connected component and finding the
			// maximum distance among all points there.
			
			// compute I(e) by counting the number of isolated eps-connected components
			
			points.clear();
			for (int j = 0; j < N; j++) 
				points.add(j);
			
			de = 0;
			boolean visited[] = new boolean[N];
			
			while (points.size() > 0) { // expand eps-components
				epscomp.clear();
				expandSet.clear();
				
				u = points.pollFirst();
				if (visited[u]) continue;
				
				epscomp.add(u);
				expandSet.add(u);
				
				while (expandSet.size() > 0) {
					v = expandSet.pollFirst();
					visited[v] = true;
					if (adjacency.get(v) != null)
						for (int k : adjacency.get(v)) {
							if (!epscomp.contains(k)) {
								visited[k] = true;
								epscomp.add(k);
								expandSet.add(k);
							}
						}
				}
				
				double dist;
				epscompList = new ArrayList<Integer>(epscomp);
				// now find out if this eps-component has
				// a diameter greater than de:
				for (int k = 0; k < epscompList.size() - 1; k++)
					for (int l = k + 1; l < epscompList.size(); l++) {
						dist = D[epscompList.get(k)][epscompList.get(l)];
						if (dist > de)
							de = dist;
					}
						
				
				if (epscomp.size() == 1)
					ie++;
			}
			
			results[i][1] = ++ce;
			results[i][2] = de;
			results[i][3] = ie;
			
			i++;
		}
		
		// compute gamma and delta coefficients:
		SimpleRegression gammaReg = new SimpleRegression();
		SimpleRegression deltaReg = new SimpleRegression();
		
		for (int j = 0; j < results.length; j++) {
			if (results[j][3] < 0.9 * N) {
				gammaReg.addData(Math.log(1.0 / results[j][0]), Math.log(results[j][1]));
				deltaReg.addData(Math.log(results[j][0]), Math.log(results[j][2]));
			}
		}
		
		results[0][4] = gammaReg.getSlope();
		results[1][4] = deltaReg.getSlope();
		
		return results;
	}
	
	private static double getMinGreaterThanZero(double[][] matrix) {
		double min = matrix[0][1];
		for (int i = 0; i < matrix.length; i++)
			for (int j = 0; j < matrix[0].length; j++)
				if(matrix[i][j] < min && matrix[i][j] > 0)
					min = matrix[i][j];
		return min;
	}
	
	/**
	 * Runs a topological analysis on the data, returning
	 * properties of its connectivity.
	 * 
	 * Mineps and maxeps are determined automatically in this version
	 * via the dissimilarity matrix.
	 * 
	 * @param data a finite point-set cloud
	 * @param nsteps the number of elements from mineps to maxeps: 
	 * 		20 is a good number usually.
	 * @return a matrix with the following columns: eps, Ce, De, Ie, Coefs
	 * The first row of the coefficients column contains the gamma coefficient, and the
	 * second row contains the delta coefficient.
	 * 
	 */
	public static double[][] analysis(double[][] data, int nsteps) {
		// create an exponential series because when its log is computed
		// we get linearly spaced points.
		//double[] epsValues = MyUtils.genExpSeries(2, nsteps);
		
		// put the exponential series in the desired range:
		//epsValues = MyUtils.rescaleToRange(epsValues, mineps, maxeps);
		
//		if (data.length < 10) { // we need to synthesize more data, otherwise estimation is very poor
//			int PADDING = 20; // how much new data we synthesize
//			double[][] newdata = new double[data.length + PADDING][data[0].length];
//			for(int i = 0; i < data.length; i++)
//				System.arraycopy(data[i], 0, newdata[i], 0, data[0].length);
//			// synthesize new data adding a small jitter
//			for(int i = data.length; i < newdata.length; i++) {
//				// pick a random object from the original data
//				int pick = RNG.nextInt(data.length);
//				for (int j = 0; j < data[0].length; j++)
//					newdata[i][j] = data[pick][j] + RNG.nextGaussian() / 100.0; // add a N(0,0.25) value
//			}
//			data = newdata;
//		}
		
		double[][] D = MyUtils.getEuclideanMatrix(data);
		//double[] epsValues = MyUtils.linspace(getMinGreaterThanZero(D), MyUtils.getMatrixMax(D), nsteps);
		double[] epsValues = MyUtils.expspace(getMinGreaterThanZero(D), MyUtils.getMatrixMax(D), nsteps);
		
		int Nx = epsValues.length;
		
		// results has the following columns:
		// x, Ce, De, Ie, Coefs
		// the Coefs column has only two elements:
		// in the first row, the gamma coefficient and, in the second row,
		// the delta coefficient
		double[][] results = new double[Nx][5];
				
		int[][] mst = MyUtils.fastPrim(data, D);
		ArrayList<Edge> cmst = new ArrayList<Edge>();
		int[][] cmstar;
		
		int i = 0;
		int N = data.length;
		ArrayDeque<Integer> points = new ArrayDeque<Integer>(N);
		HashSet<Integer> epscomp = new HashSet<Integer>();
		ArrayList<Integer> epscompList;
		ArrayDeque<Integer> expandSet = new ArrayDeque<Integer>();
		int ce, ie;
		int u, v;
		double de;
		
		for (double eps : epsValues) {
			ce = 0; de = 0; ie = 0;
			
			results[i][0] = eps;
			
			// compute C(e) by counting the number of eps-connected components, 
			// which is just one more than the number of edges with length 
			// greater than eps
			
			for (int j = 0; j < mst.length; j++)
				if (D[mst[j][0]][mst[j][1]] <= eps)
					cmst.add(new Edge(mst[j][0],mst[j][1],D[mst[j][0]][mst[j][1]]));
				else
					ce++;
			
			cmstar = new int[cmst.size()][2];
			for (int j = 0; j < cmst.size(); j++) {
				cmstar[j][0] = cmst.get(j).u;
				cmstar[j][1] = cmst.get(j).v;
			}
			
			// create an adjacency list from the cmst:
			HashMap<Integer, ArrayList<Integer>> adjacency = new HashMap<Integer, ArrayList<Integer>>();
			ArrayList<Integer> neighbors;
			for (int j = 0; j < cmstar.length; j++) {
				if (!adjacency.containsKey(cmstar[j][0])) {
					neighbors = new ArrayList<Integer>();
					neighbors.add(cmstar[j][1]);
					adjacency.put(cmstar[j][0], neighbors);
				} else {
					neighbors = adjacency.get(cmstar[j][0]);
					neighbors.add(cmstar[j][1]);
					adjacency.remove(cmstar[j][0]);
					adjacency.put(cmstar[j][0], neighbors);
				}
				if (!adjacency.containsKey(cmstar[j][1])) {
					neighbors = new ArrayList<Integer>();
					neighbors.add(cmstar[j][0]);
					adjacency.put(cmstar[j][1], neighbors);
				} else {
					neighbors = adjacency.get(cmstar[j][1]);
					neighbors.add(cmstar[j][0]);
					adjacency.remove(cmstar[j][1]);
					adjacency.put(cmstar[j][1], neighbors);
				}
			}
			
			// compute D(e) by expanding each eps-connected component and finding the
			// maximum distance among all points there.
			
			// compute I(e) by counting the number of isolated eps-connected components
			
			points.clear();
			for (int j = 0; j < N; j++) 
				points.add(j);
			
			de = 0;
			boolean visited[] = new boolean[N];
			
			while (points.size() > 0) { // expand eps-components
				epscomp.clear();
				expandSet.clear();
				
				u = points.pollFirst();
				if (visited[u]) continue;
				
				epscomp.add(u);
				expandSet.add(u);
				
				while (expandSet.size() > 0) {
					v = expandSet.pollFirst();
					visited[v] = true;
					if (adjacency.get(v) != null)
						for (int k : adjacency.get(v)) {
							if (!epscomp.contains(k)) {
								visited[k] = true;
								epscomp.add(k);
								expandSet.add(k);
							}
						}
				}
				
				double dist;
				epscompList = new ArrayList<Integer>(epscomp);
				// now find out if this eps-component has
				// a diameter greater than de:
				for (int k = 0; k < epscompList.size() - 1; k++)
					for (int l = k + 1; l < epscompList.size(); l++) {
						dist = D[epscompList.get(k)][epscompList.get(l)];
						if (dist > de)
							de = dist;
					}
						
				
				if (epscomp.size() == 1)
					ie++;
			}
			
			results[i][1] = ++ce;
			results[i][2] = de;
			results[i][3] = ie;
			
			i++;
		}
		
		// compute gamma and delta coefficients:
		SimpleRegression gammaReg = new SimpleRegression();
		SimpleRegression deltaReg = new SimpleRegression();
		
		// only add to the regression while the number of connected components
		// is changing. Once it reaches 1.0, stop adding more points.
		boolean reachedOneConnectedComponent = false;
		
		for (int j = 0; j < results.length; j++) {
			if(reachedOneConnectedComponent)
				break;
			if(results[j][1] == 1.0)
				reachedOneConnectedComponent = true;
				gammaReg.addData(results[j][0], results[j][1]);
				deltaReg.addData(results[j][0], results[j][2]);
		}
		
		results[0][4] = gammaReg.getSlope();
		results[1][4] = deltaReg.getSlope();
		
		return results;
	}
	
	public static void main(String[] args) {
		if (args.length != 1) {
			System.out.println("Usage: TopologyAnalysis input.dat");
			System.exit(0);
		}
		
		double[][] data = MyUtils.readCSVdataSet(args[0], false, ' ');
		
		double[][] results = TopologyAnalysis.analysis(data, 1e-5, 1.0, 20);
		
		MyUtils.print_matrix(results);
	}

}
