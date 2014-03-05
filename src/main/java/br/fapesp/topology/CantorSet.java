package br.fapesp.topology;

import java.util.Random;

import br.fapesp.myutils.MyUtils;

public class CantorSet {
	
	private static double[] f1_26(double[] x) {
		return new double[] {0.5 * (-x[1]+1), 0.5 * x[0] };
	}
	
	private static double[] f2_26(double[] x) {
		return new double[] {0.5 * (x[1]+1), 0.5 * x[0] };
	}
	
	private static double[] f3_26(double[] x) {
		return new double[] {0.5 * x[1], 0.5 * (-x[0] + 2) };
	}
	
	
	/**
	 * Generates a Cantor set following equation 2.6 (figure 2.12) of Vanessa Robins's
	 * PhD thesis: "Computational Topology at Multiple Resolutions: Foundations and
	 * Applications to Fractals and Dynamics".
	 * 
	 * @param N the number of points to generate
	 * @param seed the seed for the RNG
	 * @return a 2-d cantor set data set
	 */
	public static double[][] genCantorSet26(int N, long seed) {
		double[][] x = new double[N][2];
		Random r = new Random(seed); 
		int sel;
		
		x[0][0] = 0; x[0][1] = 0;
		
		for(int i = 1; i < N; i++) {
			sel = r.nextInt(3);
			switch(sel) {
			case 0:
				x[i] = f1_26(x[i-1]);
				break;
			case 1:
				x[i] = f2_26(x[i-1]);
				break;
			case 2:
				x[i] = f3_26(x[i-1]);
				break;
			}
		}
		
		return x;
	}
	
	public static void main(String[] args) {
		if (args.length != 2) {
			System.out.println("Usage: CantorSet Npoints seed");
			System.exit(0);
		}
		double[][] data = genCantorSet26(Integer.parseInt(args[0]), Long.parseLong(args[1]));
		MyUtils.print_matrix(data);
	}

}
