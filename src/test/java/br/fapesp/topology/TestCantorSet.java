package br.fapesp.topology;

import static org.junit.Assert.*;

import org.junit.Ignore;
import org.junit.Test;

import br.fapesp.myutils.MyUtils;

public class TestCantorSet {

	@Ignore
	@Test
	public void testCantor() {
		double[][] data = CantorSet.genCantorSet26(10000, 1234);
		MyUtils.print_matrix(data);
	}
	
	@Ignore
	@Test
	public void testMST() {
		double[][] data = CantorSet.genCantorSet26(20, 1234);
		int[][] mst = MyUtils.computeMinimumSpanningTreePrim(data);
		MyUtils.print_matrix(mst);
	}
	
	
	@Test
	public void testMST2() {
		double[][] data = CantorSet.genCantorSet26(1000, 1234);
//		int[][] mst = MyUtils.computeMinimumSpanningTreePrim(data);
		int[][] mst = MyUtils.fastPrim(data);
//		MyUtils.print_matrix(mst);
	}

}
