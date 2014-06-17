package br.fapesp.topology;

public class FiltrationResult {

	/**
	 * number of N dimensional holes, starting from 0 (connected components), holes, voids, etc.
	 */
	public int[] nDholes;
	
	/**
	 * maximum lifetime of an N-dimensional hole, starting from 0.
	 */
	public double[] maxHoleLifeTime;
	
	/**
	 * average lifetime of a N-dimensional hole
	 */
	public double[] averageHoleLifeTime;
	
	/**
	 * get the number of revelant holes based on a user-defined threshold
	 */
	public int[] nDrelevantHoles;

}
