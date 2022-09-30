package gaussseidel;
/**
 * @author Mahsima Madanipour
 */

import java.io.*;
import java.util.Arrays;
import java.util.StringTokenizer;

public class GaussSeidel {
	public static final int MAX = 100;
	private double[][] MATRIX;
	
	/**	
	* MATRIX is a matrix with the following form:
	* (a11 * x1) + (a12 * x2) + ... + (a1n * xn) = b1
	* (a21 * x1) + (a22 * x2) + ... + (a2n * xn) = b2
	* .
	* .
	* .
	* (an1 * xn) + (an2 * x2) + ... + (ann * xn) = bn	
	*/
	
	public GaussSeidel(double [][] matrix) { 
		MATRIX = matrix; 
	}
	
	public static void main(String[] args) throws IOException{
		int num;
		double[][] mat;

		BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
		PrintWriter writer = new PrintWriter(System.out, true);

		System.out.println("Enter the number of equations: ");
		num = Integer.parseInt(reader.readLine());
		mat = new double[num][num+1];

		System.out.println("Enter fixed numbers(e.g. if you enter 3 then you should enter each line like 'a b c d' for a*x1 + b*x2 + c*x3 = d): ");
		for (int i = 0; i < num; i++) {
			StringTokenizer strtk = new StringTokenizer(reader.readLine());

			while (strtk.hasMoreTokens()){
				for (int j = 0; j < num + 1 && strtk.hasMoreTokens(); j++){
					mat[i][j] = Integer.parseInt(strtk.nextToken());
				}
			}
		}

    
		GaussSeidel gausSeidel = new GaussSeidel(mat);

		if (!gausSeidel.makeDominant()) {
			writer.println("The system isn't diagonally dominant(The method cannot guarantee convergence).");
		}

		writer.println();
		gausSeidel.solve();
	}
	
	public void solve(){
		int iterations = 0;
		int num = MATRIX.length;
		double epsilon = 1e-10;
		double[] X = new double[num]; // Approximations
		double[] P = new double[num]; // Prev
		Arrays.fill(X, 0);

		while (true) {
			for (int i = 0; i < num; i++) {
				double sum = MATRIX[i][num]; // b_n

				for (int j = 0; j < num; j++){
					if (j != i)
						sum -= MATRIX[i][j] * X[j];
				}

				// Update xi to use in the next row calculation
				X[i] = 1/MATRIX[i][i] * sum;   
			}

			System.out.print("X" + iterations + " = {");
			for (int i = 0; i < num; i++){
				if(i == 0){
					System.out.print(X[i]);
				}else{
					System.out.print(", " + X[i]);
				}
			}
			System.out.println("}");

			iterations++;
			if (iterations == 1) continue;

			boolean stop = true;
			for (int i = 0; i < num && stop; i++)
				if (Math.abs(X[i] - P[i]) > epsilon)
					stop = false;

			if (stop || iterations == MAX) break;
			P = (double[])X.clone();
		}
	}

	public boolean transformToDominant(int r, boolean[] V, int[] R){
		int num = MATRIX.length;
		if (r == MATRIX.length) {
			double[][] T = new double[num][num+1];
			for (int i = 0; i < R.length; i++) {
				for (int j = 0; j < num + 1; j++){
					T[i][j] = MATRIX[R[i]][j];
				}
			}

			MATRIX = T;
      
			return true;
		}

		for (int i = 0; i < num; i++) {
			if (V[i]) continue;

			double sum = 0;
      
			for (int j = 0; j < num; j++)
				sum += Math.abs(MATRIX[i][j]);

			if (2 * Math.abs(MATRIX[i][r]) > sum) { // diagonally dominant?
				V[i] = true;
				R[r] = i;

				if (transformToDominant(r + 1, V, R))
					return true;

				V[i] = false;
			}
		}

		return false;
	}
  
  
  /**
   * Returns true if is possible to transform M(data member) to a diagonally
   * dominant matrix, false otherwise.
  */
	public boolean makeDominant(){
		boolean[] visited = new boolean[MATRIX.length];
		int[] rows = new int[MATRIX.length];

		Arrays.fill(visited, false);

		return transformToDominant(0, visited, rows);
	}
}
