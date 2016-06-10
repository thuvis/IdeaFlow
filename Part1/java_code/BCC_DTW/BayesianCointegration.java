package BCC_DTW;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import BCC_DTW.GaussianIntegral;

//We need three assistant classes first:
//Mixture component in Cointegration
class rhoComponent {
	public int rho;  //weight
	public double integral;
	public double phi;
	public double phi_2;
	
	public rhoComponent() {
		this.integral = 0;
		this.integral = 0.0;
		this.phi = 0.0;
		this.phi_2 = 0.0;
	}
	
	public rhoComponent(int rho, double integral, double phi, double phi_2) {
		this.rho = rho;
		this.integral = integral;
		this.phi = phi;
		this.phi_2 = phi_2;
	}
	
	public rhoComponent(rhoComponent ori) {
		this.rho = ori.rho;
		this.integral = ori.integral;
		this.phi = ori.phi;
		this.phi_2 = ori.phi_2;
	}
	
	public void addIntegral(double integral) {
		this.integral += integral;
	}
	
	public void addPhi(double phi) {
		this.phi += phi;
	}
	
	public void addPhi_2(double phi_2) {
		this.phi_2 += phi_2;
	}
}

//Mixture component bucket in Cointegration
class rhoBucket {
	public double integral_sum;
	public double phi_sum;
	public double phi_2_sum;
	public ArrayList<rhoComponent> rhoComponentList;
	
	public rhoBucket() {
		this.integral_sum = 0.0;
		this.phi_sum = 0.0;
		this.phi_2_sum = 0.0;
		this.rhoComponentList = new ArrayList<rhoComponent>();
	}
	
	public rhoBucket(double integral_sum, double phi_sum, double phi_2_sum) {
		this.integral_sum = integral_sum;
		this.phi_sum = phi_sum;
		this.phi_2_sum = phi_2_sum;
		this.rhoComponentList = new ArrayList<rhoComponent>();
	}
}

//Result of DTW alignment
class ResultOfDTW {
	public double[] xs;  //New x series.
	public double[] ys;  //New y series.
	public int[] ori_index;  //Original index
	public int[] time;  //Propagation delay.
	public int length;
	
	public ResultOfDTW() {
		
	}
	
	public ResultOfDTW(int N) {
		xs = new double[N];
		ys = new double[N];
		ori_index = new int[N];
		time = new int[N];
		length = N;
	}
}

//Here is the main class.
//Bayesian Conditional Cointegration Class.
public class BayesianCointegration {
	private double tao[];  //Transition probability.
	private double a_estimate;  //Estimation for regression parameter a in Cointegration.
	private double b_estimate;  //Estimation for regression parameter b in Cointegration.
	private double sigma_2_estimate;  //Estimation for parameter sigma_2 in Cointegration.
	private double gamma_i_1[];  //Posterior probability series for Random Walk.
	private double gamma_i_0[];  //Posterior probability series for Cointegration.
	private double expectation_phi[];  //Expectation array of phi.
	private double  expectation_phi_2[];  //Expectation array of phi^2.
	private int iteration_maximum;  //Maximum times of EM iterations.
	private double iteration_epsilon;  //Stop EM if the relative difference between new and old likelihood is less than this value. 
	private double gamma_threshold;  //Posterior threshold. If and only if posterior probability less than this value will be determined to be local cointegration.
	private int DTW_window;  //Maximum allowable shift window of DTW alignment.
	private static int component_num;  //Number of mixture components that you want to save in every iteration.
	
	public BayesianCointegration(){
		tao = new double[4];
		tao[0] = 0.9;  //Transition probability from Cointegration to Cointegration.
		tao[1] = 0.1;  //Transition probability from Cointegration to Random Walk.
		tao[2] = 0.8;  //Transition probability from Random Walk to Cointegration.
		tao[3] = 0.2;  //Transition probability from Random Walk to Random Walk.
		iteration_maximum = 20;
		iteration_epsilon = 0.0005;
		gamma_threshold = 0.5;
		DTW_window = 6;
		component_num = 20;
	}
	
	//Get transition probability.
	public double[] getTao() {
		return tao;
	}
	
	//Set transition probability. Parameters: p(0->0) and p(1->1)
	public void setTao(double tao00, double tao11) {
		tao[0] = tao00;
		tao[1] = 1 - tao00;
		tao[2] = 1 - tao11;
		tao[3] = tao11;
	}
	
	//Get maximum times of EM iterations.
	public int getIterationMaximum() {
		return iteration_maximum;
	}
	
	//Set maximum times of EM iterations.
	public void setIterationMaximum(int v) {
		iteration_maximum = v;
	}
	
	//Get iteration epsilon. Stop EM if the relative difference between new and old likelihood is less than this value. 
	public double getIterationEpsilon() {
		return iteration_epsilon;
	}
	
	//Set iteration epsilon. Stop EM if the relative difference between new and old likelihood is less than this value. 
	public void setIterationEpsilon(double v) {
		iteration_epsilon = v;
	}
	
	//Get the number of mixture components that you want to save in every iteration.
	public int getComponent_num() {
		return component_num;
	}
	
	//Set the number of mixture components that you want to save in every iteration.
	public void setComponent_num(int num) {
		component_num = num;
	}
	
	//Get posterior threshold. If and only if posterior probability less than this value will be determined to be local cointegration.
	public double getGamma_threshold() {
		return gamma_threshold;
	}
	
	//Set posterior threshold. If and only if posterior probability less than this value will be determined to be local cointegration.
	public void setGamma_threshold(double v) {
		gamma_threshold = v;
	}
	
	//Get maximum allowable shift window of DTW alignment.
	public int getDTW_window() {
		return DTW_window;
	}
	
	//Set maximum allowable shift window of DTW alignment.
	public void setDTW_window(int v) {
		DTW_window = v;
	}
	
	//Get estimation for regression parameter a in Cointegration.
	public double getA() {
		return a_estimate;
	}
	
	//Get estimation for regression parameter b in Cointegration.
	public double getB() {
		return b_estimate;
	}
	
	//Get estimation for parameter sigma_2 in Cointegration.
	public double getSigma_2() {
		return sigma_2_estimate;
	}
	
	//Get posterior probability series for Random Walk.
	public double[] getGamma_i_1() {
		return gamma_i_1;
	}
	
	//Get posterior probability series for Cointegration.
	public double[] getGamma_i_0() {
		return gamma_i_0;
	}
	
	//Get posterior probability series for Random Walk.
	public double[] getExpectation_phi() {
		return expectation_phi;
	}
	
	//Get posterior probability series for Cointegration.
	public double[] getExpectation_phi_2() {
		return expectation_phi_2;
	}
	
	public static void printList(double list[]){
		for(int i = 0; i < list.length; i++){
			System.out.print(list[i] + " ");
		}
		System.out.println();
	}
	
	public static void printList_index(double list[]){
		for(int i = 0; i < list.length; i++){
			System.out.print("["+i+"]" + list[i] + " ");
		}
		System.out.println();
	}
	
	public static void printList(int list[]){
		for(int i = 0; i < list.length; i++){
			System.out.print(list[i] + " ");
		}
		System.out.println();
	}
	
	public static void printList_index(int list[]){
		for(int i = 0; i < list.length; i++){
			System.out.print("["+i+"]" + list[i] + " ");
		}
		System.out.println();
	}
	
	public static double calculateSum(double list[]){
		double sum = 0;
		for(int i = 0; i < list.length; i++){
			sum += list[i];
		}
		return sum;
	}
	
	public static double calculateSum(double list[], int start, int end){
		double sum = 0;
		for(int i = start - 1; i < end; i++){
			sum += list[i];
		}
		return sum;
	}
	
	public static double innerProduct(double list1[], double list2[]){
		double sum = 0;
		for(int i = 0; i < list1.length; i++){
			sum += list1[i] * list2[i];
		}
		return sum;
	}
	
	//Indices of start1, end1, start2 and end2 begin with 1.
	public static double innerProduct(double list1[], int start1, int end1, double list2[], int start2, int end2){
		double sum = 0;
		for(int i = 0; i < end1-start1+1; i++){
			sum += list1[i+start1-1] * list2[i+start2-1];
		}
		return sum;
	}
	
	public static double[] dotProduct(double list1[], double list2[]){
		double result[] = new double[list1.length];
		for(int i = 0; i < list1.length; i++){
			result[i] = list1[i] * list2[i];
		}
		return result;
	}
	
	//Indices of start1, end1, start2 and end2 begin with 1.
	public static double[] dotProduct(double list1[], int start1, int end1, double list2[], int start2, int end2){
		double result[] = new double[end1-start1+1];
		for(int i = 0; i < end1-start1+1; i++){
			result[i] = list1[i+start1-1] * list2[i+start2-1];
		}
		return result;
	}
	
	public static double[] calculateResidual(double y[], double x[], double a, double b){
		double result[] = new double[x.length];
		for(int i = 0; i < x.length; i++){
			result[i] = y[i] - (a+b*x[i]);
		}
		return result;
	}
	
	public static double calculateVariance(double residual[]){
		double ex = 0, ex_2 = 0;
		for(int i=0; i < residual.length; i++){
			ex += residual[i];
			ex_2 += residual[i] * residual[i];
		}
		ex /= residual.length;
		ex_2 /= residual.length;
		return (ex_2 - ex * ex);
	}
	
	public static double calculateExpectationOfPhi(double weight_rho_t[], double f[]){
		double sum = 0;
		for(int i = 0; i < weight_rho_t.length; i++){
			sum += weight_rho_t[i] * f[i];
		}
		return sum;
	}
	
	public static double calculateExpectationOfPhi2(double weight_rho_t[], double f[], double F[]){
		double sum = 0;
		for(int i = 0; i < weight_rho_t.length; i++){
			sum += weight_rho_t[i] * (f[i]*f[i] + F[i]);
		}
		return sum;
	}
	
	public static double calculateLikelihood(double Z[]){
		double sum = 0;
		for(int i = 0; i < Z.length; i++){
			sum += Math.log10(Z[i]);
		}
		return sum;
	}
	
	public static double calculateLikelihood_RW(double residual[], double sigma_2){
		double sum = 0;
		for(int i=1; i<residual.length; i++){
			sum += Math.log10((new GaussianIntegral().f(residual[i], residual[i-1], sigma_2)));
		}
		return sum;
	}
	
	public static double[] solveEquation(double[][] A, double[] B) {
		double[] result = new double[2];
		double[][] A_inv = new double[2][2];
		double a = A[0][0], b = A[0][1], c = A[1][0], d = A[1][1], p = B[0], q = B[1];
		double denominator = a*d - b*c;
		
		A_inv[0][0] = d / denominator;
		A_inv[0][1] = -b / denominator;
		A_inv[1][0] = -c / denominator;
		A_inv[1][1] = a / denominator;
		
		result[0] = A_inv[0][0] * p + A_inv[0][1] * q;
		result[1] = A_inv[1][0] * p + A_inv[1][1] * q;
		return result;
	}
	
	public static double[] calculateInitialParameters(double x[], double y[], int N){
		double[][] A = new double[2][2];
		double[] B = new double[2];
		double[] result;
		A[0][0] = N;
		double x_sum =  calculateSum(x);
		A[0][1] = x_sum;
		A[1][0] = x_sum;
		A[1][1] = innerProduct(x, x);
		B[0] = calculateSum(y);
		B[1] = innerProduct(x, y);
		result = solveEquation(A, B);
		return result;
	}
	
	public static double eStep(int N, double sigma_2, double tao[], double residual[], double alpha_i_0[], double alpha_i_1[], double gamma_i_0[], double gamma_i_1[], double Expectation_phi[], double Expectation_phi_2[], double Z[]){
		double likelihood = 0;
		int indexN = component_num;  //Mixture component number that you want to save in every iteration.
		double f[][] = new double[N][indexN];
		double F[][] = new double[N][indexN];
		double weight_rho_t[][] = new double[N][indexN];  //Weights of components.
		int rho_t[][] = new int[N][indexN];  //Cumulative indices of components.
		double normalize_t[][] = new double[N][indexN];  //Integral Factors for normalization.
		likelihood = filtering(N, sigma_2, tao, residual, alpha_i_0, alpha_i_1, Expectation_phi, Expectation_phi_2, Z, indexN, f, F, weight_rho_t, rho_t, normalize_t);
		smoothing(N, tao, alpha_i_0, alpha_i_1, gamma_i_0, gamma_i_1, Expectation_phi, Expectation_phi_2, indexN, f, F, weight_rho_t, rho_t, normalize_t);
		return likelihood;
	}

	public static double filtering(int N, double sigma_2, double tao[] , double residual[], double alpha_i_0[], double alpha_i_1[], double Expectation_phi[], double Expectation_phi_2[], double Z[], int indexN, double f[][], double F[][], double weight_rho_t[][], int rho_t[][], double normalize_t[][]){
		
		GaussianIntegral GI = new GaussianIntegral();
		Z[0] = 1;
		
		if (residual[0] == 0.0)
			residual[0] = 1.00E-20;
		f[1][0] = residual[1] / residual[0];
		F[1][0] = sigma_2 / (residual[0]*residual[0]);
		
		alpha_i_0[0] = 0.5;  //Probability of Cointegration at the beginning
		alpha_i_1[0] = 0.5;  //Probability of Random Walk at the beginning
		alpha_i_1[1] = GI.f(residual[1], residual[0], sigma_2) * (alpha_i_0[0]*tao[1] + alpha_i_1[0]*tao[3]);
		
		weight_rho_t[1][0] = (tao[2]*alpha_i_1[0]*0.5 + tao[0]*alpha_i_0[0]*0.5) / Math.abs(residual[0]);  //0.5 for the integral factor of uniform distribution of emission term
		normalize_t[1][0] = GI.integral(f[1][0], F[1][0], -1, 1);  //Integral factor: utilize numerical integral to calculate.
		weight_rho_t[1][0] = weight_rho_t[1][0] * normalize_t[1][0];
		alpha_i_0[1] = weight_rho_t[1][0];
		
		Z[1] = alpha_i_0[1] + alpha_i_1[1];
		alpha_i_0[1] = alpha_i_0[1] / Z[1];
		alpha_i_1[1] = alpha_i_1[1] / Z[1];
		weight_rho_t[1][0] = weight_rho_t[1][0] / (alpha_i_0[1]*Z[1]);
		rho_t[1][0] = 0;
		
		for(int n = 2; n < N; n++){
			int index = indexN;
			if(n <= indexN)
				index = n-1;
			double f_temp[]= new double[index+1];
			double F_temp[]= new double[index+1];
			double weight_temp[]= new double[index+1];
			int rho_temp[] = new int[index+1];
			double normalize_temp[]= new double[index+1];
			for(int m = 0; m < index; m++){
				double denominator = sigma_2 + residual[n-1]*residual[n-1]*F[n-1][m];
				f_temp[m] = (f[n-1][m]*sigma_2 + residual[n]*residual[n-1]*F[n-1][m]) / denominator;
				F_temp[m] = sigma_2*F[n-1][m] / denominator;
				weight_temp[m] =  weight_rho_t[n-1][m]*tao[0]*alpha_i_0[n-1]*GI.f(residual[n], residual[n-1]*f[n-1][m], denominator);
				rho_temp[m] = rho_t[n-1][m] + 1;
			}
			if (residual[n-1] == 0.0)
				residual[n-1] = 1.00E-20;
			f_temp[index] = residual[n] / residual[n-1];
			F_temp[index] = sigma_2 / (residual[n-1]*residual[n-1]);
			weight_temp[index] = alpha_i_1[n-1]*tao[2] * 0.5 / Math.abs(residual[n-1]);
			rho_temp[index] = 0;

			for (int m = 0; m < index + 1; m++)
			{
				normalize_temp[m]= GI.integral(f_temp[m], F_temp[m], -1, 1);
				weight_temp[m]= weight_temp[m] * normalize_temp[m];
			}
			
			double weight_sum = 0.0;
			if(index < indexN) {
				for(int m = 0; m < index+1; m++){
					f[n][m] = f_temp[m];
					F[n][m] = F_temp[m];
					weight_rho_t[n][m] = weight_temp[m];
					rho_t[n][m] = rho_temp[m];
					normalize_t[n][m] = normalize_temp[m];
					weight_sum += weight_temp[m];
				}
				alpha_i_0[n]= weight_sum;
			}
			else {
				int flag=0;
				double min = weight_temp[0];
				for(int m = 1; m < index + 1; m++){
					if(weight_temp[m]<min){
						min = weight_temp[m];
						flag = m;
					}
				}
				for(int m=0; m < index+1; m++){
					if(m == flag)
						continue;
					int temp=0;
					if (m < flag)
						temp = m;
					else if(m > flag)
						temp = m - 1;
					f[n][temp] = f_temp[m];
					F[n][temp] = F_temp[m];
					weight_rho_t[n][temp] = weight_temp[m];
					rho_t[n][temp] = rho_temp[m];
					normalize_t[n][temp] = normalize_temp[m];
					weight_sum += weight_temp[m];
				}
				alpha_i_0[n]= weight_sum;
			}
			
			alpha_i_1[n]= GI.f(residual[n], residual[n-1], sigma_2) * (alpha_i_1[n-1]*tao[3] + alpha_i_0[n-1]*tao[1]);
			
			Z[n] = alpha_i_0[n] + alpha_i_1[n];
		    alpha_i_0[n] = alpha_i_0[n] / Z[n];
		    alpha_i_1[n] = alpha_i_1[n] / Z[n];
		    if (index < indexN)
		    	index = index+1;
		    for (int m = 0; m < index; m++)
		    	weight_rho_t[n][m] = weight_rho_t[n][m] / (Z[n]*alpha_i_0[n]);
		}
		double likelihood = calculateLikelihood(Z);
		return likelihood;
		
	}
	
	public static void smoothing(int N, double tao[], double alpha_i_0[], double alpha_i_1[], double gamma_i_0[], double gamma_i_1[], double Expectation_phi[], double Expectation_phi_2[], int indexN, double f[][], double F[][], double weight_rho_t[][], int rho_t[][], double normalize_t[][]){
		
		GaussianIntegral GI = new GaussianIntegral();
		rhoBucket[] gammaBucket = new rhoBucket[N];
		int length;
		int flag;
		rhoComponent temp = new rhoComponent();
		
		for (int m = 0; m < N; m++)
			gammaBucket[m] = new rhoBucket();
		
		gamma_i_0[N-1] = alpha_i_0[N-1];
		gamma_i_1[N-1] = alpha_i_1[N-1];
		length = N > indexN ? indexN : (N-1);
		for (int m = 0; m < length; m++)
		{
			double integral_temp;
			double phi_temp;
			double phi_2_temp;
			integral_temp = weight_rho_t[N-1][m] * alpha_i_0[N-1];
			phi_temp = GI.integral_x(f[N-1][m], F[N-1][m], -1, 1) * weight_rho_t[N-1][m] * alpha_i_0[N-1] / normalize_t[N-1][m];
			phi_2_temp = GI.integral_x_2(f[N-1][m], F[N-1][m], -1, 1) * weight_rho_t[N-1][m] * alpha_i_0[N-1] / normalize_t[N-1][m];
			gammaBucket[N-1].rhoComponentList.add(new rhoComponent(rho_t[N-1][m], integral_temp, phi_temp, phi_2_temp));
			gammaBucket[N-1].integral_sum += integral_temp;
			gammaBucket[N-1].phi_sum += phi_temp;
			gammaBucket[N-1].phi_2_sum += phi_2_temp;
		}
		Expectation_phi[N-1] = gamma_i_1[N-1] + gammaBucket[N-1].phi_sum;
		Expectation_phi_2[N-1] = gamma_i_1[N-1] + gammaBucket[N-1].phi_2_sum;
	
		for(int n = N-2; n >= 0; n--){
			double integral_temp;
			double phi_temp;
			double phi_2_temp;
			double prefactor;
			length = gammaBucket[n+1].rhoComponentList.size();
			flag = -1;
			for (int m = 0; m < length; m++) {
				temp = gammaBucket[n+1].rhoComponentList.get(m);
				if (temp.rho == 0)
					flag = m;
				else {
					gammaBucket[n].rhoComponentList.add(new rhoComponent(temp.rho-1, temp.integral, temp.phi, temp.phi_2));
					gammaBucket[n].integral_sum += temp.integral;
					gammaBucket[n].phi_sum += temp.phi;
					gammaBucket[n].phi_2_sum += temp.phi_2;
				}
			}
			if (flag == -1) {
				gamma_i_1[n] = gamma_i_1[n+1]*tao[3]*alpha_i_1[n] / (tao[3]*alpha_i_1[n] + tao[1]*alpha_i_0[n]);
			}
			else {
				temp = gammaBucket[n+1].rhoComponentList.get(flag);
				gamma_i_1[n] = temp.integral + gamma_i_1[n+1]*tao[3]*alpha_i_1[n] / (tao[3]*alpha_i_1[n] + tao[1]*alpha_i_0[n]);
			}
			
			if (n > indexN)
				length = indexN;
			else
				length = n;
			prefactor = gamma_i_1[n+1]*alpha_i_0[n]*tao[1] / (alpha_i_0[n]*tao[1] + alpha_i_1[n]*tao[3]);
			for (int m = 0; m < length; m++) {
				integral_temp = weight_rho_t[n][m] * prefactor;
				phi_temp = GI.integral_x(f[n][m], F[n][m], -1, 1) * weight_rho_t[n][m] * prefactor / normalize_t[n][m];
				phi_2_temp = GI.integral_x_2(f[n][m], F[n][m], -1, 1) * weight_rho_t[n][m] * prefactor / normalize_t[n][m];
				int len = gammaBucket[n].rhoComponentList.size();
				flag = -1;
				for (int k = 0; k < len; k++) {
					temp = gammaBucket[n].rhoComponentList.get(k);
					if (rho_t[n][m] == temp.rho)
					{
						flag = k;
						temp.addIntegral(integral_temp);
						temp.addPhi(phi_temp);
						temp.addPhi_2(phi_2_temp);
						gammaBucket[n].integral_sum += integral_temp;
						gammaBucket[n].phi_sum += phi_temp;
						gammaBucket[n].phi_2_sum += phi_2_temp;
						break;
					}
				}
				if (flag == -1) {
					gammaBucket[n].rhoComponentList.add(new rhoComponent(rho_t[n][m], integral_temp, phi_temp, phi_2_temp));
					gammaBucket[n].integral_sum += integral_temp;
					gammaBucket[n].phi_sum += phi_temp;
					gammaBucket[n].phi_2_sum += phi_2_temp;
				}
			}
			
			if (n == 0) {
				prefactor = gamma_i_1[n+1]*alpha_i_0[n]*tao[1] / (alpha_i_0[n]*tao[1] + alpha_i_1[n]*tao[3]);
				gamma_i_0[n] =  prefactor;
				Expectation_phi[n] = gamma_i_1[n];
				Expectation_phi_2[n] = gamma_i_1[n] + prefactor * 1.0/3.0;
			}
			else {
				gamma_i_0[n] = gammaBucket[n].integral_sum;
				Expectation_phi[n] = gamma_i_1[n] + gammaBucket[n].phi_sum;
				Expectation_phi_2[n] = gamma_i_1[n] + gammaBucket[n].phi_2_sum;
			}
		}
	}
	
	
	public static double[] mStep(int N ,double sigma_2, double x[],double y[],double tao[],double residual[],double alpha_i_0[],double alpha_i_1[],double gamma_i_0[],double gamma_i_1[],double Expectation_phi[],double Expectation_phi_2[],double Z[]){
		
		double a11 = (N-1) + calculateSum(Expectation_phi_2,2,N)- 2*calculateSum(Expectation_phi,2,N) + 1-Expectation_phi_2[0];
		
		double a12 = calculateSum(x,2,N)-innerProduct(x,1,N-1,Expectation_phi,2,N)-innerProduct(x,2,N,Expectation_phi,2,N) + 
				innerProduct(x,1,N-1,Expectation_phi_2,2,N) +  x[0]*(1-Expectation_phi_2[0]);
		
		double a21 = a12;
		
		double a22 = innerProduct(x,2,N,x,2,N)-2*innerProduct(dotProduct(x,1,N-1,x,2,N),1,N-1,Expectation_phi,2,N) + innerProduct(dotProduct(x,1,N-1,x,1,N-1),1,N-1,Expectation_phi_2,2,N) + x[0]*x[0] *(1-Expectation_phi_2[0]);  ;
		
		double b1 = calculateSum(y,2,N)-innerProduct(y,1,N-1,Expectation_phi,2,N) - innerProduct(y,2,N,Expectation_phi,2,N)
				+ innerProduct(y,1,N-1,Expectation_phi_2,2,N) + y[0]*(1-Expectation_phi_2[0]) ;
		
		double b2 = innerProduct(x,2,N,y,2,N) - innerProduct(dotProduct(x,1,N-1,y,2,N),1,N-1,Expectation_phi,2,N) - 
				innerProduct(dotProduct(x,2,N,y,1,N-1),1,N-1,Expectation_phi,2,N) + 
				innerProduct(dotProduct(x,1,N-1,y,1,N-1),1,N-1,Expectation_phi_2,2,N)+
				x[0]*y[0]*(1-Expectation_phi_2[0]);
		
		double[][] A = new double[2][2];
		double[] B = new double[2];
		double[] result;
		A[0][0] = a11;
		A[0][1] = a12;
		A[1][0] = a21;
		A[1][1] = a22;
		B[0] = b1;
		B[1] = b2;
		result = solveEquation(A, B);
		double a_new,b_new,sigma_2_new;
		a_new = result[0];
		b_new = result[1];
		residual = calculateResidual(y, x, a_new, b_new);
		sigma_2_new = residual[0]*residual[0]*(1-Expectation_phi_2[0]) + innerProduct(residual,2,N,residual,2,N) - 
				2*innerProduct(dotProduct(residual,1,N-1,residual,2,N),1,N-1,Expectation_phi,2,N) +
				innerProduct(dotProduct(residual,1,N-1,residual,1,N-1),1,N-1,Expectation_phi_2,2,N);
		sigma_2_new = sigma_2_new / N;
		
		double parameters[] = new double[3];
		parameters[0] = a_new;
		parameters[1] = b_new;
		parameters[2] = sigma_2_new;
		
		return parameters;
	}

	//The function to calculate whether series x and y cointegrated.
	//The output is an array of 3 values: log10(likelihood_Cointegration); log10(likelihood_RandomWalk); log10(likelihood_RandomWalk / likelihood_Cointegration).
	public double[] calculateBayesianCointegration(double x[],double y[]){
		
		if(x.length != y.length){
			System.out.println("ERROR! : The input time series have different length!");
			return null;
		}
		
		int N = x.length;
		double[] initialParameters = calculateInitialParameters(x, y, N);
		double a = initialParameters[0];
		double b = initialParameters[1];
		double residual[] = calculateResidual(y, x, a, b);
		double sigma_2 = calculateVariance(residual);
		double res[] = new double[3];
		
		if (sigma_2 <= 1.0E-10) {
			a_estimate = a;
			b_estimate = b;
			sigma_2_estimate = 0.0;
			this.gamma_i_1 = new double[N];
			this.gamma_i_0 = new double[N];
			for (int i = 0; i < N; i++) {
				this.gamma_i_0[i] = 1.0;
				this.gamma_i_1[i] = 0.0;
			}
			res[0] = Double.MAX_VALUE;
			res[1] = Double.MIN_VALUE;
			res[2] = Double.MIN_VALUE;
			return res;
		}
		
		int counter_Maximum = iteration_maximum;
		int counter = 0;
		double likelihood = 0;
		double likelihood_old = 0;
		double likelihood_record[] = new double[counter_Maximum];
		double tao[] = new double[4];
		tao[0] = this.tao[0];//tao00
		tao[1] = this.tao[1]; //tao10
		tao[2] = this.tao[2]; //tao01
		tao[3] = this.tao[3];//tao11
		
		double alpha_i_0[] = new double[N];
		double alpha_i_1[] = new double[N];
		double gamma_i_0[] = new double[N];
		double gamma_i_1[] = new double[N];
		
		double Expectation_phi[] = new double[N];
		double Expectation_phi_2[] = new double[N];
		double Z[] = new double[N];
		
		double a_new = 0;
		double b_new = 0;
		double sigma_2_new = 0;
		
		double sigma_2_RW = 0;
		double likelihood_RW = 0;
		double possibility_RW = 0;
		
		while(counter == 0 || Math.abs(( likelihood - likelihood_old ) / likelihood_old) > iteration_epsilon){
			counter++;
			if(counter > counter_Maximum )
				break;
			if(counter > 1){
				a = a_new;
				b = b_new;
				sigma_2 = sigma_2_new;
			}
			
			likelihood_old = likelihood;
			residual = calculateResidual(y, x, a, b);
			
			likelihood = eStep(N, sigma_2, tao, residual, alpha_i_0, alpha_i_1, gamma_i_0, gamma_i_1, Expectation_phi, Expectation_phi_2, Z);
			likelihood_record[counter-1] = likelihood;
			double parameters[] = new double[3]; 
			parameters = mStep(N, sigma_2, x, y, tao, residual, alpha_i_0, alpha_i_1, gamma_i_0, gamma_i_1, Expectation_phi, Expectation_phi_2, Z);
			
			a_new = parameters[0];
			b_new = parameters[1];
			if (parameters[2] > 0)
				sigma_2_new = parameters[2];
			else {
				sigma_2_new = 0.0;
				break;
			}
		}
		
		a_estimate = a_new;
		b_estimate = b_new;
		sigma_2_estimate = sigma_2_new;
		sigma_2_RW = innerProduct(residual,2,N,residual,2,N) + innerProduct(residual,1,N-1,residual,1,N-1)
				- 2 * innerProduct(residual,2,N,residual,1,N-1);
		sigma_2_RW = sigma_2_RW / (N-1);
		likelihood_RW = calculateLikelihood_RW(residual, sigma_2_RW);
		possibility_RW = likelihood_RW - likelihood;
		
		res[0] = likelihood;
		res[1] = likelihood_RW;
		res[2] = possibility_RW;
		this.gamma_i_1 = gamma_i_1;
		this.gamma_i_0 = gamma_i_0;
		this.expectation_phi = Expectation_phi;
		this.expectation_phi_2 = Expectation_phi_2;
		return (res);
	}
	
	public double[] subset_se(double[] ori, int start, int end) {
		double res[] = new double[end-start+1];
		if (start < 0)
			start = 0;
		if (end > ori.length-1)
			end = ori.length-1;
		for (int i = start, j = 0; i <= end; i++, j++)
			res[j] = ori[i];
		return res;
	}
	
	//Align series x and y with DTW.
	//If one time point aligns to several, take the mean of the several time points as the new value in new series.
	//normalizeFlag == 1 means standard normalizing the series before alignment.
	//trendFlag == 1 means additionally using the first and second derivative to calculate cost in DTW. (Only works if normalizeFlag == 1.)
	public ResultOfDTW DTW_MeanNormalized(double x[], double y[], int normalizeFlag, int trendFlag) {
		int n = x.length, m = y.length;
		double[][] DTW = new double[n+1][m+1];
		int[][] path = new int[n+1][m+1];
		int w = Math.max(DTW_window, Math.abs(n-m));
		int i, j;
		
		for (i = 0; i <= n; i++) {
			for (j = 0; j <= m; j++)
				DTW[i][j] = Double.MAX_VALUE / 2.0;
		}
		DTW[0][0] = 0;
		
		double[] xn = x.clone();
		double[] yn = y.clone();
		double[] xd = x.clone();
		double[] yd = y.clone();
		double[] xdd = x.clone();
		double[] ydd = y.clone();
		
		if (normalizeFlag == 1) {
			standard_normalize(x, xn);
			standard_normalize(y, yn);
			if (trendFlag == 1) {
				for (i = 0; i < n - 1; i++) {
					xd[i] = xn[i+1] - xn[i];
					yd[i] = yn[i+1] - yn[i];
				}
				if (n-1 >= 0)
					xd[n-1] = yd[n-1] = 0;
				for (i = 0; i < n - 2; i++) {
					xdd[i] = xd[i+1] - xd[i];
					ydd[i] = yd[i+1] - yd[i];
				}
				if (n-2 >= 0)
					xdd[n-2] = ydd[n-2] = 0;
				if (n-1 >= 0)
					xdd[n-1] = ydd[n-1] = 0;
			}
		}
		
		double weight = 1;
		for (i = 1; i <= n; i++) {
			int st = Math.max(1, i-w), ed = Math.min(m, i+w);
			for (j = st; j <= ed; j++) {
				double cost = Math.abs(xn[i-1] - yn[j-1]);
				if (normalizeFlag == 1 && trendFlag == 1) {
					double h1 = xd[i-1], h2 = yd[j-1], h3 = xdd[i-1], h4 = ydd[j-1];
					double cos_sim1 = (1 + h1*h2) / Math.sqrt((1 + h1*h1) * (1 + h2*h2));
					double cos_sim2 = (1 + h3*h4) / Math.sqrt((1 + h3*h3) * (1 + h4*h4));
					cost += 2 - (cos_sim1 + cos_sim2);
				}
				DTW[i][j] = Math.min(DTW[i-1][j]+cost, Math.min(DTW[i][j-1]+cost, DTW[i-1][j-1]+weight*cost));
				if (DTW[i][j] == DTW[i-1][j-1]+weight*cost)
					path[i][j] = -1;
				else if (DTW[i][j] == DTW[i][j-1]+cost)
					path[i][j] = 1;
				else if (DTW[i][j] == DTW[i-1][j]+cost)
					path[i][j] = 0;
			}
		}
		
		//Calculate lead lag time array.
		ArrayList<Double> xs = new ArrayList<Double>(), ys = new ArrayList<Double>();
		ArrayList<Integer> ori_index = new ArrayList<Integer>(), time = new ArrayList<Integer>();
		ArrayList<Integer> route = new ArrayList<Integer>();
		i = n; j = m;
		while (i != 1 || j != 1) {
			if (path[i][j] == -1) {
				route.add(-1);
				i--; j--;
			}
			else if (path[i][j] == 0) {
				route.add(0);
				i--;
			}
			else if (path[i][j] == 1) {
				route.add(1);
				j--;
			}
		}
		
		double flag = -1;
		int shiftType;
		int st1, ed1, st2, ed2;
		i = 0; j = 0;
		st1 = ed1 = st2 = ed2 = 0;
		ori_index.add(Math.min(i, j));
		time.add(j-i);
		for (int k = route.size()-1; k >= 0; k--) {
			shiftType = route.get(k);
			if (shiftType == -1) {
				i++; j++;
				xs.add(mean(x, st1, ed1));
				ys.add(mean(y, st2, ed2));
				ori_index.add(Math.min(i, j));
				time.add(j-i);
				flag = -1;
				st1 = ed1 = i;
				st2 = ed2 = j;
			}
			else if (shiftType == 0) {
				i++;
				if (flag == -1) {
					flag = 0;
					ed1++;
				}
				else if (flag == 1) {
					xs.add(mean(x, st1, ed1));
					ys.add(mean(y, st2, ed2));
					ori_index.add(Math.min(i, j));
					time.add(j-i);
					flag = -1;
					st1 = ed1 = i;
					st2 = ed2 = j;
				}
				else {
					ed1++;
				}
			}
			else if (shiftType == 1) {
				j++;
				if (flag == -1) {
					flag = 1;
					ed2++;
				}
				else if (flag == 0) {
					xs.add(mean(x, st1, ed1));
					ys.add(mean(y, st2, ed2));
					ori_index.add(Math.min(i, j));
					time.add(j-i);
					flag = -1;
					st1 = ed1 = i;
					st2 = ed2 = j;
				}
				else {
					ed2++;
				}
			}
		}
		xs.add(mean(x, st1, ed1));
		ys.add(mean(y, st2, ed2));
		
		int N = xs.size();
		ResultOfDTW res = new ResultOfDTW(N);
		for (i = 0; i < N; i++) {
			res.xs[i] = xs.get(i);
			res.ys[i] = ys.get(i);
			res.ori_index[i] = ori_index.get(i);
			res.time[i] = time.get(i);
		}
		
		return res;
	}
	
	public void standard_normalize(double[] x, double[] xs) {
		double sum = 0;
		double mean, std_v;
		int N = x.length;
		
		for (int i = 0; i < N; i++)
			sum += x[i];
		mean = sum / N;
		sum = 0;
		for (int i = 0; i < N; i++)
			sum += (x[i] - mean)*(x[i] - mean);
		std_v = Math.sqrt(sum / N);
		for (int i = 0; i < N; i++)
			xs[i] = (x[i] - mean) / std_v;
	}
	
	public double mean(double[] x, int st, int ed) {
		if (st > ed)
			return 0;
		
		double sum = 0;
		double mean;
		for (int i = st; i <= ed; i++)
			sum += x[i];
		mean = sum / (ed - st + 1);
		return mean;
	}

	//Augmented BCC, combined with DTW.
	//This function will implicitly use parameters DTW_window and gamma_threshold. (Please set them by member functions before calculation.)
	//Return the lead lag time array. (NaN means not Cointegrated.)
	public double[] cointegration_DTW_MeanNormalized(double x[], double y[]) {
		
		if(x.length != y.length){
			System.out.println("ERROR! : The input time series have different length!");
			return null;
		}
		
		int N = x.length;
		double res[] = new double[N];
		double gamma_i_1_1[], gamma_i_1_2[];
		ResultOfDTW result;
		int i, j;
		for (i = 0; i < N; i++)
			res[i] = Double.NaN;
		
		double xs[], ys[];
		
		xs = x.clone();
		ys = y.clone();
		result = DTW_MeanNormalized(xs, ys, 1, 1);
		
		calculateBayesianCointegration(result.ys, result.xs);
		gamma_i_1_1 = getGamma_i_1();
		calculateBayesianCointegration(result.xs, result.ys);
		gamma_i_1_2 = getGamma_i_1();
		
		//Smooth.
		for (i = 0; i < result.length; i++) {
			double gamma = (gamma_i_1_1[i] + gamma_i_1_2[i]) / 2;
			if (gamma < gamma_threshold) {
				res[result.ori_index[i]] = result.time[i];
				if (i < result.length-1 && result.ori_index[i] < result.ori_index[i+1]-1) {
					double d = (result.time[i+1] - result.time[i]) / ((double)(result.ori_index[i+1] - result.ori_index[i]));
					for (j = result.ori_index[i]+1; j < result.ori_index[i+1]; j++) {
						res[j] = (int)Math.round(result.time[i] + d * (j - result.ori_index[i]));
					}
				}
			}
		}
		return res;
	}
	
	//Augmented BCC, combined with DTW.
	//Parameters: x series, y series, DTW_window, gamma_threshold. (It won't change the member variables.)
	//Return the lead lag time array. (NaN means not Cointegrated.)
	public double[] cointegration_DTW_MeanNormalized(double x[], double y[], int DTW_win, double gamma_thr) {
		int window_member = DTW_window;
		double gamma_member = gamma_threshold;
		
		DTW_window = DTW_win;
		gamma_threshold = gamma_thr;
		double[] res = cointegration_DTW_MeanNormalized(x, y);
		
		DTW_window = window_member;
		gamma_threshold = gamma_member;
		return res;
	}
	
	//Output an index matrix according to the leadLagTime tensor between words of two social group.
	//This output file can be used as the input file for the tensor-based clustering algorithm (MATLAB code).
	public void outputIndexMatrix(String fileName, double[][][] leadLagTime, int taoMax) throws IOException {
		int wordNum1 = leadLagTime.length, wordNum2 = leadLagTime[0].length, timepointNum = leadLagTime[0][0].length;
		File file = new File(fileName);
		if (!file.exists()) {
			file.createNewFile();
		}
		FileWriter fileWritter = new FileWriter(file.getName(), false);
		BufferedWriter bufferWritter = new BufferedWriter(fileWritter);
		
		bufferWritter.write(wordNum1 + " " + wordNum2 + " " + timepointNum + " " + (2 * taoMax + 1) + "\r\n");
		for (int i = 0; i < wordNum1; i++) {
			double[][] subTensor1 = leadLagTime[i];
			for (int j = 0; j < wordNum2; j++) {
				double[] subTensor2 = subTensor1[j];
				for (int k = 0; k < timepointNum; k++) {
					if (!Double.isNaN(subTensor2[k]))
						bufferWritter.write((i + 1) + " " + (j + 1) + " " + (k + 1) + " " + ((int)subTensor2[k] + taoMax + 1) + "\r\n");
				}
			}
		}
		
    	bufferWritter.close();
    	fileWritter.close();
	}
	
	//This function gives examples about how to use the functions above.
	public static void main(String[] args) throws IOException{
		
		double x[] = {71,73,75,80,80,80,77,75,74,72,71,71,71,72,74,75,75,68,75,75,74,72,71,70,70,69,68,68,72,73,78,79,80,80,78};
		double y[] = x.clone();
		for (int i = 0; i < y.length; i++) {
			y[i] = x[i] * 2 + 25 + (Math.random()-0.5)*10;
		}
		
		BayesianCointegration bcc = new BayesianCointegration();
		bcc.setComponent_num(15);  //Set the number of mixture components that you want to save in every iteration.
		bcc.setIterationMaximum(15);  //Set maximum times of EM iterations.
		bcc.setIterationEpsilon(0.0005);  //Set iteration epsilon. Stop EM if the relative difference between new and old likelihood is less than this value. 
		bcc.setTao(0.9, 0.5);  //Set transition probability: p(Cointegration->Cointegration), p(Random Walk->Random Walk);
		bcc.calculateBayesianCointegration(x, y);  //Calculate BCC of x and y series.
		
		System.out.print("x: ");
		printList_index(x);
		System.out.print("y: ");
		printList_index(y);
		System.out.println("a = " + bcc.getA() + "; b = " + bcc.getB() + "; sigma_2 = " + bcc.getSigma_2());  //Print estimation of regression parameters.
		System.out.print("residual: ");
		printList_index(calculateResidual(y, x, bcc.getA(), bcc.getB()));  //Print final residual.
		System.out.print("expectation_phi: ");
		printList_index(bcc.getExpectation_phi());  //Print expectation of phi.
		System.out.print("gamma_i_1: ");
		printList_index(bcc.getGamma_i_1());  //Print posterior probability series for Random Walk.
		System.out.print("gamma_i_0: ");
		printList_index(bcc.getGamma_i_0());  //Print posterior probability series for Cointegration.
		
		for (int i = 0; i < y.length; i++) {
			y[i] = x[i] * 1.2 + 13;
		}
		for (int i = 2; i < 16; i++)
			y[i] = y[i+2];
		for (int i = 33; i >24; i--)
			y[i] = y[i-3];
		
		bcc.setDTW_window(6);  //Set maximum allowable shift window of DTW alignment.
		bcc.setGamma_threshold(0.5);  //Set posterior threshold. If and only if posterior probability less than this value will be determined to be local cointegration.
		double[] leadLagTime = bcc.cointegration_DTW_MeanNormalized(x, y);  //Advanced BCC with DTW. Return the lead-lag time array. (NaN for not cointegrated.)
		//bcc.cointegration_DTW_MeanNormalized(x, y, 5, 0.45);  //It will execute with these parameters but won't change the member variables DTW_window and gamma_threshold.
		
		System.out.println();
		System.out.print("x: ");
		printList_index(x);
		System.out.print("y: ");
		printList_index(y);
		System.out.println("a = " + bcc.getA() + "; b = " + bcc.getB() + "; sigma_2 = " + bcc.getSigma_2());  //Print estimation of regression parameters.
		System.out.print("lead-lag time: ");
		printList_index(leadLagTime);  //"leadLagTime[23] = 3.0" means x[23] lead y[23+3] = y[26]; "leadLagTime[2] = -1.0" means y[2] lead x[2+1] = x[3]. (Sign indicate which is leading and absolute value indicate lead-lag time. NaN for not cointegrated.)
		
		
		//To illustrate the output function for the index matrix, we first create a virtual lead-lag time tensor.
		//Imagine the first social group has 4 words, the second group has 5 words.
		int wordNum1 = 4, wordNum2 = 5, timepointNum = leadLagTime.length;
		double lltTensor[][][] = new double[wordNum1][wordNum2][timepointNum];
		
		//Imagine the first 2 words from the first group form idea1 and the first 3 words from the second form idea2.
		//Idea1 and idea2 have flow relationship, so their words have similar lead-lag times; other word pairs between these two groups don't.
		for (int i = 0; i < wordNum1; i++)
			for (int j = 0; j < wordNum2; j++)
				for (int k = 0; k < timepointNum; k++)
					lltTensor[i][j][k] = Double.NaN;
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 3; j++)
				for (int k = 0; k < timepointNum; k++)
					lltTensor[i][j][k] = leadLagTime[k];
		
		//Output the index matrix file. It can be used as the input file for the other part of our implementation (MATLAB code for tensor-based clustering).
		bcc.outputIndexMatrix("indexMatrix.txt", lltTensor, 6);
	}
	
}
