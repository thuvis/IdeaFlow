package BCC_DTW;

public class GaussianIntegral {
	
	private double EPSILON = 0.01;
	private double mean ;
	private double variance;
	
	public GaussianIntegral(double mean,double variance){
		this.mean = mean;
		this.variance = variance;
	}
	
	public GaussianIntegral(){
		this.mean = 0.0;
		this.variance = 1.0;
	}
	
	public double f(double x){
		double y = Math.exp(- (x - mean)*(x - mean) / (2*variance)) / Math.sqrt(2*Math.PI*variance);
		if (Double.isNaN(y) || y < 1.0E-300)
			return (1.0E-300);
		return y;
	}
	
	public double f(double x, double mean, double variance) {
		this.mean = mean;
		this.variance = variance;
		return f(x);
	}
	
	public double integral(double a,double b){
		double N = 5;
		double step = (b - a) / N;
		double Q1 = 0.0;
		double Q2 = 0.0;
		if (b <= a)
			return 0.0;
        for (double i = a; i < b; i += step)
        	Q2 += f(i);
        Q2 *= step;
        do {
        	Q1 = Q2;
        	Q2 = 0.0;
        	step /= 2;
        	for (double i = a; i < b; i += step)
            	Q2 += f(i);
        	Q2 *= step;
        } while (Math.abs(Q2 - Q1) > EPSILON);
        return Q2;
	}
	
	public double integral(double mean, double variance, double a, double b){
		this.mean = mean;
		this.variance = variance;
		return integral(a, b);
	}
	
	public double integral_x(double a,double b){
		double N = 5;
		double step = (b - a) / N;
		double Q1 = 0.0;
		double Q2 = 0.0;
		if (b <= a)
			return 0.0;
        for (double i = a; i < b; i += step)
        	Q2 += i * f(i);
        Q2 *= step;
        do {
        	Q1 = Q2;
        	Q2 = 0.0;
        	step /= 2;
        	for (double i = a; i < b; i += step)
            	Q2 += i * f(i);
        	Q2 *= step;
        } while (Math.abs(Q2 - Q1) > EPSILON);
        return Q2;       
	}
	
	public double integral_x(double mean, double variance, double a, double b){
		this.mean = mean;
		this.variance = variance;
		return integral_x(a, b);
	}
	
	public double integral_x_2(double a,double b){
		double N = 5;
		double step = (b - a) / N;
		double Q1 = 0.0;
		double Q2 = 0.0;
		if (b <= a)
			return 0.0;
        for (double i = a; i < b; i += step)
        	Q2 += i * i * f(i);
        Q2 *= step;
        do {
        	Q1 = Q2;
        	Q2 = 0.0;
        	step /= 2;
        	for (double i = a; i < b; i += step)
            	Q2 += i * i * f(i);
        	Q2 *= step;
        } while (Math.abs(Q2 - Q1) > EPSILON);
        return Q2;
	}
	
	public double integral_x_2(double mean, double variance, double a, double b){
		this.mean = mean;
		this.variance = variance;
		return integral_x_2(a, b);
	}
	
	public static void main(String[] args) {
		GaussianIntegral a = new GaussianIntegral(21, 33*33);
		double result = a.integral(0, 1, -2,2);
		System.out.println(result);
	}
	
	
}
