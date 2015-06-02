package az.epsilon;

public interface Gradient {

    double f(double[] x);

    void gradient(double[] gradients);
    
    double argmin(double alpha);

    double minimize();
    
    void calculate(double[] x0);
    
    double norm(double[] gradients);
}
