package az.epsilon.rn;

import az.epsilon.RnFunction;

public interface Gradient {

    public abstract void setFunction(RnFunction f);

    public abstract void gradient(double[] gradients);

    public abstract double argmin(double alpha);

    public abstract double minimize();

    public abstract void calculate(double[] x0);

    public abstract double norm(double[] x);
    
    public void setEpsilon(double epsilon);
}
