package az.epsilon;

import java.lang.reflect.Array;
import java.util.Arrays;

public class Gradient1 implements Gradient {

    private double[] x;

    @Override
    public double f(double[] x) {
        return x[0] * x[0] + x[1] * x[1];
    }

    @Override
    public void gradient(double[] gradients) {
        double h = 0.0001;
        
        for (int i = 0; i < gradients.length; i++) {
            
            double _x = x[i];
            
            x[i] = x[i] + h;
            double f1 = f(x);
            x[i] = x[i] - 2 * h;
            double f2 = f(x);
            x[i] = x[i] + h;
            
            x[i] = _x;

            gradients[i] = (f1 - f2) / (2.0 * h);
        }
    }

    @Override
    public double argmin(double alpha) {
        return alpha;
    }

    @Override
    public double minimize() {
        return argmin(0.0);
    }

    @Override
    public void calculate(double[] x0) {

        x = Arrays.copyOf(x0, x0.length);
        double[] g = new double[x0.length];


        do {

            gradient(g);

            if (norm(g) < 0.00001) {
                break;
            }

            double alpha = minimize();

            for (int i = 0; i < x.length; i++) {
                x[i] = x[i] - alpha * g[i];
            }

        } while (true);
    }

    @Override
    public double norm(double[] gradients) {
        return 0.0;
    }
}
