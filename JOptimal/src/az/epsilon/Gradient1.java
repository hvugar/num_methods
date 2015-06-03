package az.epsilon;

import java.util.Arrays;

public class Gradient1 implements Gradient {

    private double[] x;
    private double[] g;

    @Override
    public double f(double[] x) {
        double x1 = x[0];
        double x2 = x[1];
        return ((1 - x1) * (1 - x1)) + 100 * (x2 - x1 * x1) * (x2 - x1 * x1);
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

        for (int i = 0; i < x.length; i++) {
            x[i] = x[i] + alpha * g[i];
        }

        return alpha;
    }

    @Override
    public double minimize() {
        return argmin(0.0);
    }

    @Override
    public void calculate(double[] x0) {

        x = null;
        x = Arrays.copyOf(x0, x0.length);
        g = new double[x0.length];

        do {

            gradient(g);

            double g_norm = norm(g);
            if (g_norm < 0.00001) {
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
        double g_norm = 0.0;
        for (int i = 0; i < x.length; i++) {
            g_norm = g_norm + x[i]*x[i];
        }
        g_norm = Math.sqrt(g_norm);        
        return g_norm;
    }
}
