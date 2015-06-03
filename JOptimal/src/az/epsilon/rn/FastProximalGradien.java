package az.epsilon.rn;

import az.epsilon.R1Function;
import az.epsilon.r1.R1Minimize;
import az.epsilon.RnFunction;
import java.util.Arrays;

public class FastProximalGradien implements Gradient {

    private double[] x;
    private double[] g;
    private double epsilon;
    private RnFunction f;

    @Override
    public void setFunction(RnFunction f) {
        this.f = f;
    }

    @Override
    public void gradient(double[] gradients) {
        double h = 0.000001;

        for (int i = 0; i < gradients.length; i++) {

            double _x = x[i];

            x[i] = x[i] + h;
            double f1 = f.f(x);
            x[i] = x[i] - 2 * h;
            double f2 = f.f(x);
            x[i] = x[i] + h;

            x[i] = _x;

            gradients[i] = (f1 - f2) / (2.0 * h);
        }
    }

    @Override
    public double argmin(double alpha) {
        double[] x1 = Arrays.copyOf(x, x.length);
        for (int i = 0; i < x.length; i++) {
            x1[i] = x[i] - alpha * g[i];
        }
        return f.f(x1);
    }

    @Override
    public double minimize() {

        R1Minimize minimize = new R1Minimize();
        minimize.setFunction(new R1Function() {
            @Override
            public double f(double alpha) {
                return argmin(alpha);
            }
        });
        minimize.setX(0.0);
        minimize.setStep(0.1);
        minimize.setEpsilon(0.00001);
        minimize.straightLineSearch();
        return minimize.goldenSectionSearch();
    }

    @Override
    public void calculate(double[] x0) {

        int n = x0.length;
        x = null;
        x = Arrays.copyOf(x0, n);
        g = new double[n];
        double[] x1 = new double[n];

        System.out.printf("x1 = %f x2 = %f\n", x[0], x[1]);

        double grad_norm = 0.0;
        do {

            gradient(g);

            grad_norm = norm(g);
            if (grad_norm < 0.00001) {
                break;
            }

            double alpha = minimize();

            x1 = Arrays.copyOf(x, n);
            for (int i = 0; i < x.length; i++) {
                x[i] = x[i] - alpha * g[i];
            }

            System.out.printf("x1 = %.10f x2 = %.10f\n", x[0], x[1]);

        } while (grad_norm > epsilon && distance(x1, x, n) > epsilon);
    }

    @Override
    public double norm(double[] x) {
        double g_norm = 0.0;
        for (int i = 0; i < x.length; i++) {
            g_norm = g_norm + x[i] * x[i];
        }
        g_norm = Math.sqrt(g_norm);
        return g_norm;
    }

    protected double distance(double[] x, double[] y, int n) {
        double dist = 0.0;
        for (int i = 0; i < x.length; i++) {
            dist = dist + (x[i] - y[i]) * (x[i] - y[i]);
        }
        dist = Math.sqrt(dist);
        return dist;
    }

    @Override
    public void setEpsilon(double epsilon) {
        this.epsilon = epsilon;
    }
}
