package az.epsilon.main;

import az.epsilon.RnFunction;
import az.epsilon.rn.Gradient;
import az.epsilon.rn.FastProximalGradient;

public class Main {

    public static void main(String[] args) {

        double[] x0 = {-1.0, +1.2};
        Gradient g = new FastProximalGradient();
        g.setFunction(new RnFunction() {
            @Override
            public double f(double[] x) {
                double x1 = x[0];
                double x2 = x[1];
                return ((1 - x1) * (1 - x1)) + 100 * (x2 - x1 * x1) * (x2 - x1 * x1);
            }
        });
        g.setEpsilon(0.000001);
        g.calculate(x0);

//        R1Minimize min = new R1Minimize();
//        min.setFunction(new R1Function() {
//            @Override
//            public double f(double x) {
//                return x * x;
//            }
//        });
//        min.setEpsilon(0.000001);
//        min.setStep(0.1);
//        min.setX(2.0);
//        min.straightLineSearch();
//        System.out.println(min.getA() + ", " + min.getB() + ", " + min.goldenSectionSearch());
    }
}
