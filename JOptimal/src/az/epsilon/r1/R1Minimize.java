package az.epsilon.r1;

import az.epsilon.R1Function;

public class R1Minimize {

    private double x;
    private double epsilon;
    private double step;
    private double a;
    private double b;
    private R1Function f;

    public void setFunction(R1Function f) {
        this.f = f;
    }

    public double straightLineSearch() {

        if (step == 0.0) {
            return Double.NaN;
        }

        double y0 = f.f(x);
        double y1 = f.f(x - step);
        double y2 = f.f(x + step);

        // if y1 and y2 are both greater than y0 then minimum point is inside x1 and x2
        if (y1 >= y0 && y0 <= y2) {
            a = x - Math.abs(step);
            b = x + Math.abs(step);
            double c = (a + b) / 2.0;
            return c;
        }

        // if y1 and y2 are both lesser than y0 then there is not minimum. function is not unimodal
        if (y1 <= y0 && y0 >= y2) {
            System.err.println("Function is not unimodal");
            return Double.NaN;
        }

        if (y1 >= y0 && y0 >= y2) {
            while (y2 < y0) {
                x = x + Math.abs(step);
                y1 = y0;
                y0 = y2;
                y2 = f.f(x + Math.abs(step));
            }
            a = x - Math.abs(step);
            b = x + Math.abs(step);
            double c = (a + b) / 2.0;
            return c;
        }

        if (y1 <= y0 && y0 <= y2) {
            while (y1 < y0) {
                x = x - Math.abs(step);
                y2 = y0;
                y0 = y1;
                y1 = f.f(x - Math.abs(step));
            }
            a = x - Math.abs(step);
            b = x + Math.abs(step);
            double c = (a + b) / 2.0;
            return c;
        }

        a = x - Math.abs(step);
        b = x + Math.abs(step);
        double c = (a + b) / 2.0;
        return c;
    }

    public double goldenSectionSearch() {
        //double sqrt_5 = 2.2360679774997896964091736687313
        //double phi = (sqrt(5) + 1.0) / 2.0;
        //double phi = (sqrt(5) - 1) / 2.0;
        double phi = 1.6180339887498948482045868343656;

        double x1 = Double.NaN;
        double x2 = Double.NaN;

        double y1 = 0.0;
        double y2 = 0.0;

        // Lazimi epsilon deqiqliyini alana qeder
        // iterasiyalari davam edirik
        while (Math.abs(b - a) > epsilon) {
            if (Double.isNaN(x1)) {
                x1 = b - Math.abs(b - a) / phi;
                y1 = f.f(x1);
            }

            if (Double.isNaN(x2)) {
                x2 = a + Math.abs(b - a) / phi;
                y2 = f.f(x2);
            }

            if (y1 >= y2) {
                a = x1;
                x1 = x2;    // Tapilmish x2 noqtesi ve bu noqtede funkisiyanin qiymeti
                y1 = y2;    // sonraki iterasiyada x1 qiymeti kimi istifade olunacaq.
                x2 = Double.NaN;   // x2 novbeti iterasiyada axtarilacaq
            } else {
                b = x2;
                x2 = x1;    // Tapilmish x1 noqtesi ve bu noqtede funkisiyanin qiymeti
                y2 = y1;    // sonraki iterasiyada x2 qiymeti kimi istifade olunacaq.
                x1 = Double.NaN;   // x1 novbeti iterasiyada axtarilacaq
            }
        }

        double c = (a + b) / 2.0;

        if (f.f(a) < f.f(b)) {
            c = a;
        }
        if (f.f(a) > f.f(b)) {
            c = b;
        }

        return c;
    }

    public double halphIntervalMethod() {
        return 0.0;
    }

    public double getX() {
        return x;
    }

    public void setX(double x) {
        this.x = x;
    }

    public double getEpsilon() {
        return epsilon;
    }

    public void setEpsilon(double epsilon) {
        this.epsilon = epsilon;
    }

    public double getStep() {
        return step;
    }

    public void setStep(double step) {
        this.step = step;
    }

    public double getA() {
        return a;
    }

    public void setA(double a) {
        this.a = a;
    }

    public double getB() {
        return b;
    }

    public void setB(double b) {
        this.b = b;
    }
}
