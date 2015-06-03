package az.epsilon;

public class Main {

    public static void main(String[] args) {
        
        double[] x0 = {-1.0, +1.2};
        
        Gradient g = new Gradient1();
        g.calculate(x0);
    }
}
