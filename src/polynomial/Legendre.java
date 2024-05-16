package polynomial;

public class Legendre extends Polynomial {
    public Legendre(int degree) throws IllegalArgumentException {
        super();
        this.degree = degree;
    }

    /**
     * Generates the Legendre polynomial of degree 'n'.
     *
     * @param n The degree of the Legendre polynomial to generate.
     * @return The Legendre polynomial of degree 'n'.
     */
    public static Polynomial getLegendrePolynomial(int n) {
        Polynomial p0 = new Polynomial(0, new double[]{1.});
        if (n == 0)
            return p0;

        Polynomial p1 = new Polynomial(1, new double[]{0., 1.});
        if (n == 1)
            return p1;

        Polynomial x = new Polynomial(1, new double[]{0., 1.});
        for (int i = 2; i <= n; ++i) {
            Polynomial p = Polynomial.substract(Polynomial.multiply(Polynomial.multiply(p1, x), 2 * i - 1),
                                                Polynomial.multiply(p0, i - 1));
            p.multiplyInvert(i);
            p0 = p1;
            p1 = p;
        }

        return  p1;
    }

    /**
     * Calculate the roots of a polynomial using the Newton method,
     * starting from the last known roots.
     * This method iteratively refines the roots of the polynomial by applying
     * the Newton method, starting from the previously known roots.
     *
     * @param lastRoots The array containing the last known roots of the polynomial.
     * @param polynomial The polynomial for which roots are being calculated.
     * @return An array containing the refined roots of the polynomial.
     */
    private double[] getRootsFromPreviousRoots(double[] lastRoots, Polynomial polynomial) {
        double tol = 1e-5;
        double[] roots = new double[lastRoots.length + 1];
        double last = -1.0;
        double Xn; // variable to store the actual value of Newton suit
        double temp;

        for (int j = 0; j < lastRoots.length; j++) {
            Xn = (last + lastRoots[j]) / 2;
            temp = Polynomial.evaluate(Xn, polynomial);
            while ((Math.abs(temp)) > tol) {
                double derivative = polynomial.evaluateDerivative(Xn);
                if (derivative == 0) {
                    Xn = 0;
                    temp = polynomial.evaluate(Xn);
                    break;
                }
                Xn = Xn - temp / polynomial.evaluateDerivative(Xn);
                temp = polynomial.evaluate(Xn);
            }
            last = lastRoots[j];
            roots[j] = Xn;
            if (j == lastRoots.length - 1) {
                Xn = (1.0 + last) / 2;
                temp = polynomial.evaluate(Xn);

                while (Math.abs(temp) > tol) {
                    Xn = Xn - temp / polynomial.evaluateDerivative(Xn);
                    temp = polynomial.evaluate(Xn);
                }
                roots[j + 1] = Xn;
            }
        }
        return roots;
    }

    /**
     * Computes the roots of Legendre polynomials up to degree 'n'.
     * In this approach we progressively compute higher degree polynomials as we want.
     *
     * @return An array containing the roots of Legendre polynomials up to degree 'n'.
     */
    public double[] getRoots() {
        double[] roots = new double[]{0};

        if (this.degree != 1) {
            Polynomial p0 = new Polynomial(0, new double[]{1.});
            Polynomial x = new Polynomial(1, new double[]{0., 1.});
            Polynomial p1 = new Polynomial(1, new double[]{0., 1.});

            for (int i = 2; i <= this.degree; ++i) {
                // computing actual i-th legendre polynomial
                Polynomial p = Polynomial.substract(Polynomial.multiply(Polynomial.multiply(p1, x), 2 * i - 1),
                                                    Polynomial.multiply(p0, i - 1));
                p.multiplyInvert(i);

                // get new roots
                double[] nRoots = getRootsFromPreviousRoots(roots, p);
                roots = new double[nRoots.length];
                System.arraycopy(nRoots, 0, roots, 0, roots.length);

                // update further polynomials
                p0 = p1;
                p1 = p;
            }
        }

        return roots;
    }

    /**
     * This function calculates the weight of Gauss approximation from roots of Legendre polynomial
     *
     * @param roots roots of Legendre polynomial
     * @return array of weight
     */
    public double[] getWeight(double[] roots) {
        double[] weights = new double[roots.length];
        int n = roots.length;
        Polynomial der = getLegendrePolynomial(this.degree).derivative();
        for (int j = 0; j < n; j++) {
            double temp = Math.pow(der.evaluate(roots[j]), 2);
            weights[j] = 2 / ((1 - roots[j] * roots[j]) * temp);
        }

        return weights;
    }
}
