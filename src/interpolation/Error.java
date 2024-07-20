package interpolation;

import polynomial.Legendre;
import polynomial.Polynomial;

import java.util.function.Function;

/**
 * Represents the error calculation methods for polynomial interpolation.
 */
public class Error {
    private final int degree;
    private final double lowerBound;
    private final double upperBound;
    private final double[] legendreRoots;
    private final double[] legendreWeights;
    private final Function<Double, Double> f;

    /**
     * Constructs an Error object.
     *
     * @param lowerBound  the lower bound of the interval
     * @param upperBound  the upper bound of the interval
     * @param degree      the degree of the polynomial
     * @param f           the function to interpolate
     */
    public Error(double lowerBound, double upperBound, int degree, Function<Double, Double> f) {
        this.f = f;
        this.degree = degree;
        this.lowerBound = lowerBound;
        this.upperBound = upperBound;

        Legendre legendre = new Legendre(degree);
        this.legendreRoots = legendre.getRoots();
        this.legendreWeights = legendre.getWeight(this.legendreRoots);

        for (Double x : this.legendreWeights)
            System.out.print(x + " ");
        System.out.println();
    }

    /**
     * Calculates the Lagrangian error of the interpolation.
     *
     * @param p the interpolating polynomial
     * @return the Lagrangian error
     */
    public double lagrangian(Polynomial p) {
        double sum = .0;
        double m = (lowerBound + upperBound) / 2.;
        double sr = (upperBound - lowerBound) / 2.;
        for (int  i = 0; i < legendreRoots.length; ++i) {
            double toEval = sr*legendreRoots[i] + m;
            sum += legendreWeights[i]*Math.pow(f.apply(toEval) - p.evaluate(toEval), 2.);
        }

        return  Math.sqrt(sr*sum);
    }

    /**
     * Calculates any spline(cubic or linear) interpolation error.
     *
     * @param polynomials the array of interpolating polynomials
     * @param partition   the array representing the partition of the interval
     * @return the linear error
     */
    public double spline(Polynomial[] polynomials, double[] partition) {
        double sum = .0;
        double m = (lowerBound + upperBound) / 2.;
        double sr = (upperBound - lowerBound) / 2.;

        for (int i = 0; i < degree; ++i) {
            for (int j = 0; j < Math.min(partition.length, polynomials.length) - 1; ++j) {
                double toEval = sr*legendreRoots[i] + m;
                if (toEval >= partition[j] && toEval <= partition[j + 1]) {
                    sum += legendreWeights[i]*Math.pow(f.apply(toEval) - polynomials[j].evaluate(toEval), 2.);
                    break;
                }
            }
        }

        return Math.sqrt(sr*sum);
    }


    /* getters */

    public Function<Double, Double> getF() {
        return f;
    }

    public double getLowerBound() {
        return lowerBound;
    }

    public double getUpperBound() {
        return upperBound;
    }
}