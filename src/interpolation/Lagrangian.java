package interpolation;

import java.io.FileReader;
import java.io.IOException;
import java.io.BufferedReader;

import commons.Point;
import polynomial.Polynomial;

/**
 * Provides methods for performing polynomial interpolation.
 */
public class Lagrangian extends Interpolation {

    /**
     * Constructs an Interpolation object with the provided sample points.
     *
     * @param samplePoints an array of Point objects representing sample points
     */
    public Lagrangian(Point[] samplePoints) {
        super(samplePoints);
    }

    /**
     * Constructs an Interpolation Object with the provided sample points from an input file.
     *
     * @param inputFilePath sample points input file path.
     */
    Lagrangian(String inputFilePath) {
        super(inputFilePath);
    }

    /**
     * Constructs an Interpolation object with arrays of x and y values.
     *
     * @param xValues an array of x values
     * @param yValues an array of y values
     * @throws IllegalArgumentException if the sizes of xValues and yValues arrays are different
     */
    Lagrangian(double[] xValues, double[] yValues) throws IllegalArgumentException {
        super(xValues, yValues);
    }

    /**
     * Evaluates the y value corresponding to the given x value using polynomial interpolation. Based on Neville Algorithm
     *
     * @param x the x value for interpolation
     * @return the interpolated y value
     */
    @Override
    public double evaluate(double x) {
        double[][] p = new double[this.samplePoints.length][this.samplePoints.length];

        for (int i = 0; i < this.samplePoints.length; ++i)
            p[i][0] = this.samplePoints[i].getY();

        for (int i = 1; i < this.samplePoints.length; ++i)
            for (int j = 1; j <= i; ++j)
                p[i][j] = ((x - this.samplePoints[i - j].getX()) * p[i][j - 1] - (x - this.samplePoints[i].getX()) * p[i - 1][j - 1]) / (this.samplePoints[i].getX() - this.samplePoints[i - j].getX());

        return p[this.samplePoints.length - 1][this.samplePoints.length - 1];
    }

    /**
     * Computes the polynomial coefficients for Lagrange interpolation.
     *
     * @return a Polynomial object representing the lagrangian polynomial
     */
    public Polynomial computePolynomial() {
        double[] tmpExpr = new double[this.samplePoints.length];
        double[] lagrangianCoeffs = new double[this.samplePoints.length];

        for (int i = 0; i < this.samplePoints.length; ++i) {
            double prod = 1.0;

            for (int j = 0; j < this.samplePoints.length; ++j)
                tmpExpr[j] = .0;

            for (int j = 0; j < this.samplePoints.length; ++j)
                if (j != i)
                    prod *= (this.samplePoints[i].getX() - this.samplePoints[j].getX());
            prod = this.samplePoints[i].getY() / prod;

            tmpExpr[0] = prod;
            for (int j = 0; j < this.samplePoints.length; ++j) {
                if (j != i)
                    for (int k = this.samplePoints.length - 1; k > 0; --k) {
                        tmpExpr[k] += tmpExpr[k - 1];
                        tmpExpr[k - 1] *= (-this.samplePoints[j].getX());
                    }
            }

            for (int j = 0; j < this.samplePoints.length; ++j)
                lagrangianCoeffs[j] += tmpExpr[j];
        }

        return new Polynomial(lagrangianCoeffs.length - 1, lagrangianCoeffs);
    }
}