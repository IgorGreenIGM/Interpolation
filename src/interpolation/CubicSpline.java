package interpolation;

import commons.Point;
import polynomial.Polynomial;

import java.util.Arrays;
import java.util.Comparator;

public class CubicSpline extends Interpolation {
    private final Polynomial[] splines;

    /**
     * Constructs a CubicSpline object with the given sample values and corresponding sample data.
     *
     * @param xValues an array of sample values representing the x-coordinates of the data points.
     * @param yValues an array of sample data representing the y-coordinates of the data points.
     * @throws IllegalArgumentException if the length of sampleValue is less than 2,
     *                                  or if the lengths of sampleValue and sampleData are different.
     */
    public CubicSpline(double[] xValues, double[] yValues) throws IllegalArgumentException {
        super(xValues, yValues);
        Arrays.sort(this.samplePoints, Comparator.comparingDouble(Point::getX));
        this.splines = this.getSplines();
    }

    /**
     * Constructs a CubicSpline object with the given array of Point samples.
     *
     * @param samples an array of Point objects representing the data points.
     */
    public CubicSpline(Point[] samples) {
        super(samples);
        Arrays.sort(this.samplePoints, Comparator.comparingDouble(Point::getX));
        this.splines = this.getSplines();
    }

    /**
     * Constructs an CubicSpline Object with the provided sample points from an input file.
     *
     * @param inputFilePath input file path
     */
    CubicSpline(String inputFilePath) {
        super(inputFilePath);
        Arrays.sort(this.samplePoints, Comparator.comparingDouble(Point::getX));
        this.splines = this.getSplines();
    }

    /**
     * Computes the height between adjacent sample points along the x-axis.
     *
     * @return An array containing the height between adjacent sample points.
     */
    private double[] heightPoint() {
        int n = this.samplePoints.length - 1;
        double[] height = new double[n];
        for (int i = 0; i < n; i++) {
            height[i] = this.samplePoints[i + 1].getX() - this.samplePoints[i].getX();
        }

        return height;
    }

    /**
     * Computes the sum of heights between adjacent sample points.
     *
     * @param height An array containing the height between adjacent sample points.
     * @return An array containing the sum of heights between adjacent sample points.
     */
    private double[] heightSum(double[] height) {
        double[] heiSum = new double[height.length - 1];
        for (int i = 0; i < height.length - 1; i++) {
            heiSum[i] = height[i] + height[i + 1];
        }

        return heiSum;
    }

    /**
     * Constructs a tridiagonal matrix for cubic spline interpolation.
     *
     * @param heiSum An array containing the sum of heights between adjacent sample points.
     * @param height An array containing the height between adjacent sample points.
     * @return A 2D array representing the constructed tridiagonal matrix.
     * @throws IllegalArgumentException If the sizes of the arrays do not match.
     */
    public double[][] matrix(double[] heiSum, double[] height) throws IllegalArgumentException {
        if (heiSum.length != height.length - 1) {
            throw new IllegalArgumentException("Sizes of table passed do not match with data. They should be equal");
        }

        double[][] matrix = new double[2][heiSum.length];
        matrix[1][0] = 2 * heiSum[0];
        for (int i = 1; i < heiSum.length - 1; i++) {
            matrix[0][i] = height[i];
            matrix[1][i] = 2 * heiSum[i];
        }
        matrix[0][heiSum.length - 1] = height[heiSum.length - 1];
        matrix[1][heiSum.length - 1] = 2 * heiSum[heiSum.length - 1];

        return matrix;
    }

    /**
     * Computes the constant vector for cubic spline interpolation matrix system resolution.
     *
     * @param height An array containing the height between adjacent sample points.
     * @return An array representing the constant vector for cubic spline interpolation.
     * @throws IllegalArgumentException If the size of the height array does not match the expected size.
     */
    public double[] constantVector(double[] height) throws IllegalArgumentException {
        int n = this.samplePoints.length;
        if (height.length != n - 1) {
            throw new IllegalArgumentException("Image and height tables sizes do not match with data");
        }

        double[] constVector = new double[n - 2];
        double[] beta = new double[n - 1];
        beta[0] = (this.samplePoints[1].getY() - this.samplePoints[0].getY()) / height[0];
        for (int i = 1; i < n - 1; i++) {
            beta[i] = (this.samplePoints[i + 1].getY() - this.samplePoints[i].getY()) / height[i];
            constVector[i - 1] = 6*(beta[i] - beta[i - 1]);
        }

        return constVector;
    }

    /**
     * Solves a tridiagonal symmetric matrix equation.
     * This method solves a tridiagonal symmetric matrix equation using the Thomas algorithm.
     *
     * @param sub The sub-diagonal elements of the tridiagonal matrix.
     * @param diag The diagonal elements of the tridiagonal matrix.
     * @param constVector The constant vector in the matrix equation.
     * @return The solution vector for the matrix equation.
     */
    private double[] solveTridiagonalSymMatrix(double[] sub, double[] diag, double[] constVector) {
        double[] x = new double[diag.length];
        double[] cp = new double[sub.length];
        double[] dp = new double[sub.length];

        cp[0] = sub[1] / diag[0];
        for (int i = 1; i < diag.length - 1; ++i) {
            cp[i] = sub[i + 1] / (diag[i] - cp[i - 1] * sub[i]);
        }

        dp[0] = constVector[0] / diag[0];
        for (int i = 1; i < diag.length; ++i) {
            dp[i] = (constVector[i] - dp[i - 1] * sub[i]) / (diag[i] - cp[i - 1] * sub[i]);
        }

        x[x.length - 1] = dp[x.length - 1];
        for (int i = x.length - 2; i >= 0; --i) {
            x[i] = dp[i] - cp[i] * x[i + 1];
        }

        return x;
    }

    /**
     * Computes the second derivatives at each sample point for cubic spline interpolation.
     *
     * @param matrix A 2D array representing the tridiagonal matrix used in cubic spline interpolation.
     * @param constVector An array representing the constant vector used in cubic spline interpolation.
     * @return An array containing the second derivatives at each sample point.
     */
    private double[] secondDerivativePoint(double[][] matrix, double[] constVector) {
        double[] secDeriv = new double[matrix[0].length + 2];
        double[] result = solveTridiagonalSymMatrix(matrix[0], matrix[1], constVector);
        System.arraycopy(result, 0, secDeriv, 1, matrix[0].length);

        return secDeriv;
    }

    /**
     * Computes the coefficients of cubic spline interpolation.
     *
     * @param height An array containing the height between adjacent sample points.
     * @param secDeriv An array containing the second derivatives at each sample point.
     * @return A 2D array containing the coefficients of cubic spline interpolation for each interval.
     * @throws IllegalArgumentException If the size of the height array does not match the expected size.
     */
    private double[][] coefficients(double[] height, double[] secDeriv) throws IllegalArgumentException {
        int n = this.samplePoints.length;
        if (height.length != n - 1) {
            throw new IllegalArgumentException("The size of magnitude table does not match with data");
        }

        double[][] coefficients = new double[n - 1][4];
        for (int i = 0; i < n - 1; i++) {
            coefficients[i][0] = this.samplePoints[i].getY();
            coefficients[i][1] = (this.samplePoints[i + 1].getY() - this.samplePoints[i].getY()) / height[i] - (2 * secDeriv[i] + secDeriv[i + 1]) * height[i] / 6;
            coefficients[i][2] = secDeriv[i];
            coefficients[i][3] = (secDeriv[i + 1] - secDeriv[i]) / height[i];
        }

        return coefficients;
    }

    /**
     * Converts coefficients of a cubic polynomial to the canonical base.
     *
     * @param coeffs An array containing the coefficients of a cubic polynomial.
     * @param a The value at which the polynomial is evaluated.
     * @return An array containing the coefficients of the polynomial in the canonical base.
     * @throws IllegalArgumentException If the size of the coefficients array is not 4.
     */
    private double[] toCanonicalBase(double[] coeffs, double a) throws IllegalArgumentException{
        if (coeffs.length != 4)
            throw new IllegalArgumentException("only coefficients array of degree-3(array length = 4) polynomials are accepted");

        double[] nCoeffs = new double[4];
        nCoeffs[0] = coeffs[0] + (coeffs[1]*a) + (0.5*coeffs[2]*Math.pow(a, 2.0)) + ((coeffs[3]*Math.pow(a, 3.0)) / 6.0);
        nCoeffs[1] = coeffs[1] + (a*coeffs[2]) + (.5*Math.pow(a, 2.0)*coeffs[3]);
        nCoeffs[2] = .5*coeffs[2] + (.5*a*coeffs[3]);
        nCoeffs[3] = coeffs[3]/6.;

        return nCoeffs;
    }

    /**
     * Computes cubic spline polynomials for interpolation.
     *
     * @return An array of cubic spline polynomials for interpolation.
     */
    public Polynomial[] getSplines() {
        double[] height = heightPoint();
        double[][] coeffs = coefficients(height, secondDerivativePoint(matrix(heightSum(height), height), constantVector(height)));

        Polynomial[] result = new Polynomial[this.samplePoints.length - 1];
        for (int i = 0; i < this.samplePoints.length - 1; i++) {
            result[i] = new Polynomial(3, toCanonicalBase(coeffs[i], -this.samplePoints[i].getX()));
        }

        return result;
    }

    /**
     * Evaluates the cubic spline interpolation at a given x-coordinate.
     *
     * @param x The x-coordinate at which to evaluate the interpolation.
     * @return The interpolated y-coordinate corresponding to the given x-coordinate.
     * @throws IllegalArgumentException If there are no sample points or if the provided x-coordinate
     *         is outside the range of the sample points.
     */
    @Override
    public double evaluate(double x) throws IllegalArgumentException {

        if (this.samplePoints.length == 0) {
            throw new IllegalArgumentException("You need  to provide at least one point for evaluation.");
        }

        if (x < this.samplePoints[0].getX() || x > this.samplePoints[this.samplePoints.length - 1].getX()) {
            throw new IllegalArgumentException("The value should be in the correct interval");
        }

        int sub = -1;
        for (int i = 0; i < this.samplePoints.length - 1; ++i) {
            if (x >= samplePoints[i].getX() && x <= samplePoints[i + 1].getX()) {
                sub = i;
                break;
            }
        }

        return this.splines[sub].evaluate(x);
    }
}