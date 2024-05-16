package interpolation;

import commons.Point;
import polynomial.Polynomial;

import java.util.Arrays;
import java.util.Comparator;

public class LinearSpline extends Interpolation {
    /**
     * Constructs a LinearSpline object with the given sample values and corresponding sample data.
     *
     * @param xValues an array of sample values representing the x-coordinates of the data points.
     * @param yValues  an array of sample data representing the y-coordinates of the data points.
     * @throws IllegalArgumentException if the length of sampleValue is less than 2,
     *                                  or if the lengths of sampleValue and sampleData are different.
     */

    public LinearSpline(double[] xValues, double[] yValues) throws IllegalArgumentException {
        super(xValues, yValues);
        Arrays.sort(this.samplePoints, Comparator.comparingDouble(Point::getX));
    }

    /**
     * Constructs a LinearSpline object with the given array of Point samples.
     *
     * @param samples an array of Point objects representing the data points.
     */
    public LinearSpline(Point[] samples) {
        super(samples);
        Arrays.sort(this.samplePoints, Comparator.comparingDouble(Point::getX));
    }

    /**
     * Constructs an LinearSpline Object with the provided sample points from an input file.
     *
     * @param inputFilePath input file path
     */
    LinearSpline(String inputFilePath) {
        super(inputFilePath);
        Arrays.sort(this.samplePoints, Comparator.comparingDouble(Point::getX));
    }

    /**
     * Evaluate the interpolated function
     *
     * @throws IllegalArgumentException if the number is not within the range of points represented
     */
    @Override
    public double evaluate(double x) throws IllegalArgumentException {
        Point maxInterval = null;
        Point minInterval = null;

        if (this.samplePoints.length == 0)
            throw new IllegalArgumentException("You need to define sample points of at least two points to obtain an approximation");

        if (x < this.samplePoints[0].getX() || x > this.samplePoints[this.samplePoints.length - 1].getX())
            throw new IllegalArgumentException("This number does not belong to the sampling interval");

        // Find the interval containing the given x
        for (int i = 0; i < samplePoints.length - 1; i++) {
            if (x >= samplePoints[i].getX() && x <= samplePoints[i + 1].getX()) {
                minInterval = this.samplePoints[i];
                maxInterval = this.samplePoints[i + 1];
                break;
            }
        }

        return minInterval.getY() + ((maxInterval.getY() - minInterval.getY()) / (maxInterval.getX() - minInterval.getX())) * (x - minInterval.getX());
    }

    /**
     * Calculates linear splines between consecutive points and returns an array of results.
     *
     * @return an array of results containing the linear splines polynomials.
     */
    public Polynomial[] getSplines() {
        Polynomial[] results = new Polynomial[this.samplePoints.length - 1];
        for (int i = 0; i < this.samplePoints.length - 1; ++i) {
            double m = (this.samplePoints[i + 1].getY() - this.samplePoints[i].getY()) / (this.samplePoints[i + 1].getX() - this.samplePoints[i].getX());
            double b = this.samplePoints[i].getY() - (m * this.samplePoints[i].getX());

            results[i] = new Polynomial(1, new double[]{b, m});
        }

        return results;
    }

    @Override
    public String toString() {
        StringBuilder stringBuilder = new StringBuilder();
        for (Point sample : this.samplePoints) {
            stringBuilder.append("x: ").append(sample.getX()).append(" y: ").append(sample.getY()).append("\n");
        }
        return stringBuilder.toString();
    }
}