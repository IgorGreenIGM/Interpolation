package commons;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.function.Function;

public class Point {

    private double x;
    private double y;

    public Point(double x, double y) {
        this.x = x;
        this.y = y;
    }

    public double getX() {
        return x;
    }

    public void setX(double x) {
        this.x = x;
    }

    public double getY() {
        return y;
    }

    public void setY(double y) {
        this.y = y;
    }

    /**
     * Generates the Chebyshev abscissa for the interval [a, b].
     *
     * @param a The lower bound of the interval.
     * @param b The upper bound of the interval.
     * @param n The number of points to generate.
     * @return An array containing the Chebyshev abscissa in the interval [a, b].
     */
    public static double[] tchebychevAbscissas(double a, double b, int n) {
        double[] abscissas = new double[n];

        double t = ((a + b) / 2.0);
        double tf = ((b - a) / 2.0);
        for (int i = 0; i < n; ++i)
            abscissas[i] = t + tf * Math.cos(((2 * i + 1) * Math.PI) / (2 * n));

        return abscissas;
    }

    /**
     * Generate Tchebyshev samples points for the interval [a, b]
     *
     * @param a The lower bound of the interval.
     * @param b The upper bound of the interval.
     * @param n The number of points to generate.
     * @param f The function used to evaluate Tchebyshev Abscissa
     * @return An array containing the uniform abscissa in the interval [a, b].
     */
    public static Point[] tchebychevSamples(double a, double b, int n, Function<Double, Double> f) {
        Point[] samples = new Point[n];

        double t = (a + b) / 2.0;
        double tf = (b - a) / 2.0;
        for (int i = 1; i <= n; ++i) {
            double x = t + tf * Math.cos(((2.0 * i - 1) * Math.PI) / (2.0 * n));
            samples[i - 1] = new Point(x, f.apply(x));
        }

        return samples;
    }

    /**
     * Generate uniform distributed abscissa for the interval [a, b]
     *
     * @param a The lower bound of the interval.
     * @param b The upper bound of the interval.
     * @param n The number of points to generate.
     * @return An array containing the uniform abscissa in the interval [a, b].
     */
    public static double[] uniformAbscissas(double a, double b, int n) {
        double x = a;
        double range = b - a;
        double[] sampleAbscissas = new double[n];
        sampleAbscissas[0] = a;
        sampleAbscissas[sampleAbscissas.length - 1] = b;

        for (int i = 1; i < n - 1; ++i) {
            x += (range / (n - 1));
            sampleAbscissas[i] = x;
        }

        return sampleAbscissas;
    }

    /**
     * Generate uniform distributed sample points for the interval [a, b]
     *
     * @param a The lower bound of the interval.
     * @param b The upper bound of the interval.
     * @param n The number of points to generate.
     * @param f The function used to evaluate generated abscissas
     * @return An array containing the uniform abscissa in the interval [a, b].
     */
    public static Point[] uniformSamples(double a, double b, int n, Function<Double, Double> f) {
        double x = a;
        double range = b - a;
        Point[] samplePoints = new Point[n];
        samplePoints[0] = new Point(a, f.apply(a));
        samplePoints[samplePoints.length - 1] = new Point(b, f.apply(b));

        for (int i = 1; i < n - 1; ++i) {
            x += (range / (n - 1));
            samplePoints[i] = new Point(x, f.apply(x));
        }

        return samplePoints;
    }

    public static void save(Point[] samplePoints, String path) throws RuntimeException {
        FileWriter fw = null;
        File file = new File(path);
        try {
            file.delete();
            file.createNewFile();
            fw = new FileWriter(file);
            for (Point p : samplePoints)
                fw.append(Double.toString(p.getX())).append(" ").append(Double.toString(p.getY())).append("\n");

        } catch (IOException e) {
            System.out.println("An error occurred while creating the file: " + e.getMessage());
        } finally {
            try {
                assert fw != null;
                fw.close();
            } catch (IOException ignored) {
            }
        }
    }
}
