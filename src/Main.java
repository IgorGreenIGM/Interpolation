import commons.Point;
import interpolation.Error;
import interpolation.Lagrangian;
import interpolation.CubicSpline;
import interpolation.LinearSpline;

import java.io.File;
import java.util.Arrays;
import java.io.FileWriter;
import java.io.IOException;
import java.io.BufferedWriter;
import java.util.function.Function;

public class Main {

    private static final int nbSamples = 9;
    private static final double lowerBound = -1.;
    private static final double upperBound = 1.;
    private static final int nbPlotPoints = 200;
    private static final int errorDegree = 10;

    public static double f(double x) {
        return 1.0 / ( 1 + (25*Math.pow(x, 2.)));
    }

    public static void main(String[] args) {
        // generating uniform sample points
        Point[] samples = Point.uniformSamples(lowerBound, upperBound, nbSamples, Main::f);

        // generating tchebyshev sample points
        Point[] tchebyshevSamples = Point.tchebychevSamples(lowerBound, upperBound, nbSamples, Main::f);

        // save generated sample Points
        Point.save(samples, "C:\\Users\\pc\\Desktop\\python\\Interpolation\\samples.txt");
        Point.save(tchebyshevSamples, "C:\\Users\\pc\\Desktop\\python\\Interpolation\\samplesTchebyshev.txt");

        // computing several interpolations
        Lagrangian lagrangian = new Lagrangian(samples);
        Lagrangian lagrangianTchebyshev = new Lagrangian(tchebyshevSamples);
        LinearSpline linearSpline = new LinearSpline(samples);
        CubicSpline cubicSpline = new CubicSpline(samples);

        // generate and save Plotting data
        lagrangian.savePlottingDatas(lowerBound, upperBound, nbPlotPoints, "C:\\Users\\pc\\Desktop\\python\\Interpolation\\lagrangian.txt");
        lagrangianTchebyshev.savePlottingDatas(lowerBound, upperBound, nbPlotPoints, "C:\\Users\\pc\\Desktop\\python\\Interpolation\\lagrangianTchebyshev.txt");
        linearSpline.savePlottingDatas(lowerBound, upperBound, nbPlotPoints, "C:\\Users\\pc\\Desktop\\python\\Interpolation\\linearSplines.txt");
        cubicSpline.savePlottingDatas(lowerBound, upperBound, nbPlotPoints, "C:\\Users\\pc\\Desktop\\python\\Interpolation\\cubicSplines.txt");

        // Generate and save function plotting data
        saveFunctionPlottingData(lowerBound, upperBound, nbPlotPoints, "C:\\Users\\pc\\Desktop\\python\\Interpolation\\function.txt", Main::f);

        // compute interpolation errors
        Error error = new Error(lowerBound, upperBound, errorDegree, Main::f);
        double lagrangianError = error.lagrangian(lagrangian.computePolynomial());
        double lagrangianTchebyshevError = error.lagrangian(lagrangianTchebyshev.computePolynomial());
        double linearSplineError = error.spline(linearSpline.getSplines(), Arrays.stream(samples).mapToDouble(Point::getX).toArray());
        double cubicSplineError = error.spline(cubicSpline.getSplines(), Arrays.stream(samples).mapToDouble(Point::getX).toArray());

        // save interpolations errors
        saveErrors(lagrangianError, lagrangianTchebyshevError, linearSplineError, cubicSplineError, "C:\\Users\\pc\\Desktop\\python\\Interpolation\\errors.txt");
    }

    private static void saveErrors(double lagrangian, double lagrangianTchebyshev, double linearSpline, double cubicSpline, String path) {
        BufferedWriter bw = null;
        try {
            bw = new BufferedWriter(new FileWriter(path));
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        try {
            bw.write("Lagrangian Error : " + lagrangian + "\n");
            bw.write("Lagrangian tchebyshev Error : " + lagrangianTchebyshev + "\n");
            bw.write("Linear spline Error : " + linearSpline + "\n");
            bw.write("Cubic Spline Error : " + cubicSpline + "\n");
            bw.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public static void saveFunctionPlottingData(double a, double b, int nbPoints, String path, Function<Double, Double> f) {

        double range = b - a;
        FileWriter fw = null;
        File file = new File(path);
        try {
            file.delete();
            file.createNewFile();
            fw = new FileWriter(file);

            double x = a;
            fw.append(Double.toString(a)).append(" ").append(Double.toString(f.apply(a))).append("\n");
            for (int i = 1; i < nbPoints - 1; ++i) {
                x += (range / nbPoints);
                fw.append(Double.toString(x)).append(" ").append(Double.toString(f.apply(x))).append("\n");
            }
            fw.append(Double.toString(b)).append(" ").append(Double.toString(f.apply(b))).append("\n");
        } catch (IOException e) {
            System.err.println("An error occurred while creating the file: " + e.getMessage());
        } finally {
            try {
                assert fw != null;
                fw.close();
            } catch (IOException ignored) {
            }
        }
    }
}