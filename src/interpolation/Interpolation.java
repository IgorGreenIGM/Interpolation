package interpolation;

import commons.Point;

import java.io.*;

public abstract class Interpolation {

    protected Point[] samplePoints;

    /**
     * Default constructor
     */
    Interpolation() {
        this.samplePoints = new Point[0];
    }

    /**
     * Constructs an Interpolation object with the provided sample points.
     *
     * @param samplePoints an array of Point objects representing sample points
     */
    Interpolation(Point[] samplePoints) {
        this.samplePoints = samplePoints;
    }

    /**
     * Constructs an Interpolation Object with the provided sample points from an input file.
     *
     * @param inputFilePath sample points input file path.
     */
    Interpolation(String inputFilePath) {

        BufferedReader buffr = null;
        BufferedReader buffrcpt = null;
        try {
            // compute and set the number of sample points inside the input file
            int nbPoints = 0;
            buffrcpt = new BufferedReader(new FileReader(inputFilePath));
            while (buffrcpt.readLine() != null) ++nbPoints;
            this.samplePoints = new Point[nbPoints];

            buffr = new BufferedReader(new FileReader(inputFilePath));
            int cpt = 0;
            while (true) {
                String[] xy = buffr.readLine().split(" ");
                this.samplePoints[cpt] = new Point(Double.parseDouble(xy[0]), Double.parseDouble(xy[1]));
                if (xy == null)
                    break;
            }

        } catch (IOException e) {
            System.out.println("An error occurred while creating the file: " + e.getMessage());
        } finally {
            try {
                assert buffr != null;

                buffr.close();
                buffrcpt.close();
            } catch (IOException ignored) {
            }
        }
    }

    /**
     * Constructs an Interpolation object with arrays of x and y values.
     *
     * @param xValues an array of x values
     * @param yValues an array of y values
     * @throws IllegalArgumentException if the sizes of xValues and yValues arrays are different
     */
    Interpolation(double[] xValues, double[] yValues) throws IllegalArgumentException {
        if (xValues.length != yValues.length)
            throw new IllegalArgumentException("Error, xValues array size must be the same as yValues");

        this.samplePoints = new Point[xValues.length];
        for (int i = 0; i < xValues.length; ++i)
            this.samplePoints[i] = new Point(xValues[i], yValues[i]);
    }

    public abstract double evaluate(double x);

    public void savePlottingDatas(double a, double b, int nbPoints, String path) {

        double range = b - a;
        FileWriter fw = null;
        File file = new File(path);
        try {
            file.delete();
            file.createNewFile();
            fw = new FileWriter(file);

            double x = a;
            fw.append(Double.toString(a)).append(" ").append(Double.toString(evaluate(a))).append("\n");
            for (int i = 1; i < nbPoints - 1; ++i) {
                x += (range / (double) (nbPoints - 1));
                fw.append(Double.toString(x)).append(" ").append(Double.toString(evaluate(x))).append("\n");
            }
            fw.append(Double.toString(b)).append(" ").append(Double.toString(evaluate(b))).append("\n");
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