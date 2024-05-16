package polynomial;

import commons.Utils;

/**
 * Represents a polynomial with integer degree and double coefficients.
 * The polynomial is represented by its degree and an array of coefficients.
 */
public class Polynomial {

    protected int degree;             // the degree of the polynomial
    protected double[] coefficients;  // the array of coefficients

    /**
     * Default constructor
     */
    public Polynomial() {
        this.degree = -1;
        this.coefficients = new double[0];
    }

    /**
     * Constructs a polynomial with the given degree and coefficients.
     *
     * @param degree       The degree of the polynomial.
     * @param coefficients The coefficients of the polynomial.
     * @throws IllegalArgumentException if the degree does not match the size of the coefficients array.
     */
    public Polynomial(int degree, double[] coefficients) throws IllegalArgumentException {
        if (degree != coefficients.length - 1)
            throw new IllegalArgumentException(
                    "coefficients array size doesn't matches the degree. Check if degree == coefficients.length - 1");

        this.degree = degree;
        this.coefficients = coefficients;
    }

    /**
     * Constructs a polynomial with the given degree and initializes coefficients to zero.
     *
     * @param degree The degree of the polynomial.
     */
    public Polynomial(int degree) {
        this.degree = degree;
        this.coefficients = new double[degree + 1];
    }

    /**
     * Copy constructor
     * @param p polynomial to copy
     */
    public Polynomial(Polynomial p) {
        this.degree = p.degree;
        this.coefficients = new double[p.coefficients.length];
        System.arraycopy(p.coefficients, 0, this.coefficients, 0, this.coefficients.length);
    }

    public static void copy(Polynomial src, Polynomial dest) {
        src.degree = dest.degree;
        src.coefficients = new double[dest.coefficients.length];
        System.arraycopy(src.coefficients, 0, dest.coefficients, 0, src.coefficients.length);
    }

    /**
     * Negates all coefficients of the polynomial.
     */
    public void negate() {
        for (int i = 0; i < this.degree + 1; ++i)
            this.coefficients[i] = -this.coefficients[i];
    }

    /**
     * Returns a new polynomial which is the negation of the given polynomial.
     *
     * @param p The polynomial to be negated.
     * @return The negated polynomial.
     */
    public static Polynomial negate(Polynomial p) {
        double[] newCoeffs = new double[p.degree + 1];
        for (int i = 0; i < p.degree + 1; ++i)
            newCoeffs[i] = -p.coefficients[i];

        return new Polynomial(p.degree, newCoeffs);
    }

    /**
     * Adds another polynomial to this polynomial.
     *
     * @param other The polynomial to be added.
     */
    public void add(Polynomial other) {
        if (this.degree >= other.degree) {
            for (int i = 0; i < other.degree + 1; ++i)
                this.coefficients[i] += other.coefficients[i];
        } else {
            double[] newCoeffs = new double[other.degree + 1];
            System.arraycopy(other.coefficients, 0, newCoeffs, 0, other.degree + 1);
            for (int i = 0; i < this.degree + 1; ++i)
                newCoeffs[i] += this.coefficients[i];

            this.degree = other.degree;
            this.coefficients = newCoeffs;
        }
    }

    /**
     * Returns a new polynomial which is the sum of the given polynomials.
     *
     * @param p1 The first polynomial.
     * @param p2 The second polynomial.
     * @return The sum of the two polynomials.
     */
    public static Polynomial add(Polynomial p1, Polynomial p2) {
        int newDegree = Math.max(p1.degree, p2.degree);
        double[] newCoeffs = new double[newDegree + 1];

        if (p1.degree >= p2.degree) {
            System.arraycopy(p1.coefficients, 0, newCoeffs, 0, p1.degree + 1);
            for (int i = 0; i < p2.degree + 1; ++i)
                newCoeffs[i] += p2.coefficients[i];
        } else {
            System.arraycopy(p2.coefficients, 0, newCoeffs, 0, p2.degree + 1);
            for (int i = 0; i < p1.degree + 1; ++i)
                newCoeffs[i] += p1.coefficients[i];
        }

        return new Polynomial(newDegree, newCoeffs);
    }

    /**
     * Subtracts another polynomial from this polynomial.
     *
     * @param other The polynomial to be subtracted.
     */
    public void substract(Polynomial other) {
        if (this.degree >= other.degree) {
            for (int i = 0; i < other.degree + 1; ++i)
                this.coefficients[i] -= other.coefficients[i];
        } else {
            double[] newCoeffs = new double[other.degree + 1];
            System.arraycopy(this.coefficients, 0, newCoeffs, 0, this.degree + 1);
            for (int i = 0; i < other.degree + 1; ++i)
                newCoeffs[i] -= this.coefficients[i];

            this.degree = other.degree;
            this.coefficients = newCoeffs;
        }
    }

    /**
     * Returns a new polynomial which is the difference of the given polynomials.
     *
     * @param p1 The first polynomial.
     * @param p2 The second polynomial.
     * @return The difference of the two polynomials.
     */
    public static Polynomial substract(Polynomial p1, Polynomial p2) {
        int newDegree = Math.max(p1.degree, p2.degree);
        double[] newCoeffs = new double[newDegree + 1];

        System.arraycopy(p1.coefficients, 0, newCoeffs, 0, p1.degree + 1);
        for (int i = 0; i < p2.degree + 1; ++i)
            newCoeffs[i] -= p2.coefficients[i];

        return new Polynomial(newDegree, newCoeffs);
    }

    /**
     * Returns a new polynomial which is the product of the given polynomials.
     *
     * @param other the polynomial to multiply with.
     */
    public void multiply(Polynomial other) {
        double[] newCoeffs = new double[this.degree + other.degree + 1];

        for (int i = 0; i < this.degree + 1; ++i)
            for (int j = 0; j < other.degree + 1; ++j)
                newCoeffs[i + j] += this.coefficients[i] * other.coefficients[j];

        this.degree += other.degree;
        this.coefficients = newCoeffs;
    }

    /**
     * Returns a new polynomial which is the product of the given polynomials.
     *
     * @param p1 The first polynomial.
     * @param p2 The second polynomial.
     */
    public static Polynomial multiply(Polynomial p1, Polynomial p2) {
        int newDegree = p1.degree + p2.degree;
        double[] newCoeffs = new double[newDegree + 1];

        for (int i = 0; i < p1.degree + 1; ++i)
            for (int j = 0; j < p2.degree + 1; ++j)
                newCoeffs[i + j] += p1.coefficients[i] * p2.coefficients[j];

        return new Polynomial(newDegree, newCoeffs);
    }

    /**
     * Multiply the actual polynomial by a double
     * @param x input double
     */
    public void multiply(double x) {
        for (int i = 0; i < this.coefficients.length; ++i) {
            this.coefficients[i] *= x;
        }
    }

    /**
     * Returns a new polynomial which is the product of a polynomial and a double
     * @param p input polynomial
     * @param x input double
     * @return Polynomial
     */
    public static Polynomial multiply(Polynomial p, double x) {
        double[] newCoeffs = new double[p.degree + 1];

        System.arraycopy(p.coefficients, 0, newCoeffs, 0, p.degree + 1);
        for (int i = 0; i < p.degree + 1; ++i)
            newCoeffs[i] *= x;

        return new Polynomial(p.degree, newCoeffs);
    }

    /**
     * divide the actual polynomial by a double
     * @param x input double
     */
    public void multiplyInvert(double x) {
        for (int i = 0; i < this.coefficients.length; ++i) {
            this.coefficients[i] /= x;
        }
    }

    /**
     * Returns a new polynomial which is the division of a polynomial and a double
     * @param p input polynomial
     * @param x input double
     * @return Polynomial
     */
    public static Polynomial multiplyInvert(Polynomial p, double x) {
        double[] newCoeffs = new double[p.degree + 1];

        System.arraycopy(p.coefficients, 0, newCoeffs, 0, p.degree + 1);
        for (int i = 0; i < p.degree + 1; ++i)
            newCoeffs[i] /= x;

        return new Polynomial(p.degree, newCoeffs);
    }

    /**
     * Evaluates the polynomial at a given value of x. Using horner scheme.
     *
     * @param x The value at which the polynomial is evaluated.
     * @return The result of evaluating the polynomial at x.
     */
    public double evaluate(double x) {
        double y = this.coefficients[this.degree];

        for (int i = this.degree - 1; i >= 0; --i)
            y = x * y + this.coefficients[i];
        return y;
    }

    /**
     * Evaluates the value of a polynomial at a given point.
     *
     * @param x The point at which to evaluate the polynomial.
     * @param polynomial The polynomial to evaluate.
     * @return The value of the polynomial at the specified point 'x'.
     */
    public static double evaluate(double x, Polynomial polynomial) {
        double y = polynomial.coefficients[polynomial.degree];
        for (int i = polynomial.degree - 1; i >= 0; --i)
            y = x * y + polynomial.coefficients[i];
        return y;
    }

    /**
     * Derivate the actual polynomial
     *
     * @return Polynomial
     */
    public Polynomial derivative() {
        double[] newCoefficients = new double[this.degree];
        for (int i = 0; i < this.degree; ++i) {
            newCoefficients[i] = (i + 1) * this.coefficients[i + 1];
        }
        return new Polynomial(this.degree - 1, newCoefficients);
    }

    /**
     * Derivate input polynomial
     *
     * @param polynomial input polynomial
     * @return derivated polynomial
     */
    public static Polynomial derivative(Polynomial polynomial) {
        double[] newCoefficients = new double[polynomial.degree];
        for (int i = 0; i < polynomial.degree; ++i) {
            newCoefficients[i] = (i + 1) * polynomial.coefficients[i];
        }
        return new Polynomial(polynomial.degree - 1, newCoefficients);
    }

    /**
     * Evaluate the derivate of this polynomial at a certain point.
     *
     * @param x point where to evaluate
     * @return double
     */
    public double evaluateDerivative(double x) {
        Polynomial polynomial = this.derivative();
        return polynomial.evaluate(x);
    }

    /**
     * Evaluate the derivate of the input polynomial at a certain point.
     *
     * @param x point where to evaluate
     * @return double
     */
    public static double evaluateDerivative(double x, Polynomial polynomial) {
        Polynomial polynomial1 = derivative(polynomial);
        return polynomial1.evaluate(x);
    }

    /* Setters and Getters */

    public int getDegree() {
        return degree;
    }

    public void setDegree(int degree) {
        this.degree = degree;
    }

    public double[] getCoefficients() {
        return coefficients;
    }

    public void setCoefficients(double[] coefficients) {
        this.coefficients = coefficients;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append(this.coefficients[0]).append(" + ");
        for (int i = 1; i < this.degree + 1; ++i) {
            sb.append(Double.toString(this.coefficients[i]));
            sb.append("x");
            sb.append(Utils.toUpperScript(i));
            sb.append(" + ");
        }

        return sb.substring(0, sb.length() - 2); // deleting last '+ '.
    }
}