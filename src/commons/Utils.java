package commons;

/**
 * Utility class providing methods.
 */
public class Utils {

    /**
     * Converts an integer to its upper-script representation.
     *
     * @param n The integer to be converted.
     * @return The upper-script representation of the integer.
     */
    public static String toUpperScript(int n) {

        String digits = "0123456789";
        String superScriptDigits = "⁰¹²³⁴⁵⁶⁷⁸⁹";

        String numberString = Integer.toString(n);
        StringBuilder output = new StringBuilder();
        for (int i = 0; i < numberString.length(); ++i)
            output.append(superScriptDigits.charAt(digits.indexOf(numberString.charAt(i))));

        return output.toString();
    }
}