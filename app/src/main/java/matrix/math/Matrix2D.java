//todo>> tests
package matrix.math;

import java.io.Serial;
import java.io.Serializable;
import java.util.Random;
import java.util.function.BiFunction;
import java.util.function.Function;

public class Matrix2D implements Serializable {
    @Serial
    private static final long serialVersionUID = 1L;
    
    private int rowsNumber;
    private int columnsNumber;
    private double[][] matrixArray;
    
    public void setColumn(int columnIndex, double[] columnValues) {
        if (columnValues.length != rowsNumber) {
            throw new IllegalArgumentException("Wrong length of the values array. Should be equal to the number of rows");
        }
        for (int row = 0; row < rowsNumber; row++) {
            matrixArray[row][columnIndex] = columnValues[row];
        }
    }
    
    public static class Element {
        public int row;
        public int column;
        public double value;
    }
    
    public Matrix2D(int rowsNumber, int columnsNumber) {
        this.rowsNumber = rowsNumber;
        this.columnsNumber = columnsNumber;
        this.matrixArray = new double[rowsNumber][columnsNumber];
    }
    
    public Matrix2D(double[][] matrixArray) {
        setMatrixArray(matrixArray);
    }
    
    public int getColumnsNumber() {
        return columnsNumber;
    }
    
    public int getRowsNumber() {
        return rowsNumber;
    }
    
    public void setMatrixArray(double[][] matrixArray) {
        if (matrixArray == null || matrixArray.length == 0) {
            this.rowsNumber = 0;
            this.columnsNumber = 0;
            this.matrixArray = new double[0][0];
        } else {
            this.rowsNumber = matrixArray.length;
            this.columnsNumber = matrixArray[0].length;
            this.matrixArray = matrixArray;
        }
    }
    
    public double[][] getMatrixArray() {
        return matrixArray;
    }
    
    public void setValue(int row, int column, double value) {
        matrixArray[row][column] = value;
    }
    
    public double getValue(int row, int column) {
        return matrixArray[row][column];
    }
    
    public Element getMaxValueElement() {
        Element max = new Element();
        max.value = Double.NEGATIVE_INFINITY;
        for (int r = 0; r < rowsNumber; r++) {
            for (int c = 0; c < columnsNumber; c++) {
                double value = getValue(r, c);
                if (value > max.value) {
                    max.value = value;
                    max.row = r;
                    max.column = c;
                }
            }
        }
        return max;
    }
    
    public Element getMaxMagnitudeElement() {
        Element max = new Element();
        max.value = 0;
        for (int r = 0; r < rowsNumber; r++) {
            for (int c = 0; c < columnsNumber; c++) {
                double magnitude = Math.abs(getValue(r, c));
                if (magnitude > max.value) {
                    max.value = magnitude;
                    max.row = r;
                    max.column = c;
                }
            }
        }
        return max;
    }
    
    public double getColumnSum(int column) {
        double result = 0;
        for (int row = 0; row < rowsNumber; row++) {
            result += matrixArray[row][column];
        }
        return result;
    }
    
    public double getRowSum(int row) {
        double result = 0;
        for (int column = 0; column < columnsNumber; column++) {
            result += matrixArray[row][column];
        }
        return result;
    }
    
    public void fillWithGaussianRandoms(double[] matrix) {
        Random rnd = new Random();
        for (int c = 0; c < columnsNumber; c++) {
            matrix[c] = rnd.nextGaussian();
        }
    }
    
    public void fillWithGaussianRandoms(double[][] matrix) {
        Random rnd = new Random();
        for (int r = 0; r < rowsNumber; r++) {
            for (int c = 0; c < columnsNumber; c++) {
                matrix[r][c] = rnd.nextGaussian();
            }
        }
    }
    
    public void fillWithGaussianRandoms() {
        fillWithGaussianRandoms(this.matrixArray);
    }
    
    public void fillWithUniformRandoms(double[] matrix) {
        Random rnd = new Random();
        for (int c = 0; c < columnsNumber; c++) {
            matrix[c] = rnd.nextDouble();
        }
    }
    
    public void fillWithUniformRandoms(double[][] matrix) {
        Random rnd = new Random();
        for (int r = 0; r < rowsNumber; r++) {
            for (int c = 0; c < columnsNumber; c++) {
                matrix[r][c] = rnd.nextDouble();
            }
        }
    }
    
    public void fillWithUniformRandoms() {
        fillWithUniformRandoms(matrixArray);
    }
    
    public void fillWithZeros() {
        matrixArray = new double[rowsNumber][columnsNumber];
    }
    
    public void fill(double value) {
        for (int r = 0; r < rowsNumber; r++) {
            for (int c = 0; c < columnsNumber; c++) {
                matrixArray[r][c] = value;
            }
        }
    }
    
    @SuppressWarnings("DuplicatedCode")
    public double[] flat() {
        double[] vector = new double[rowsNumber * columnsNumber];
        for (int r = 0; r < rowsNumber; r++) {
            if (columnsNumber > 0) {
                System.arraycopy(matrixArray[r], 0, vector, r * columnsNumber, columnsNumber);
            }
        }
        return vector;
    }
    
    public Matrix2D inverse() {
        return new Matrix2D(inverse(this.matrixArray));
    }
    
    public double determinant() {
        return determinant(matrixArray);
    }
    
    public double cofactor(int row, int column) {
        return cofactor(matrixArray, row, column);
    }
    
    public Matrix2D cofactors() {
        return new Matrix2D(cofactors(matrixArray));
    }
    
    public double minor(int row, int column) {
        return minor(matrixArray, row, column);
    }
    
    public Matrix2D minors() {
        return new Matrix2D(minors(matrixArray));
    }
    
    public Matrix2D verticalFlip() {
        return new Matrix2D(verticalFlip(matrixArray));
    }
    
    public Matrix2D horizontalFlip() {
        return new Matrix2D(horizontalFlip(matrixArray));
    }
    
    public Matrix2D antiTranspose() {
        return new Matrix2D(antiTranspose(matrixArray));
    }
    
    public Matrix2D transpose() {
        return new Matrix2D(transpose(matrixArray));
    }
    
    public Matrix2D multiply(Matrix2D m) {
        return new Matrix2D(multiply(matrixArray, m.getMatrixArray()));
    }
    
    public Matrix2D scale(double c) {
        return new Matrix2D(scale(c, matrixArray));
    }
    
    public Matrix2D add(Matrix2D m) {
        return new Matrix2D(add(matrixArray, m.getMatrixArray()));
    }
    
    public Matrix2D subtract(Matrix2D m) {
        return new Matrix2D(subtract(matrixArray, m.getMatrixArray()));
    }
    
    public Matrix2D dotProduct(Matrix2D m) {
        return new Matrix2D(dotProduct(matrixArray, m.getMatrixArray()));
    }
    
    public Matrix2D getRowMatrix(int rowIndex) {
        return new Matrix2D(getRowMatrixArray(rowIndex));
    }
    
    public Matrix2D getColumnMatrix(int columnIndex) {
        return new Matrix2D(getColumnMatrixArray(columnIndex));
    }
    
    public double[] getRowVectorArray(int rowIndex) {
        return matrixArray[rowIndex];
    }
    
    public double[][] getRowMatrixArray(int rowIndex) {
        return new double[][]{matrixArray[rowIndex]};
    }
    
    public double[] getColumnVectorArray(int columnIndex) {
        double[] column = new double[rowsNumber];
        for (int row = 0; row < rowsNumber; row++) {
            column[row] = matrixArray[row][columnIndex];
        }
        return column;
    }
    
    public double[][] getColumnMatrixArray(int columnIndex) {
        double[][] column = new double[rowsNumber][1];
        for (int row = 0; row < rowsNumber; row++) {
            column[row][0] = matrixArray[row][columnIndex];
        }
        return column;
    }
    
    public static double[][] inverse(double[][] m) {
        return scale(1 / determinant(m), transpose(cofactors(m)));
    }
    
    public static double determinant(double[][] m) {
        if (m == null || m.length == 0 ||
            m.length != m[0].length) {
            throw new IllegalArgumentException("The argument of the getDeterminant method should be a nonempty square matrix");
        }
        if (m.length == 1) {
            return m[0][0];
        }
        if (m.length == 2) {
            return m[0][0] * m[1][1] - m[0][1] * m[1][0];
        }
        double result = 0;
        for (int r = 0; r < m.length; r++) {
            result += m[r][0] * cofactor(m, r, 0);
        }
        return result;
    }
    
    public static double cofactor(double[][] m, int row, int col) {
        // row + col + 2 -> 2 to get only positive reminders
        return (row + col + 2) % 2 == 0 ? minor(m, row, col) : -minor(m, row, col);
    }
    
    public static double minor(double[][] m, int row, int col) {
        double[][] minorMatrix = new double[m.length - 1][m.length - 1];
        for (int r = 0; r < m.length; r++) {
            int minorRow;
            if (r < row) {
                minorRow = r;
            } else if (r > row) {
                minorRow = r - 1;
            } else {
                continue;
            }
            for (int c = 0; c < m.length; c++) {
                int minorCol;
                if (c < col) {
                    minorCol = c;
                } else if (c > col) {
                    minorCol = c - 1;
                } else {
                    continue;
                }
                minorMatrix[minorRow][minorCol] = m[r][c];
            }
        }
        return determinant(minorMatrix);
    }
    
    public static double[][] verticalFlip(double[][] m) {
        if (m == null || m.length == 0) {
            return null;
        }
        double[][] result = new double[m.length][m[0].length];
        double middle = (double) m.length / 2;
        for (int r = 0; r < middle; r++) {
            for (int c = 0; c < m[0].length; c++) {
                result[r][c] = m[m.length - r - 1][c];
                result[m.length - r - 1][c] = m[r][c];
            }
        }
        return result;
    }
    
    public static double[][] horizontalFlip(double[][] m) {
        if (m == null || m.length == 0) {
            return null;
        }
        double[][] result = new double[m.length][m[0].length];
        double middle = (double) m[0].length / 2;
        for (int r = 0; r < m.length; r++) {
            for (int c = 0; c < middle; c++) {
                result[r][c] = m[r][m[0].length - c - 1];
                result[r][m[0].length - c - 1] = m[r][c];
            }
        }
        return result;
    }
    
    public static double[][] antiTranspose(double[][] m) {
        if (m == null || m.length == 0) {
            return null;
        }
        double[][] result = new double[m[0].length][m.length];
        for (int r = 0; r < result.length; r++) { // result rows
            for (int c = 0; c < result[0].length; c++) { // result columns
                result[r][c] = m[m.length - c - 1][m[0].length - r - 1];
            }
        }
        return result;
    }
    
    public static double[][] transpose(double[][] m) {
        if (m == null || m.length == 0) {
            return null;
        }
        double[][] result = new double[m[0].length][m.length];
        for (int r = 0; r < result.length; r++) { // result rows
            for (int c = 0; c < result[0].length; c++) { // result columns
                result[r][c] = m[c][r];
            }
        }
        return result;
    }
    
    // Scalar matrix operations
    public static double[][] generalScalarMatrixOperator(double scalar, double[][] matrix,
                                                         BiFunction<Double, Double, Double> operator) {
        return generalScalarMatrixOperator(scalar, matrix, operator, true);
    }
    
    public static double[][] generalScalarMatrixOperator(double scalar, double[][] matrix,
                                                         BiFunction<Double, Double, Double> operator,
                                                         boolean throwError) {
        
        if (matrix == null || matrix.length == 0) {
            
            if (throwError) {
                throw new IllegalArgumentException("Scalar matrix operation error: wrong matrices' dimensions or null pointer exception");
            }
            
            return null;
        }
        
        double[][] result = new double[matrix.length][matrix[0].length];
        
        for (int r = 0; r < matrix.length; r++) {
            for (int c = 0; c < matrix[0].length; c++) {
                result[r][c] = operator.apply(scalar, matrix[r][c]);
            }
        }
        
        return result;
        
    }
    
    public static double[][] scale(double multiplier, double[][] matrix) {
        return generalScalarMatrixOperator(multiplier, matrix, (c, m) -> c * m);
    }
    
    public static double[][] divide(double divider, double[][] matrix) {
        return generalScalarMatrixOperator(divider, matrix, (d, m) -> m / d);
    }
    
    // Pairwise matrix operations
    public static double[][] generalPairwiseMatrixOperator(double[][] matrix1, double[][] matrix2,
                                                           BiFunction<Double, Double, Double> operator) {
        return generalPairwiseMatrixOperator(matrix1, matrix2, operator, true);
    }
    
    public static double[] generalPairwiseMatrixOperator(double[] matrix1, double[] matrix2,
                                                         BiFunction<Double, Double, Double> operator) {
        return generalPairwiseMatrixOperator(matrix1, matrix2, operator, true);
    }
    
    public static double[][] generalPairwiseMatrixOperator(double[][] matrix1, double[][] matrix2,
                                                           BiFunction<Double, Double, Double> operator,
                                                           boolean throwError) {
        
        if (matrix1 == null || matrix1.length == 0
            || matrix1.length != matrix2.length
            || matrix1[0].length != matrix2[0].length) {
            
            if (throwError) {
                throw new IllegalArgumentException("Pairwise matrix operation error: wrong matrices' dimensions or null pointer exception");
            }
            
            return null;
        }
        
        double[][] result = new double[matrix1.length][matrix1[0].length];
        
        for (int r = 0; r < matrix1.length; r++) {
            for (int c = 0; c < matrix1[0].length; c++) {
                result[r][c] = operator.apply(matrix1[r][c], matrix2[r][c]);
            }
        }
        
        return result;
        
    }
    
    public static double[] generalPairwiseMatrixOperator(double[] matrix1, double[] matrix2,
                                                         BiFunction<Double, Double, Double> operator,
                                                         boolean throwError) {
        
        if (matrix1 == null || matrix1.length == 0
            || matrix1.length != matrix2.length) {
            
            if (throwError) {
                throw new IllegalArgumentException("Pairwise matrix operation error: wrong matrices' dimensions or null pointer exception");
            }
            
            return null;
        }
        
        double[] result = new double[matrix1.length];
        
        for (int c = 0; c < matrix1.length; c++) {
            result[c] = operator.apply(matrix1[c], matrix2[c]);
        }
        
        return result;
        
    }
    
    public static double[][] add(double[][] m1, double[][] m2) {
        return generalPairwiseMatrixOperator(m1, m2, Double::sum);
    }
    
    public static double[][] subtract(double[][] m1, double[][] m2) {
        return generalPairwiseMatrixOperator(m1, m2, (a, b) -> a - b);
    }
    
    public static double[] subtract(double[] m1, double[] m2) {
        return generalPairwiseMatrixOperator(m1, m2, (a, b) -> a - b);
    }
    
    public static double[][] dotProduct(double[][] m1, double[][] m2) {
        return generalPairwiseMatrixOperator(m1, m2, (a, b) -> a * b);
    }
    
    public static double[] dotProduct(double[] m1, double[] m2) {
        return generalPairwiseMatrixOperator(m1, m2, (a, b) -> a * b);
    }
    
    public static double[][] generalCrossMatrixOperator(double[][] matrix1, double[][] matrix2,
                                                        BiFunction<Double, Double, Double> operator) {
            return generalCrossMatrixOperator(matrix1, matrix2, operator, true);
    }
    
    public static double[][] generalCrossMatrixOperator(double[][] matrix1, double[][] matrix2,
                                                        BiFunction<Double, Double, Double> operator,
                                                        boolean throwError) {
        
        if (matrix1 == null || matrix2 == null
            || matrix1.length == 0 || matrix2.length == 0
            || matrix1[0].length != matrix2.length) {
            
            if (throwError) {
                throw new IllegalArgumentException("Cross matrix operation error: wrong matrices' dimensions or null pointer exception");
            }
            
            return null;
            
        }
        
        double[][] result = new double[matrix1.length][matrix2[0].length];
        
        for (int r = 0; r < matrix1.length; r++) {
            for (int c = 0; c < matrix2[0].length; c++) {
                for (int i = 0; i < matrix2.length; i++) {
                    result[r][c] += operator.apply(matrix1[r][i], matrix2[i][c]);
                }
            }
        }
        
        return result;
        
    }
    
    public static double[][] multiply(double[][] m1, double[][] m2) {
        return generalCrossMatrixOperator(m1, m2, (a, b) -> a * b);
    }
    
    // Matrix modifications
    public static Matrix2D generalMatrixCellsModification(Matrix2D matrix,
                                                          Function<Double, Double> operation) {
        return new Matrix2D(generalMatrixCellsModification(matrix.getMatrixArray(), operation, true));
    }
    
    public static double[][] generalMatrixCellsModification(double[][] matrix,
                                                            Function<Double, Double> operation) {
        return generalMatrixCellsModification(matrix, operation, true);
    }
    
    public static Matrix2D generalMatrixCellsModification(Matrix2D matrix,
                                                          Function<Double, Double> operation,
                                                          boolean throwError) {
        return new Matrix2D(generalMatrixCellsModification(matrix.getMatrixArray(), operation, true));
    }
    
    public static double[][] generalMatrixCellsModification(double[][] matrix,
                                                            Function<Double, Double> operation,
                                                            boolean throwError) {
        if (matrix == null) {
            if (throwError) {
                throw new IllegalArgumentException("Matrix modification operation error: matrix can't be Null");
            }
            return null;
        }
        
        if (matrix.length == 0) {
            return matrix;
        }
        
        double[][] result = new double[matrix.length][matrix[0].length];
        
        for (int r = 0; r < matrix.length; r++) {
            for (int c = 0; c < matrix[0].length; c++) {
                result[r][c] = operation.apply(matrix[r][c]);
            }
        }
        
        return result;
        
    }
    
    // Calculations of submatrices
    public static double[][] generalSubmatricesCalculations(double[][] matrix, // I'm not sure about name for this method
                                                            Function<double[][], BiFunction<Integer, Integer, Double>> operation,
                                                            boolean throwError) throws IllegalArgumentException{
        if (matrix == null) {
            if (throwError) {
                throw new IllegalArgumentException("Submatrix calculations error: matrix can't be Null");
            }
            return null;
        }
        
        if (matrix.length == 0) {
            return matrix;
        }
        
        double[][] result = new double[matrix.length][matrix[0].length];
        
        for (int r = 0; r < matrix.length; r++) {
            for (int c = 0; c < matrix[0].length; c++) {
                result[r][c] = operation.apply(matrix).apply(r, c);
            }
        }
        
        return result;
        
    }
    
    public static double[][] cofactors(double[][] matrix) {
        double[][] result = new double[matrix.length][matrix[0].length];
        for (int r = 0; r < matrix.length; r++) {
            for (int c = 0; c < matrix[0].length; c++) {
                result[r][c] = cofactor(matrix, r, c);
            }
        }
        return result;
    }
    
    public static double[][] minors(double[][] matrix) {
        double[][] result = new double[matrix.length][matrix[0].length];
        for (int r = 0; r < matrix.length; r++) {
            for (int c = 0; c < matrix[0].length; c++) {
                result[r][c] = minor(matrix, r, c);
            }
        }
        return result;
    }
    
    @Override
    public String toString() {
        StringBuilder result = new StringBuilder();
        
        if (matrixArray == null || matrixArray.length == 0) {
            return "";
        }
        
        for (int r = 0; r < rowsNumber; r++) {
            result.append("[");
            for (int c = 0; c < columnsNumber; c++) {
                result.append(String.format("%+3.2f ", matrixArray[r][c])).append(" ");
            }
            result = new StringBuilder(result.toString().trim());
            result.append("]\n");
        }
        
        return result.toString();
    }
}