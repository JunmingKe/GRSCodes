# GRSCodes

This is the implementation of Generalized Reed-Solomon (GRS) Codes by Java.

Please contact junmingke1994@gmail.com if you have any question.

## How to use

> // decoding
> // The first line contains one prime number  which defines the field .
> // The second line contains two integers  and , parameters of the GRS code.
> // The third line contains  numbers separated by the spaces: . Each number is a non-zero element of  represented as an integer in the range .
> // The fourth line contains  numbers separated by the spaces: . Each number is a non-zero element of  represented as an integer in the range .
> // The fifth line contains  numbers separated by the spaces: . Each number is an element of  represented as an integer in the range .

```
        int p = 7;
        int n = 6, k = 2;
        int d = n - k + 1;
        int[] alpha = new int[]{1,2,3,4,5,6};
        int[] v = new int[]{1,1,1,1,1,1};
        int[] y = new int[]{1,2,5,4,5,1};
        // expect output {1,2,3,4,5,6}
        int[][] H = generateGRS(alpha, v, n, k, p);
        int[] Sx = multiplyHY(H, y, p);
        int[][] RTx = euclidPoly(Arrays.copyOfRange(Sx, 0, n - k), d, p);
        int[][] errorLocator = polyRoot(RTx[1], p);
        int[] derivativeTx = derivate(RTx[1], p);
        int[] errorEvaluator = errorEvaluate(errorLocator, derivativeTx, RTx[0], p);
        int[] codeword = polyMinus(p, y.length - 1, y, errorEvaluator.length - 1, errorEvaluator);
```
