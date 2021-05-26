package kjm.java;

import java.util.Arrays;
import java.util.HashMap;

public class GRScodes {


    public static void main(String[] args) {

        // polynomial plus
//        int primeNum = 3;
//        int degreeA = 4;
//        int[] polyA = {2,1,2,1,2};
//        int degreeB = 5;
//        int[] polyB = {1,0,2,0,0,1};
//        int[] ans = polyPlus(primeNum, degreeA, polyA, degreeB, polyB);
//        HashMap<Integer, String> manageAns = manage(ans);
//        System.out.println(manageAns.get(1));
//        System.out.println(manageAns.get(2));

        // polynomial minus
//        ans = polyMinus(primeNum, degreeA, polyA, degreeB, polyB);
//        manageAns = manage(ans);
//        System.out.println(manageAns.get(1));
//        System.out.println(manageAns.get(2));

        // polynomial multiplication
        // parameter
//        primeNum = 2;
//        degreeA = 4;
//        polyA = new int[]{1,0,1,0,1};
//        degreeB = 2;
//        polyB = new int[]{1,0,0};
//        // calculation
//        ans = polyMulti(primeNum, degreeA, polyA, degreeB, polyB);
//        manageAns = manage(ans);
//        System.out.println(primeNum + " " + manageAns.get(1));
//        System.out.println(manageAns.get(2));

        // polynomial remainder
        // parameter
//        primeNum = 7;
//        degreeA = 6;
//        polyA = new int[]{3,6,5,1,0,1,3};
//        degreeB = 6;
//        polyB = new int[]{3,6,5,1,2,1,3};
//        // calculation
//        int[][] ansRemain = polyRemainder(primeNum, degreeA, polyA, degreeB, polyB);
//        HashMap<Integer, String> manage2Ans = manage2(ansRemain, primeNum);
//        System.out.println(manage2Ans.get(0));
//        System.out.println(manage2Ans.get(1));
//        System.out.println(manage2Ans.get(2));

        // gcd
        // parameter
//        primeNum = 7;
//        degreeA = 6;
//        polyA = new int[]{3,6,5,1,0,1,3};
//        degreeB = 3;
//        polyB = new int[]{2,6,5,5};
//        // calculation
//        ansRemain = polyGCD(primeNum, degreeA, polyA, degreeB, polyB);
//        String ans1 = arrayToString(ansRemain[0], "para");
//        String ans2 = arrayToString(ansRemain[1], "poly");
//        System.out.println(ans1);
//        System.out.println(ans2);
//        System.out.println("==============================================");

        // decoding
        // The first line contains one prime number  which defines the field .
        // The second line contains two integers  and , parameters of the GRS code.
        // The third line contains  numbers separated by the spaces: . Each number is a non-zero element of  represented as an integer in the range .
        // The fourth line contains  numbers separated by the spaces: . Each number is a non-zero element of  represented as an integer in the range .
        // The fifth line contains  numbers separated by the spaces: . Each number is an element of  represented as an integer in the range .
        int p = 7;
        int n = 6, k = 2;
        int d = n - k + 1;
        int[] alpha = new int[]{1,2,3,4,5,6};
        int[] v = new int[]{1,1,1,1,1,1};
        int[] y = new int[]{4,1,0,2,6,3};
        // expect output {1,2,3,4,5,6}
        int[][] H = generateGRS(alpha, v, n, k, p);
        int[] Sx = multiplyHY(H, y, p);
        int[][] RTx = euclidPoly(Arrays.copyOfRange(Sx, 0, n - k), d, p);
        int[][] errorLocator = polyRoot(alpha, RTx[1], p);
        int[] derivativeTx = derivate(RTx[1], p);
        int[] errorEvaluator = errorEvaluate(v, alpha, errorLocator, derivativeTx, RTx[0], p);
        errorEvaluator = Arrays.copyOfRange(errorEvaluator, 0, y.length);
        int[] codeword = polyMinus(p, y.length - 1, y, errorEvaluator.length - 1, errorEvaluator);
        System.out.println(manage(codeword).get(2));

        p = 11;
        n = 10;
        k = 5;
        d = n - k + 1;
        alpha = new int[]{4, 9, 10, 7, 6, 3, 8, 2, 5, 1};
        v = new int[]{6, 10, 7, 8, 2, 4, 9, 5, 1, 3};
        y = new int[]{5, 7, 8, 5, 5, 0, 0, 0, 0, 1};
        // expect output {5 7 2 5 4 0 0 0 0 1}
        H = generateGRS(alpha, v, n, k, p);
        Sx = multiplyHY(H, y, p);
        RTx = euclidPoly(Arrays.copyOfRange(Sx, 0, n - k), d, p);
        errorLocator = polyRoot(alpha, RTx[1], p);
        derivativeTx = derivate(RTx[1], p);
        errorEvaluator = errorEvaluate(v, alpha, errorLocator, derivativeTx, RTx[0], p);
        errorEvaluator = Arrays.copyOfRange(errorEvaluator, 0, y.length);
        codeword = polyMinus(p, y.length - 1, y, errorEvaluator.length - 1, errorEvaluator);
        System.out.println(Arrays.toString(codeword));

        p = 29;
        n = 20;
        k = 7;
        d = n - k + 1;
        alpha = new int[]{19, 5, 4, 6, 10, 3, 16, 27, 2, 15, 26, 28, 23, 14, 20, 25, 22, 18, 24, 13};
        v = new int[]{9, 25, 21, 3, 8, 18, 17, 23, 4, 28, 20, 2, 27, 15, 22, 11, 14, 13, 24, 16};
        y = new int[]{19, 10, 17, 17, 9, 20, 1, 8, 22, 18, 14, 21, 12, 17, 15, 25, 7, 22, 6, 5};
        // expect output {19 18 9 9 9 20 1 8 22 18 14 21 12 17 15 25 7 22 6 5}
        H = generateGRS(alpha, v, n, k, p);
        Sx = multiplyHY(H, y, p);
        RTx = euclidPoly(Arrays.copyOfRange(Sx, 0, n - k), d, p);
        errorLocator = polyRoot(alpha, RTx[1], p);
        derivativeTx = derivate(RTx[1], p);
        errorEvaluator = errorEvaluate(v, alpha, errorLocator, derivativeTx, RTx[0], p);
        errorEvaluator = Arrays.copyOfRange(errorEvaluator, 0, y.length);
        codeword = polyMinus(p, y.length - 1, y, errorEvaluator.length - 1, errorEvaluator);
        System.out.println(Arrays.toString(codeword));

    }

    private static int[] derivate(int[] rTx, int p) {
        int[] derivativeTx = new int[rTx.length - 1];
        for (int i = 0; i < rTx.length - 1; i++) {
            derivativeTx[i] = (rTx[i + 1] * (i + 1)) % p;
        }
        return derivativeTx;
    }

    private static int[] errorEvaluate(int[] v, int[] alpha, int[][] errorLocator, int[] derivativeTx, int[] rTx, int p) {
        int[] errorEvaluator = new int[errorLocator[1].length];
        for (int i = 0; i < errorLocator[1].length; i++) {
            if (errorLocator[1][i] != 0) {
                int para = -alpha[errorLocator[0][i] - 1];
                int gamma = calculateFx(i, rTx, p);
                int lambda = calculateFx(i, derivativeTx, p);
                int errorEva = para * calculateDivision(gamma, lambda, p);
                if (errorEva < 0) {
                    errorEva = (errorEva % p) + p;
                }
                for (int j = 0; j < p; j++) {
                    int evaValue = ((j * v[errorLocator[0][i] - 1]) % p);
                    if (evaValue == errorEva) {
                        errorEvaluator[errorLocator[0][i]] = j;
                    }
                }
            }
        }
        errorEvaluator = turnToPositive(errorEvaluator, p);
        return Arrays.copyOfRange(errorEvaluator, 1, errorLocator[1].length);
    }

    private static int calculateDivision(int gamma, int lambda, int p) {
        int divisionResult = 0;
        for (int i = 0; i < p; i++) {
            if (((i * lambda) % p) == (gamma % p)) {
                divisionResult = i;
            }
        }
        return divisionResult;
    }

    private static int calculateFx(int i, int[] rTx, int p) {
        int valueFx = 0;
        for (int j = 0; j < rTx.length; j++) {
            valueFx = valueFx + (int) (rTx[j] * (Math.pow(i,j)) % p) % p;
        }
        return valueFx % p;
    }

    private static double mathPow(int i, int j, int p) {
        double mathPower = 1;
        if (j == 0) {
            if (i == 0) {
                mathPower = 0;
            } else {
                mathPower = 1;
            }
        } else {
            for (double ii = 1; ii < j + 1; ii++) {
                mathPower = (mathPower * i) % p;
            }
        }
        return mathPower;
    }

    private static int[][] polyRoot(int[] alpha, int[] rTx, int p) {
        int[][] rootLocator = new int[2][p];
        int[] rootPoly = new int[p];
        int[] polyResult = new int[p];
        for (int i = 1; i < p; i++) {
            double sumPoly = 0;
            for (int j = 0; j < rTx.length; j++) {
                double mathPower = mathPow(i, j, p);
                double rTxMath = (rTx[j] * mathPower) % ((double) p);
                sumPoly = (sumPoly + rTxMath) % ((double) p);
            }
            polyResult[i] = (int) sumPoly;
            if ((sumPoly % ((double) p)) == 0) {
                rootPoly[i] = 1;
            }
        }
        int[] errorLocator = new int[p];
        for (int i = 0; i < p; i++) {
            if (rootPoly[i] != 0) {
                for (int j = 0; j < p; j++) {
                    if ((i * j) % p == 1) {
                        errorLocator[i] = j;
                        break;
                    }
                }
            }
        }
        for (int i = 0; i < errorLocator.length; i++) {
            if (errorLocator[i] != 0) {
                for (int j = 0; j < alpha.length; j++) {
                    if (alpha[j] == errorLocator[i]) {
                        rootPoly[i] = j + 1;
                        break;
                    }
                }
            }
        }
        rootLocator[0] = rootPoly;
        rootLocator[1] = errorLocator;
        return rootLocator;
    }

    private static int[][] euclidPoly(int[] sx, int d, int p) {
        int[][] rtx = new int[2][sx.length];
        // r_{-1}
        int[] rPre = new int[d];
        rPre[d - 1] = 1;
        // r_0
        int[] rCur = sx;
        rCur = polyManage2(rCur);
        // t_{-1}
        int[] tPre = new int[d];
        tPre = polyManage2(tPre);
        // t_{0}
        int[] tCur = new int[d];
        tCur[0] = 1;
        tCur = polyManage2(tCur);
        while (rCur.length > ((d - 1) / 2))
        {
            int[][] ansRemain = polyRemainder(p, rPre.length - 1, rPre, rCur.length - 1, rCur);
            int[] qxPoly = ansRemain[0];
            qxPoly = polyManage2(qxPoly);
            rPre = rCur;
            rCur = polyManage2(ansRemain[1]);
            int[] tMid1 = polyMulti(p, qxPoly.length - 1, qxPoly, tCur.length - 1, tCur);
            tMid1 = polyManage2(tMid1);
            int[] tMid2 = polyMinus(p, tPre.length - 1, tPre, tMid1.length - 1, tMid1);
            tPre = tCur;
            tCur = polyManage2(tMid2);
        }
        rtx[0] = rCur;
        rtx[1] = tCur;
        return rtx;
    }

    private static int[] polyManage2(int[] rCur) {
        int[] ansMat;
        int degree = -1;
        for (int i = 0; i < rCur.length; i++) {
            if (rCur[i] != 0) {
                degree = i;
            }
        }
        if (degree == -1) {
            ansMat = new int[]{0};
        } else {
            ansMat = new int[degree + 1];
            for (int i = 0; i < ansMat.length; i++)
            {
                ansMat[i] = rCur[i];
            }
        }
        return ansMat;
    }

    private static int[] multiplyHY(int[][] h, int[] y, int p) {
        int[] S = new int[h[0].length];
        for (int i = 0; i < h[1].length; i++)
        {
            for (int j = 0; j < h.length; j++)
            {
                S[i] = (S[i] + h[j][i] * y[j]) % p;
            }
        }
        return S;
    }

    private static int[][] generateGRS(int[] alpha, int[] v, int n, int k, int p) {
        int[][] H = new int[alpha.length][n - k];
        for (int i = 0; i < n - k; i++)
        {
            for (int j = 0; j < alpha.length; j++)
            {
                H[j][i] = (int) ((Math.pow(alpha[j], i) % p) * v[j] % p);
            }
        }
        return H;
    }

    private static String arrayToString(int[] ansRemain, String classification) {
        String polyStr = String.valueOf(ansRemain[0]);
        if (classification.equals("para")) {
            for (int i = 1; i < 2; i++)
            {
                polyStr = polyStr + " " + ansRemain[i];
            }
        } else {
            for (int i = 1; i < ansRemain.length; i++)
            {
                polyStr = polyStr + " " + ansRemain[i];
            }
        }

        return polyStr;
    }

    private static int[][] polyGCD(int primeNum, int degreeA, int[] polyA, int degreeB, int[] polyB) {
        int[][] ans;
        if (degreeA == -1 || degreeB == -1) {
            ans = normalizeGCD(degreeA == -1 ? polyB : polyA, primeNum);
        } else {
            int[][] qx = new int[2][Math.min(degreeA, degreeB) + 1];
            ans = calculateGCD(primeNum, degreeA, polyA, degreeB, polyB, qx);
        }
        return ans;
    }

    private static int[][] calculateQx(int primeNum, int degreeA, int[] polyA, int degreeB, int[] polyB, int[][] qx) {
        int[] rx = new int[0];
        int[] qxPoly = new int[degreeA + 1];
        while (degreeA >= degreeB) {
            int[] midPoly = new int[degreeA + 1];
            for (int i = 0; i <= degreeB; i++)
            {
                midPoly[i + degreeA - degreeB] = polyB[i];
            }
            int midDegree = degreeA;
            qxPoly[degreeA - degreeB] = qxPoly[degreeA - degreeB] + 1;
            int[] ansMinusMat = polyMinus(primeNum, degreeA, polyA, midDegree, midPoly);
            ansMinusMat = polyManage2(ansMinusMat);
            int degreeMinus = ansMinusMat.length - 1;
            if (degreeMinus < degreeB) {
                rx = ansMinusMat;
            }
            degreeA = degreeMinus;
            polyA = ansMinusMat;
        }
        qx[1] = rx;
        qx[0] = qxPoly;
        return qx;
    }

    private static int[][] calculateGCD(int primeNum, int degreeA, int[] polyA, int degreeB, int[] polyB, int[][] qx) {
        int[][] ans = new int[2][];
        int[] polyAStore = polyA;
        int[] polyBStore = polyB;
        int degreeAStore = degreeA;
        int degreeBStore = degreeB;
        while (Math.min(degreeB,degreeA) > 0) {
            int[] midPoly = Math.min(degreeA, degreeB) == degreeA ? (int[])Arrays.copyOf(polyA, degreeB + 1) : (int[])Arrays.copyOf(polyB, degreeA + 1);
            int midDegree = Math.max(degreeA, degreeB) == degreeA ? degreeA : degreeB;
            int[] ansMinusMat = polyMinus(primeNum, Math.min(degreeA, degreeB) == degreeA ? degreeB: degreeA, Math.min(degreeA, degreeB) == degreeA ? polyB: polyA, midDegree, midPoly);
            ansMinusMat = polyManage(ansMinusMat);
            int degreeMinus = ansMinusMat.length - 1;
            if (degreeMinus < midDegree) {
                midPoly = Math.min(degreeAStore, degreeBStore) == degreeAStore ? polyAStore : polyBStore;
                midDegree = Math.min(degreeAStore, degreeBStore) == degreeAStore ? degreeAStore : degreeBStore;
                polyAStore = midPoly;
                degreeAStore = midDegree;
                polyBStore = ansMinusMat;
                degreeBStore = degreeMinus;
            }
            degreeA = midDegree;
            polyA = midPoly;
            degreeB = degreeMinus;
            polyB = ansMinusMat;
        }
        int[] polyAns;
        if (Math.min(degreeB,degreeA) == 0 && polyB[0] > 0) {
            polyAns = new int[]{1};
        } else {
            polyAns = polyManage(Math.min(degreeA, degreeB) == degreeA ? polyB : polyA);
        }
        ans = normalizeGCD(polyAns, primeNum);
        return ans;
    }

    private static int[] polyManage(int[] ansMinusMat) {
        int[] ansMat;
        int degree = -1;
        for (int i = 0; i < ansMinusMat.length; i++) {
            if (ansMinusMat[i] != 0) {
                degree = i;
                break;
            }
        }
        if (degree == -1) {
            ansMat = new int[]{0};
        } else {
            ansMat = new int[ansMinusMat.length - degree];
            for (int i = degree; i < ansMinusMat.length; i++)
            {
                ansMat[i - degree] = ansMinusMat[i];
            }
        }
        return ansMat;
    }

    private static int[][] normalizeGCD(int[] polyGCD, int primeNum) {
        int[][] ans = new int[2][Math.max(polyGCD.length, 2)];
        if (polyGCD.length == 0 && polyGCD[0] == 0) {
            int[][] ansPoly = new int[2][1];
            ansPoly[0][0] = primeNum;
            ansPoly[0][1] = -1;
            ansPoly[1][0] = 0;
            ans = ansPoly;
        } else {
            int[] ansPoly = new int[polyGCD.length];
            ansPoly[0] = 1;
            for (int i = 1; i < polyGCD.length; i++) {
                for (int j = 0; j <primeNum; j++) {
                    if ((j * polyGCD[0] % primeNum) == polyGCD[i]) {
                        ansPoly[i] = j;
                        break;
                    }
                }
            }
            ans[0][0] = primeNum;
            ans[0][1] = polyGCD.length - 1;
            ans[1] = ansPoly;

        }
        return ans;
    }

    static HashMap<Integer, String> manage2(int[][] ansRemain, int primeNum) {
        HashMap<Integer, String> Sites = new HashMap<Integer, String>();
        int degQ, degR;
        int degree = -1;
        for (int i = ansRemain[0].length - 1; i >= 0; i--) {
            if (ansRemain[0][i] != 0) {
                degree = i;
                break;
            }
        }
        if (degree == -1) {
            Sites.put(1, "0");
            degQ = -1;
        } else {
            String polyStr = String.valueOf((ansRemain[0][degree] < 0 ? ansRemain[0][degree] + primeNum : ansRemain[0][degree]) % primeNum);
            for (int i = degree - 1; i >= 0; i--)
            {
                polyStr = polyStr + " " + ((ansRemain[0][i] < 0 ? ansRemain[0][i] + primeNum : ansRemain[0][i]) % primeNum);
            }
            Sites.put(1, polyStr);
            degQ = degree;
        }
        degree = -1;
        for (int i = ansRemain[1].length - 1; i >= 0; i--) {
            if (ansRemain[1][i] != 0) {
                degree = i;
                break;
            }
        }
        if (degree == -1) {
            Sites.put(2, "0");
            degR = -1;
        } else {
            String polyStr = String.valueOf((ansRemain[1][0] < 0 ? ansRemain[1][0] + primeNum : ansRemain[1][0]) % primeNum);
            for (int i = 1; i < ansRemain[1].length; i++)
            {
                polyStr = polyStr + " " + ((ansRemain[1][i] < 0 ? ansRemain[1][i] + primeNum : ansRemain[1][i]) % primeNum);
            }
            Sites.put(2, polyStr);
            degR = ansRemain[1].length - 1;
        }
        Sites.put(0, primeNum + " " + degQ + " " + degR);
        return Sites;
    }

    private static int[][] polyRemainder(int primeNum, int degreeA, int[] polyA, int degreeB, int[] polyB) {
        // define the output
        int[][] ans;
        if (degreeA == -1) {
            ans = new int[2][1];
        } else if (degreeA < degreeB) {
            ans = new int[2][degreeA + 1];
            ans[1] = polyA;
        } else {
            int[][] qx = new int[2][degreeA + 1];
            ans = calculateQx(primeNum, degreeA, polyA, degreeB, polyB, qx);
        }
        return ans;
    }
    static int[] polyPlus(int primeNum, int degreeA, int[] polyA, int degreeB, int[] polyB) {
        int[] ans;
        if (degreeA == -1 || degreeB == -1)
        {
            ans = degreeA < 0 ? polyB : polyA;
        } else {
            int minDegree = Math.min(degreeA, degreeB);
            int maxDegree = Math.max(degreeA, degreeB);
            int[] midPoly = new int[maxDegree + 1];
            for (int i = 0; i <= maxDegree; i++) {
                midPoly[i] = i < maxDegree - minDegree ? 0 : (minDegree == degreeA ? polyA[i - maxDegree + minDegree] : polyB[i - maxDegree + minDegree]);
            }
            int[] ansPoly = new int[maxDegree + 1];
            for (int i = 0; i <= maxDegree; i++) {
                ansPoly[i] = (minDegree == degreeA) ? (polyB[i] + midPoly[i]) % primeNum : (polyA[i] + midPoly[i]) % primeNum;
            }
            ans = turnToPositive(ansPoly, primeNum);
        }
        return ans;
    }

    static int[] polyMinus(int primeNum, int degreeA, int[] polyA, int degreeB, int[] polyB) {
        // define the output
        int[] ans;
        // judge the abnormal input
        if (degreeA == -1) {
            ans = turnToPositive(polyB, primeNum);
        } else if (degreeB == -1) {
            ans = polyA;
        } else {
            int minDegree = Math.min(degreeA, degreeB);
            int maxDegree = Math.max(degreeA, degreeB);
            int[] ansPoly = new int[maxDegree + 1];
            for (int i = 0; i <= minDegree; i++) {
                ansPoly[i] = (polyA[i] - polyB[i]) % primeNum;
            }
            for (int i = minDegree + 1; i <= maxDegree; i++) {
                ansPoly[i] = (minDegree == degreeA) ? (-polyB[i] % primeNum) : (polyA[i] % primeNum);
            }
            ans = turnToPositive(ansPoly, primeNum);
        }
        return ans;
    }

    static int[] forcePolyMinus(int primeNum, int degreeA, int[] polyA, int degreeB, int[] polyB) {
        // define the output
        int[] ans;
        // judge the abnormal input
        if (degreeA == -1) {
            ans = turnToPositive(polyB, primeNum);
        } else if (degreeB == -1) {
            ans = polyA;
        } else {
            int minDegree = Math.min(degreeA, degreeB);
            int maxDegree = Math.max(degreeA, degreeB);
            int[] ansPoly = new int[maxDegree + 1];
            for (int i = 0; i <= minDegree; i++) {
                ansPoly[i] = (polyA[i] - polyB[i]) % primeNum;
            }
            for (int i = minDegree + 1; i <= maxDegree; i++) {
                ansPoly[i] = (minDegree == degreeA) ? (-polyB[i] % primeNum) : (polyA[i] % primeNum);
            }
            ans = turnToPositive(ansPoly, primeNum);
        }
        return ans;
    }

    static int[] polyMulti(int primeNum, int degreeA, int[] polyA, int degreeB, int[] polyB) {
        // define the output
        int[] ans;
        // judge the abnormal input
        if (degreeA == -1 || degreeB == -1) {
            ans = new int[]{0};
        } else {
            int[] ansPoly = new int[degreeA + degreeB + 1];
            for (int i = 0; i < polyA.length; i++) {
                for (int j = 0; j < polyB.length; j++) {
                    ansPoly[i + j] = (ansPoly[i+j] + polyA[i] * polyB[j]) % primeNum;
                }
            }
            ans = ansPoly;
        }
        return ans;
    }

    static int[] turnToPositive(int[] oriMatrix, int primeNum) {
        int[] positiveMatrix = oriMatrix;
        for (int i = 0; i < oriMatrix.length; i++) {
            if (oriMatrix[i] < 0) {
                int maxMinus = oriMatrix[i] % primeNum;
                positiveMatrix[i] = maxMinus + primeNum;
            }
        }
        return positiveMatrix;
    }

    static HashMap<Integer, String> manage(int[] ans) {
        HashMap<Integer, String> Sites = new HashMap<Integer, String>();;
        if (ans.length == 1 && ans[0] == 0) {
            Sites.put(1, "-1");
            Sites.put(2, "0");
        } else {
            int degree = 0;
            for (int i = 0; i < ans.length; i++) {
                if (ans[i] != 0) {
                    degree = ans.length - i;
                    break;
                }
            }
            if (degree == 0) {
                Sites.put(1, ans[degree] == 0 ? "-1" : "0");
                Sites.put(2, ans[degree] == 0 ? "0" : String.valueOf(ans[degree]));
            } else {
                int[] newAns = Arrays.copyOfRange(ans, ans.length - degree, ans.length);
                String degreeStr = String.valueOf(degree - 1);
                Sites.put(1, degreeStr);
                String polyStr = String.valueOf(newAns[0]);
                for (int i = 1; i < newAns.length; i++)
                {
                    polyStr = polyStr + " " + newAns[i];
                }
                Sites.put(2, polyStr);
            }
        }
        return Sites;
    }


}






