#include <iostream>
#include <ctime>
#include <cmath>

const int aesSbox[] = {99, 124, 119, 123, 242, 107, 111, 197, 48, 1, 103, 43, 254, 215, 171, 118, 202, 130, 201, 125, 250,
                       89, 71, 240, 173, 212, 162, 175, 156, 164, 114, 192, 183, 253, 147, 38, 54, 63, 247, 204, 52, 165,
                       229, 241, 113, 216, 49, 21, 4, 199, 35, 195, 24, 150, 5, 154, 7, 18, 128, 226, 235, 39, 178, 117,
                       9, 131, 44, 26, 27, 110, 90, 160, 82, 59, 214, 179, 41, 227, 47, 132, 83, 209, 0, 237, 32, 252,
                       177, 91, 106, 203, 190, 57, 74, 76, 88, 207, 208, 239, 170, 251, 67, 77, 51, 133, 69, 249, 2, 127,
                       80, 60, 159, 168, 81, 163, 64, 143, 146, 157, 56, 245, 188, 182, 218, 33, 16, 255, 243, 210, 205,
                       12, 19, 236, 95, 151, 68, 23, 196, 167, 126, 61, 100, 93, 25, 115, 96, 129, 79, 220, 34, 42, 144,
                       136, 70, 238, 184, 20, 222, 94, 11, 219, 224, 50, 58, 10, 73, 6, 36, 92, 194, 211, 172, 98, 145,
                       149, 228, 121, 231, 200, 55, 109, 141, 213, 78, 169, 108, 86, 244, 234, 101, 122, 174, 8, 186, 120,
                       37, 46, 28, 166, 180, 198, 232, 221, 116, 31, 75, 189, 139, 138, 112, 62, 181, 102, 72, 3, 246, 14,
                       97, 53, 87, 185, 134, 193, 29, 158, 225, 248, 152, 17, 105, 217, 142, 148, 155, 30, 135, 233, 206,
                       85, 40, 223, 140, 161, 137, 13, 191, 230, 66, 104, 65, 153, 45, 15, 176, 84, 187, 22};

int raiseToPower(int num, int pow);

void FisherYates(int *arr, int n);

int *SBoxGeneratingDec(int n, int m, int counter);

int *valueToBinary(int i, int rank);

int *binaryElementsApprox(int *arr, int size, int count);

int *SBoxApprox(int *sbox, int size, int count);

int *elemsForN(int size);

int LATMax(int *sbox, int size, int count);

int myModulusDec(int number, int mod);

int *binaryElements(int *arr, int size, int count);

int *SBoxToBooleanFunc(int *sbox, int size, int count);

int *linearCombinations(const int *arr, int size, int count);

int myModulus(int number, int mod);

void bubbleSort(int *data, int size);

int *WHTSpectrumForLinearComb(const int *arr, int size, int count);

int *toPolarTable(const int *function, int size);

int *autoCorrelation(int *func, int size, int count);

int *ACForLinearComb(const int *arr, int size, int count);

int linearRedundancy(int *sbox, int size, int count);

void buildOneRow(int *arr, int *monomials);

int rankCalculation(int arr[256][697]);

int numOfCombinations(int n, int d);

int algebraicImmunity(const int *sbox, int size, int count);

int deltaUniformity(const int *arr, int size, int count);

int *particleSwarmOptimization(int size, int count, int N, int maxIter, int mode, int *finalIter);


int main() {
    int finalIter = 0;
    int mode = 1;
    int N = 4;
    int maxIter = 5;

    int *ar = particleSwarmOptimization(256, 8, N, maxIter, mode, &finalIter);
    system("PAUSE");
    return 0;
}

//Функція зведення до ступеня

int raiseToPower(int num, int pow) {
    int res = 1;
    for (int i = 0; i < pow; ++i) {
        res *= num;
    }
    return res;
}

//Перемішування Фішера-Йейтса

void FisherYates(int *arr, int n) {
    int i, j, tmp;

    for (i = n - 1; i > 0; i--) {
        j = rand() % (i + 1);
        tmp = arr[j];
        arr[j] = arr[i];
        arr[i] = tmp;
    }
}

//Функція генерації десяткового S-Box'у

int *SBoxGeneratingDec(int n, int m, int counter) {
    int size = raiseToPower(2, n);
    int *dec = (int *) calloc(size, sizeof(int));
    srand((counter*counter) %size);
    for (int i = 0; i < size;) {
        dec[i] = rand() % size;
        int contains = 0;
        for (int j = 0; j < i; ++j) {
            if (dec[i] == dec[j]) {
                contains = 1;
                break;
            }
        }
        if (!contains) {
            i++;
        }
    }
    /*printf("Generated s-box: ");
    for (int i = 0; i < size; ++i) {
        printf("%d, ", dec[i]);
    }
    printf("\n");*/
    FisherYates(dec,size);
    return dec;
}

//Функція перетворення числа з десяткової СЧ у двійкову СЧ

int *valueToBinary(int i, int rank) {
    int *res = (int *) calloc(rank, sizeof(int));
    for (int j = 0; j < rank; ++j) {
        res[rank - 1 - j] = i >> j & 1;
    }
    return res;
}

//Функція генерації двійкових елементів у зворотньому порядку

int *binaryElementsApprox(int *arr, int size, int count) {
    int *result = (int *) calloc(size * count, sizeof(int));
    for (int i = 0; i < size; ++i) {
        int *bin = valueToBinary(arr[i], count);
        for (int j = count - 1; j >= 0; j--) {
            result[j * size + i] = bin[j];
        }
        free(bin);
    }
    return result;
}

//Функція генерації S-Box'у у зворотньому порядку

int *SBoxApprox(int *sbox, int size, int count) {
    int *result = binaryElementsApprox(sbox, size, count);
    return result;
}

//Функція генерації чисел для вхідних векторів ступеня N

int *elemsForN(int size) {
    int *result = (int *) calloc(size, sizeof(int));
    for (int i = 0; i < size; ++i) {
        result[i] = i;
    }
    return result;
}

//Функція знаходження LAT та її максимуму

int LATMax(int *sbox, int size, int count) {
    int *ar = SBoxApprox(sbox, size, count);
    int *elems = elemsForN(size);
    int *binelems = binaryElementsApprox(elems, size, count);
    int *temp = (int *) calloc(size, sizeof(int));
    int *temp2 = (int *) calloc(size, sizeof(int));
    int *coefficients = (int *) calloc(size * size, sizeof(int));
    for (int i = 0; i < size; ++i) {
        int *bin1 = valueToBinary(i, count);
        for (int k = count - 1; k >= 0; k--) {
            if (bin1[k]) {
                //printf("K===%d ",k);
                //printf("X == \n ");
                for (int l = 0; l < size; ++l) {
                    temp[l] = temp[l] ^ binelems[k * size + l];
                    //printf("%d ",temp[l]);
                }
                //printf("\n ");
            }
        }
        //printf("\n ");
        for (int j = 0; j < size; ++j) {
            int *bin2 = valueToBinary(j, count);
            for (int q = count - 1; q >= 0; q--) {
                if (bin2[q]) {
                    //printf("K===%d ",k);
                    //printf("\nY [%d]== \n ", j);
                    for (int w = 0; w < size; ++w) {
                        temp2[w] = temp2[w] ^ ar[q * size + w];
                        //printf("%d ",temp2[l]);
                    }
                }
            }
            //printf("\n ");
            int calc = 0;
            for (int r = 0; r < size; ++r) {
                temp2[r] = temp2[r] ^ temp[r];
                //printf("%d ", temp2[l]);
                if (temp2[r] == 0) {
                    ++calc;
                }
                temp2[r] = 0;
            }
            int result = 0;
            result = calc - (size / 2);
            //printf("COEFFS = %d ", result);
            coefficients[i * size + j] = result;
            free(bin2);
        }
        for (int t = 0; t < size; ++t) {
            temp[t] = 0;
        }
        //printf("\n ");
        free(bin1);
    }
    for (int n = 0; n < size; ++n) {
        for (int m = 0; m < size; ++m) {
            //printf("%d ", coefficients[n*size+m]);
        }
        //printf("\n");
    }
    int result = 0;
    for (int p = 1; p < size * size; p++) {
        if (abs(coefficients[p]) > result)
            result = abs(coefficients[p]);
    }
    free(ar);
    free(elems);
    free(binelems);
    free(temp);
    free(temp2);
    free(coefficients);
    return result;
}

//Функція приведення числа за модулем дуякого числа

int myModulusDec(int number, int mod) {
    if (number < 0) {
        while (number < 0) {
            number = number + mod;
        }
    }
    return number % mod;
}

//Функція перетворення елементів з десяткової СЧ у двійкову СЧ, для певного ступеня N

int *binaryElements(int *arr, int size, int count) {
    int *result = (int *) calloc(size * count, sizeof(int));
    for (int i = 0; i < size; ++i) {
        int *bin = valueToBinary(arr[i], count);
        for (int j = 0, k = count - 1; j < count; ++j, k--) {
            result[j * size + i] = bin[k];
        }
        free(bin);
    }
    return result;
}

//Функція перетворення S-Box'у на набір булевих функцій

int *SBoxToBooleanFunc(int *sbox, int size, int count) {
    //printf("\nS-BOX\n");
    /*for (int i = 0; i < size; ++i) {
        printf("%d ", sbox[i]);
    }*/
    //printf("\n");
    //printf("\nS-BOX IN BOOLEAN FUNCTIONS REPRESENTATION\n");
    int *result = binaryElements(sbox, size, count);
    /*for (int i = 0; i < count; ++i) {
        printf("Function %d = ", i + 1);
        for (int j = 0; j < size; ++j) {
            printf("%d ", result[i * size + j]);
        }
        printf("\n");
    }*/

    //printf("\n");

    /*for (int i = 0; i < count; ++i) {
         int *temp = calloc(size, sizeof(int));
         //printf("Function %d", i);
         for (int j = 0; j < size; ++j) {
             temp[j] = result[i * size + j];
         }
         int weight = HammingWeight(temp, size);
         int flag = funcIsBalanced(weight, count);
         //printf("\n");
         free(temp);
     }*/
    return result;
}

//Функція знаходження лінійних комбінацій для булевих функцій S-Box'у

int *linearCombinations(const int *arr, int size, int count) {
    int *result = (int *) calloc(size*(size-1), sizeof(int));
    int *calc = (int *) calloc(size, sizeof(int));
    for (int i = 1; i < size; ++i) {
        int *bin = valueToBinary(i, count);
        for (int j = 0, k = count - 1; j < count, k >= 0; ++j, k--) {
            if (bin[k] == 1) {
                for (int w = 0; w < size; ++w) {
                    calc[w] = calc[w] ^ arr[j * size + w];
                    //printf(" %d", arr[j*size]);
                    //printf(" %d", j * size + k);
                    //printf("calc =  %d", calc[k]);
                    //result[(i-1)*size+k] = calc[k];
                    //printf("result  =  %d", (i-1)*size+k);
                }
                //printf("\n");
            }
            for (int r = 0; r < size; ++r) {
                result[(i - 1) * size + r] = calc[r];
            }
            //printf(" %d", bin[j]);
        }
        for (int l = 0; l < size; ++l) {
            //printf("calc =  %d", calc[l]);
            //result[(i-1) * size + l] = calc[l];
            calc[l] = 0;
        }
        //printf("\n");
        free(bin);
    }
    free(calc);
    return result;
}

//Функція приведення числа за модулем дуякого числа

int myModulus(int number, int mod) {
    if (number < 0) {
        while (number < 0) {
            number = number + mod;
        }
    }
    return number % 2;
}

//Функція визначення коефіцієнтів перетворення Уолдша-Адамара

int *HadamardCoefficients(const int *func, int size, int count) {
    int *result = (int *) calloc(size, sizeof(int));
    int *test = (int *) calloc(size * count, sizeof(int));
    int *functions2 = elemsForN(size);
    /*for (int i = 0; i < size; ++i) {
        printf(" %d",functions2 [i]);
    }*/
    //printf("\n");
    for (int i = 0; i < size; ++i) {
        int *bin = valueToBinary(functions2[i], count);
        for (int j = 0; j < count; ++j) {
            //printf(" bin j = %d", bin[j]);
            //*(functions + i * cols + j) = (i >> cols - j - 1) & 1u;
            test[i * count + j] = bin[j];
            //printf(" %d",test [i * count + j]);
        }
        //printf("\n");
        free(bin);
    }
    int *w = (int *) calloc(count, sizeof(int));
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < count; ++j) {
            w[j] = test[i * count + j];
            //printf("w = %d", w[j]);
        }
        int res = 0;
        for (int j = 0; j < size; ++j) {
            int r = 0;
            for (int k = 0; k < count; ++k) {
                r += myModulus(w[k] * test[j * count + k], 2);
            }
            res += raiseToPower(-1, myModulus(func[j] + r, 2));
        }
        result[i] = res;
    }
    free(test);
    free(functions2);
    free(w);
    return result;
}

//Функція бульбашкового сортування

void bubbleSort(int *data, int size) {
    int i, j;
    for (i = 0; i < size; ++i) {
        for (j = size - 1; j > i; --j) {
            if (data[j] < data[j-1]) {
                int t = data[j - 1];
                data[j - 1] = data[j];
                data[j] = t;
            }
        }
    }
}

//Функція знаходження спектру Уолдша-Адамара для кожної з лінійних комбінацій

int *WHTSpectrumForLinearComb(const int *arr, int size, int count) {
    int *result = (int *) calloc(size*(size-1), sizeof(int));
    int *temp = (int *) calloc(size, sizeof(int));
    for (int i = 0; i < size - 1; ++i) {
        //printf("\nCombination %d", i + 1);
        for (int j = 0; j < size; ++j) {
            temp[j] = arr[i * size + j];
            //printf("%d ", temp[j]);
        }
        int *fxarr = HadamardCoefficients(temp, size, count);
        for (int g = 0; g< size; ++g){
            fxarr[g] = abs(fxarr[g]);
        }
        bubbleSort(fxarr,size);
        for (int g = 0; g< size; ++g){
            result[i*size+g] = fxarr[g];
        }
        //printf("\nHADAMARD COEFFICIENTS");
        //printf("\n");
        /*for (int q = 0; q < size; ++q) {
            //printf("%d ", fxarr[q]);
            result[i*size+q] = fxarr[q];
        }*/
        free(fxarr);
    }
    free(temp);
    return result;
}

//Функція представлення таблиці істиності в полярному вигляді

int *toPolarTable(const int *function, int size) {
    int *res = (int *) calloc(size, sizeof(int));
    for (int i = 0; i < size; ++i) {
        res[i] = raiseToPower(-1, function[i]);
    }
    return res;
}

//Функція обчислення автокореляційної функції

int *autoCorrelation(int *func, int size, int count) {
    int temp = 0;
    int *acFunc = (int *) calloc(size, sizeof(int));
    int *polFunc = toPolarTable(func, size);
    int *polFunc2 = (int *) calloc(size, sizeof(int));
    for (int i =0, k = size-1; i<size, k>=0; ++i, k--) {
        polFunc2[i] = polFunc[k];
        //printf("\npf2= %d", polFunc2[i]);
    }
    acFunc[0] = raiseToPower(2, count);
    for (int f = 1; f < size; ++f) {
        for (int j = 0; j < size; ++j) {
            temp = polFunc2[j] * polFunc2[j ^ f];
            //printf("\nj = %d", j);
            //printf("\nj^i = %d", j^i);
            //printf("\ntemp= %d", temp);
            acFunc[f] = acFunc[f] + temp;
        }
        //printf("ac i = %d", acFunc[i]);
    }
    free(polFunc);
    free(polFunc2);
    return acFunc;
}

//Функція знаходження автокореляційної функції для кожної з лінійних комбінацій

int *ACForLinearComb(const int *arr, int size, int count) {
    int *result = (int *) calloc(size*(size-1), sizeof(int));
    int *temp = (int *) calloc(size, sizeof(int));
    for (int i = 0; i < size - 1; ++i) {
        //printf("\nCombination %d", i + 1);
        for (int j = 0; j < size; ++j) {
            temp[j] = arr[i * size + j];
            //printf("%d ", temp[j]);
        }
        int *ar = autoCorrelation(temp, size, count);

        for (int g = 0; g< size; ++g){
            ar[g] = abs(ar[g]);
        }

        bubbleSort(ar,size);
        for (int q = 0; q < size; ++q) {
            //printf("%d ", fxarr[q]);
            result[i*size+q] = ar[q];
        }
        free(ar);
    }
    free(temp);
    return result;
}

//Функція знаходження лінійної збитковості S-Box'у

int linearRedundancy(int *sbox, int size, int count) {
    int result;
    int sp [size-1][size];
    int ac [size-1][size];
    int *ar1 = SBoxToBooleanFunc(sbox, size, count);
    //printf("\n1");
    int *ar2 = linearCombinations(ar1, size, count);
    /*for (int i = 0; i < size-1; ++i){
        for(int j = 0; j<size; ++j){
            printf("%d ", sp[i][j]);
        }
        printf("%d ",counter);
        printf("\n");
        counter++;
    }*/
    free(ar1);
    //printf("\n2");
    //int *ar5 = calloc((size-1),sizeof(int));
    //int ar33[size-1][size];
    int *ar3 = WHTSpectrumForLinearComb(ar2, size, count);
    //printf("\n3");
    /*for (int i = 0; i < size*(size-1); ++i){
        printf("%d ", ar3[i]);
    }*/
    int *ar4 = ACForLinearComb(ar2, size, count);
    free(ar2);
    //printf("\n4");
    //ar5 = DegForLinearComb(ar2,size,count);
    for (int b = 0; b < size - 1; ++b) {
        for (int q = 0; q < size; ++q) {
            sp[b][q] = ar3[b * size + q];
        }
    }
    for (int z = 0; z < size - 1; ++z) {
        for (int n = 0; n < size; ++n) {
            ac[z][n] = ar4[z * size + n];
        }
    }
    free(ar3);
    free(ar4);
    /*int dg [size-1];
    for (int i = 0; i < size - 1; ++i) {
        dg[i] = ar5[i];
    }*/
    /*FILE *file;
    fopen_s(&file, "Linear comb and WHT Spectrum.txt", "w");
    if (file == NULL) {
        printf("ERROR: Can't save sbox to file!\n");
        for (;;);
    }
    fprintf(file, "\n");
    for (int i = 0; i < size-1; ++i){
        fprintf(file,"\nLINEAR COMBINATION\n");
        for(int j = 0; j<size; ++j){
            ar33[i][j] = ar2[i*size+j];
            fprintf(file,"%d, ", ar33[i][j]);;
        }
        fprintf(file,"\n");
        fprintf(file,"\nHADAMARD SPECTRUM\n");
        for(int k = 0; k<size; ++k){
            fprintf(file,"%d ", sp[i][k]);
        }
        fprintf(file,"\n");
    }
    fprintf(file, "\n");
    fclose(file);*/

    /*printf("\nHADAMARD SPECTRUM\n");
    for (int i = 0; i < size-1; ++i){
        for(int j = 0; j<size; ++j){
            printf("%d ", sp[i][j]);
        }
        printf("\n");
    }
    printf("\nAUTO CORRELATION FUNCTIONS\n");
    for (int i = 0; i < size-1; ++i){
        for(int j = 0; j<size; ++j){
            printf("%d ", ac[i][j]);
        }
        printf("\n");
    }
    printf("\nDEGREES\n");
    for (int i = 0; i < size-1; ++i){
        printf("%d ", dg[i]);
    }
    printf("\n");*/
    int innerCounter = 0;
    int OuterCounter = 0;
    int finalCounter = 0;
    for (int i = 0; i < size - 1; ++i) {
        for (int h = i + 1; h < size - 1; ++h) {
            for (int j = 0; j < size; ++j) {
                if (i != h) {
                    if (sp[i][j] == sp[h][j] && sp[h][j] != -999 &&
                        (ac[i][j] == ac[h][j] && ac[h][j] != -999) /*(dg[i] == dg[h] && dg[h] != -999)*/) {
                        innerCounter++;
                    }
                }
            }
            //printf("\nInner counter = %d %d %d", i, h,innerCounter);
            /*printf("\n");
            printf("\n");
            for (int j = 0; j < size; ++j){
                printf("%d ", sp[i][j]);
            }
            printf("\n");
            printf("\n");
            for (int j = 0; j < size; ++j) {
                printf("%d ", sp[h][j]);
            }*/
            if (innerCounter == size) {
                //if (i!=h) {
                //if (ar4[i]==ar4[h] && ar5[i] == ar5[h]) {
                //printf("\nSTRINGS ARE EQUAL");
                //dg[h] = -999;
                OuterCounter++;
                for (int d = 0; d < size; ++d) {
                    sp[h][d] = -999;
                    //ac[h][j] = -999;
                    //dg[h] = -999;
                }
                //}
                //}
            }
            innerCounter = 0;
        }
        /*for (int j = 0; j < size; ++j){
            sp[i][j] = -999;
            ac[i][j] = -999;
            dg[i] = -999;
        }*/
        //printf("\nOC ==%d ", OuterCounter);
    }
    /*int counter = 0;
    printf("\nHADAMARD SPECTRUM AFTER\n");
    for (int i = 0; i < size-1; ++i){
        for(int j = 0; j<size; ++j){
            printf("%d ", sp[i][j]);
        }
        printf("%d ",counter);
        printf("\n");
        counter++;
    }*/
    /*printf("\nAUTO CORRELATION FUNCTIONS AFTER\n");
    for (int i = 0; i < size-1; ++i){
        for(int j = 0; j<size; ++j){
            printf("%d ", ac[i][j]);
        }
        printf("\n");
    }
    printf("\nDEGREES AFTER\n");
    for (int i = 0; i < size-1; ++i){
        printf("%d ", dg[i]);
    }*/
    finalCounter = finalCounter + OuterCounter;
    //printf("\nFINAL ==%d ", finalCounter);
    result = (size - 1) - finalCounter;
    return result;
}

//Функція побудування одного рядка матриці, що описує S-Box

void buildOneRow(int *arr, int *monomials) {
    monomials[0] = 1;
    //monomials x1,x8,y1,...,y8
    for (int i = 1; i <= 16; i++)
        monomials[i] = arr[i - 1];
    int pos = 17;
    //monomials x1x2
    for (int i = 1; i < 16; i++) {
        for (int j = i + 1; j <= 16; j++) {
            monomials[pos] = monomials[i] & monomials[j];
            pos++;
        }
    }
    //monomials x1x2x3
    for (int i = 1; i < 15; i++) {
        for (int j = i + 1; j <= 16; j++) {
            for (int k = j + 1; k <= 16; k++) {
                monomials[pos] = monomials[i] & monomials[j] & monomials[k];
                pos++;
            }
        }
    }
}

//Функція обчислення рангу матриці

int rankCalculation(int arr[256][697]) {
    int m = 697;
    int n = 256;

    int rank = 697;
    int  line_used[697] = { 0, };
    for (int i = 0;i < 697;i++)
        line_used[i] = 0;
    for (int i = 0; i < m; ++i) {
        int j;
        for (j = 0; j < n; ++j)
            if (!line_used[j] && arr[j][i])
                break;
        if (j == n)
            --rank;
        else {
            line_used[j] = 1;
            for (int k = 0; k < n; ++k)
                if (k != j && arr[k][i])
                    for (int p = i + 1; p < m; ++p)
                        arr[k][p] ^= arr[j][p] & arr[k][i];
        }
    }
    return rank;
}

//Функція знаходження кількості сполучень

int numOfCombinations(int n, int d) {
    if (n == d)
        return 1;
    if (d == 1)
        return n;
    if (d == 0)
        return 1;
    return numOfCombinations(n - 1, d - 1) + numOfCombinations(n - 1, d);
}

//Функція обчислення алгебраїчного імунітету

int algebraicImmunity(const int *sbox, int size, int count) {
    int mat[256][697];
    int tmp[697];
    int values[16];
    int *input_values = (int *) calloc(size * count, sizeof(int));
    for (int i = 0; i < size; ++i) {
        int *bin = valueToBinary(i, count);
        for (int j = 0; j < count; ++j) {
            input_values[i * count + j] = bin[j];
            //printf("%d ", input_values[i*count+j]);
        }
        //printf("\n");
        free(bin);
    }
    for (int i = 0; i < 256; i++) {
        int y = sbox[i];
        values[0] = input_values[i * count + 0];
        values[1] = input_values[i * count + 1];
        values[2] = input_values[i * count + 2];
        values[3] = input_values[i * count + 3];
        values[4] = input_values[i * count + 4];
        values[5] = input_values[i * count + 5];
        values[6] = input_values[i * count + 6];
        values[7] = input_values[i * count + 7];
        values[8] = input_values[y * count + 0];
        values[9] = input_values[y * count + 1];
        values[10] = input_values[y * count + 2];
        values[11] = input_values[y * count + 3];
        values[12] = input_values[y * count + 4];
        values[13] = input_values[y * count + 5];
        values[14] = input_values[y * count + 6];
        values[15] = input_values[y * count + 7];
        buildOneRow((int *) &values, (int *) &mat[i]);
    }
    int rank = rankCalculation(mat);
    free(input_values);
    //printf("%d", rank);
    return rank == 256 ? 3 : 2;
}

//Функція знаходження дельта-рівномірності

int deltaUniformity(const int *arr, int size, int count) {
    int result;
    int max = 0;
    for (int a = 1; a < size; ++a) {
        for (int b = 0; b < size; ++b) {
            result = 0;
            for (int x = 0; x < size; ++x) {
                if ((arr[x] ^ arr[x ^ a]) == b) {
                    ++result;
                }
            }
            if (result > max) {
                max = result;
            }
        }

    }
    return max;
}

//Функція генерації S-Box'у за допомогою методу Рою Часток

int *particleSwarmOptimization(int size, int count, int N, int maxIter, int mode, int *finalIter){
    int flag102 = 0;
    clock_t start = clock();
    srand(time(NULL));
    int iter = 0;
    int flag = rand()%size;
    int population[2*N][size];
    for (int i = 0; i < size; ++i){
        population[0][i] = aesSbox[i];
    }
    for (int q = 1; q < N; ++q){
        int *ar1 = SBoxGeneratingDec(count,count,q+flag);
        for(int w = 0; w < size; ++w) {
            population[q][w] = ar1[w];
        }
        free(ar1);
    }
    int arrNL[N];
    for (int q = 0; q < N; ++q){
        for(int w = 0; w < size; ++w){
            printf("%d ",population[q][w]);
        }
        int LAT = LATMax(population[q],size,count);
        int NL = raiseToPower(2, count - 1) - LAT;
        printf( "\nNon-linearity from LAT = %d \n", NL);
        printf("\n");
        arrNL[q] = NL;
    }
    int g[size];
    for (int i = 0; i < N; ++i) {
        for (int j = N - 1; j > i; --j) {
            if (arrNL[j] > arrNL[j - 1]) {
                int t = arrNL[j - 1];
                arrNL[j - 1] = arrNL[j];
                arrNL[j] = t;
                for (int k = 0; k < size; ++k) {
                    g[k] = population[j - 1][k];
                    population[j - 1][k] = population[j][k];
                    population[j][k] = g[k];
                }
            }
        }
    }
    printf("\n");
    for (int q = 0; q < N; ++q){
        printf( "\n%d ", arrNL[q]);
    }
    printf("\nSORTED BY Non-Linearity\n");
    for (int q = 0; q < N; ++q){
        for(int w = 0; w < size; ++w){
            printf("%d, ",population[q][w]);
        }
        printf("\n\n");
    }
    int gBest[size];
    for (int m = 0; m < size; ++m){
        gBest[m] = population[0][m];
    }
    int pBest[N][size];
    for (int i = 1; i < N; ++i){
        for (int j = 0; j < size; ++j){
            pBest[i-1][j] = population[i][j];
        }
    }
    printf("\n\n");
    printf("\ngBest\n");
    for (int m = 0; m < size; ++m){
        printf("%d, ",gBest[m]);
    }
    printf("\npBest\n");
    for (int q = 0; q < N-1; ++q){
        for(int w = 0; w < size; ++w){
            printf("%d, ",pBest[q][w]);
        }
        printf("\n\n");
    }
    /*double weight1 = 0.1;
    double weight2 = 1.6;
    double weightCur;
    int curIter = 0;
    int Vel[N][size];
    int arrNLSorted[size];
    while(maxIter > 0) {
        weightCur = weight1+(curIter-1)*((weight2-weight1)/maxIter);
        printf("\nWHILE\n");
        int Q = 100;
        int rd1 = rand() % (Q);
        double xr1 = (double) rd1 / Q;
        double c1 = 2 * xr1;
        int rd2 = rand() % (Q);
        double xr2 = (double) rd2 / Q;
        double c2 = 2 * xr2;
        int rd3 = rand() % (Q);
        double xr3 = (double) rd3 / Q;
        double r1 = xr3;
        int rd4 = rand() % (Q);
        double xr4 = (double) rd4 / Q;
        double r2 = xr4;
        printf("c1 = %lf ", c1);
        printf("c2 = %lf ", c2);
        printf("r1 = %lf ", r1);
        printf("r2 = %lf \n", r2);
        printf("\n\n");
        for (int b = 0; b < N; ++b) {
            arrNLSorted[b] = arrNL[b];
        }
        int tempSbox[size];
        int tempSbox2[size];
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < size;) {
                if (mode == 1){
                    Vel[i][j] = ceil(weightCur * Vel[i][j] + c1 * r1 * (gBest[j] - population[i][j] +
                                                                        c2 * r2 * (gBest[j] - population[i][j])));
                }
                if (mode == 0) {
                    Vel[i][j] = ceil(weightCur * Vel[i][j] + c1 * r1 * (pBest[i][j] - population[i][j] +
                                                                        c2 * r2 * (gBest[j] - population[i][j])));
                }
                if (Vel[i][j] < 0) {
                    Vel[i][j] = myModulusDec((Vel[i][j] + 256), 256);
                }
                //printf("Vel[%d][%d] = %d ", i,j,Vel[i][j]);
                int X = myModulusDec((population[i][j] + Vel[i][j]), 256);
                int contains;
                if (contains == 0) {
                    tempSbox[j] = X;
                } else {
                    tempSbox[j] = myModulusDec((tempSbox[j] + rand()), 256);
                    if (tempSbox[j] == 0){
                        tempSbox[j] = myModulusDec((tempSbox[j]+rand()),256);
                    }
                };
                //printf("\n%d temp box b4 cycle [%d]", tempSbox[j],j);
                contains = 0;
                for (int k = 0; k < j; ++k) {
                    if (tempSbox[k] == tempSbox[j]) {
                        //printf("\n%d k = \n", k);
                        //printf("\n%d j = \n", j);
                        //tempSbox[j] = myModulusDec((tempSbox[j]+rand()),256);
                        //tempSbox2[j] = tempSbox[j];
                        //printf("\n%d j", tempSbox[j]);
                        //printf("\n%d k", tempSbox[k]);
                        contains = 1;
                        break;
                    }
                }
                if (!contains) {
                    j++;
                }
            }
            if (i == 0) {
                int LAT = LATMax(tempSbox, size, count);
                int NL = raiseToPower(2, count - 1) - LAT;
                if (NL > 98) {
                    for (int v = 0; v < 15; ++v) {
                        srand(tempSbox[v] * (curIter * v) % 256);
                        int coeff = rand() % 256;
                        printf("coeff1 %d", coeff);
                        int coeff2 = rand() % 256;
                        printf("coeff2 %d", coeff2);
                        int temp = tempSbox[coeff];
                        tempSbox[coeff] = tempSbox[coeff2];
                        tempSbox[coeff2] = temp;
                    }
                }
            }
            for (int k = 0; k < size; ++k) {
                population[N + i][k] = tempSbox[k];
                //printf("%d ", tempSbox[k]);
            }
            //printf("\n");
        }
        printf("\nNEW Arrays\n");
        for (int q = 0; q < 2*N; ++q){
            for(int w = 0; w < size; ++w){
                printf("%d, ",population[q][w]);
            }
            int LAT = LATMax(population[q],size,count);
            int NL = raiseToPower(2, count - 1) - LAT;
            printf( "\nNon-linearity from LAT = %d \n", NL);
            printf("\n");
            printf("\n\n");
        }
        int arrNL2[2 * N];
        for (int q = 0; q < 2 * N; ++q) {
            for (int w = 0; w < size; ++w) {
                printf("%d ", population[q][w]);
            }
            int LAT = LATMax(population[q], size, count);
            int NL = raiseToPower(2, count - 1) - LAT;
            printf("\nNon-linearity from LAT = %d \n", NL);
            printf("\n");
            arrNL2[q] = NL;
        }
        int g2[size];
        for (int i = 0; i < 2 * N; ++i) {
            for (int j = (2 * N - 1); j > i; --j) {
                if (arrNL2[j] > arrNL2[j - 1]) {
                    int h = arrNL2[j - 1];
                    arrNL2[j - 1] = arrNL2[j];
                    arrNL2[j] = h;
                    for (int k = 0; k < size; ++k) {
                        g2[k] = population[j - 1][k];
                        population[j - 1][k] = population[j][k];
                        population[j][k] = g2[k];
                    }
                }
            }
        }
        printf("\n");
        for (int q = 0; q < 2 * N; ++q) {
            printf("\n%d ", arrNL2[q]);
        }
        for (int q = N; q < 2 * N; ++q) {
            for (int w = 0; w < size; ++w) {
                population[q][w] = 0;
            }
        }
        printf("\nSORTED BY Non-Linearity\n");
        for (int q = 0; q < N; ++q) {
            for (int w = 0; w < size; ++w) {
                printf("%d, ", population[q][w]);
            }
            printf("\n\n");
        }
        if (curIter == 0) {
            for (int m = 0; m < size; ++m) {
                gBest[m] = population[0][m];
            }
            for (int i = 1; i < N; ++i) {
                for (int j = 0; j < size; ++j) {
                    pBest[i - 1][j] = population[i][j];
                }
            }
        }
        else {
            for (int m = 0; m < size; ++m) {
                gBest[m] = population[1][m];
            }
            for (int i = 2; i < N; ++i) {
                for (int j = 0; j < size; ++j) {
                    pBest[i - 1][j] = population[i][j];
                }
            }
        }
        printf("\n\n");
        printf("\ngBest\n");
        for (int m = 0; m < size; ++m){
            printf("%d, ",gBest[m]);
        }
        printf("\npBest\n");
        for (int q = 0; q < N-1; ++q){
            for(int w = 0; w < size; ++w){
                printf("%d, ",pBest[q][w]);
            }
            printf("\n\n");
        }
        for (int h = 1; h < N; ++h){
            int LAT3 = LATMax(population[h], size, count);
            int NL3 = raiseToPower(2, count - 1) - LAT3;
            printf("\nNon-linearity from LAT = %d \n", NL3);
            printf("\n");
            int ucCheck3 = linearRedundancy(population[h], 256, 8);
            int ai3 = algebraicImmunity(population[h],256,8);
            int du3 = deltaUniformity(population[h],256,8);
            printf( "\n\nImmunity   = %d \n", ai3);
            printf("\nDelta-Uniformity  = %d \n", du3);
            int lr;
            if (ucCheck3 == 1) {
                printf("\nLinear redundancy = %d \n", (256) - ucCheck3);
                lr = 256-ucCheck3;
            } else {
                printf("\nLinear redundancy = %d \n", (256 - 1) - ucCheck3);
                lr = (256 - 1) - ucCheck3;
            }
            if (NL3 >= 104 && ai3 == 3 && lr == 0){
                printf("current iter = %d ",curIter+1);
                *finalIter = curIter+1;
                maxIter = -99;
                break;
            }
            if (flag102 == 0 && NL3 == 102 && ai3 == 3 && lr == 0){
                FILE *fileLocal;
                fopen_s(&fileLocal, "Table results N, MaxIter, IterToFind, Time.txt", "a");
                if (fileLocal == NULL) {
                    printf("ERROR: Can't save sbox to file!\n");
                    for (;;);
                }
                clock_t finish = clock();
                double a = (double) (finish - start) / CLOCKS_PER_SEC;
                fprintf(fileLocal, "\n%d       %d       %d              %d        %d,%d", NL3,N,500,curIter+1,(int)a,(int)((-1)*(floor(a)-a)*1000000));
                printf("current iter = %d ",curIter+1);
                fclose(fileLocal);
                flag102 = 1;
            }
        }
        maxIter = maxIter-1;
        mode = 0;
        ++curIter;
    }*/
    for (int h = 0; h < N; ++h) {
        int LAT3 = LATMax(population[h], size, count);
        int NL3 = raiseToPower(2, count - 1) - LAT3;
        printf("\nNon-linearity from LAT = %d \n", NL3);
        int ucCheck3 = linearRedundancy(population[h], 256, 8);
        int ai3 = algebraicImmunity(population[h], 256, 8);
        int du3 = deltaUniformity(population[h], 256, 8);
        printf("\nImmunity   = %d \n", ai3);
        printf("\nDelta-Uniformity  = %d \n", du3);
        int lr;
        if (ucCheck3 == 1) {
            printf("\nLinear redundancy = %d \n", (256) - ucCheck3);
            lr = 256 - ucCheck3;
        } else {
            printf("\nLinear redundancy = %d \n", (256 - 1) - ucCheck3);
            lr = (256 - 1) - ucCheck3;
        }
        printf("\n");
    }
    //printf("\n\nFinal data\n\n");
    int *result = (int *) calloc((N*size), sizeof(int));
    //printf("final current iter = %d ",*finalIter);
    for (int q = 0; q < N; ++q) {
        for (int w = 0; w < size; ++w) {
            result[q * size + w] = population[q][w];
            //printf("%d, ", result[q * size + w]);
        }
        //printf("\n\n");
    }
    return result;
}

