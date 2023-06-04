typedef struct IUnknown IUnknown;
#include "fstream"
#include "iostream"
#include "iomanip"
#include "cmath"
#include "ctime"
#include "clocale"
#include "algorithm"

float** normalize(float** alt, int nA, int nC, bool* up, int type) {
    float** norm = new float* [nC];
    for (int i = 0; i < nC; i++)
        norm[i] = new float[nA]();

    for (int i = 0; i < nC; i++) {
        if (type == 1) {
            float* maxA = std::max_element(alt[i], alt[i] + nA);
            float* minA = std::min_element(alt[i], alt[i] + nA);
            if (up[i])
                for (int j = 0; j < nA; j++)
                    norm[i][j] = (alt[i][j] - *minA) / (*maxA - *minA);
            else
                for (int j = 0; j < nA; j++)
                    norm[i][j] = (*maxA - alt[i][j]) / (*maxA - *minA);
        }
        else if (type == 2) {
            if (up[i]) {
                float* maxA = std::max_element(alt[i], alt[i] + nA);
                float* minA = std::min_element(alt[i], alt[i] + nA);
                float delta = *minA <= 0 ? std::abs(*minA) + 1 : 0;
                for (int j = 0; j < nA; j++)
                    norm[i][j] = (alt[i][j] + delta) / (*maxA + delta);
            }
            else {
                float* minA = std::min_element(alt[i], alt[i] + nA);
                float delta = *minA <= 0 ? std::abs(*minA) + 1 : 0;
                for (int j = 0; j < nA; j++) {
                    norm[i][j] = (*minA + delta) / (alt[i][j] + delta);
                }
            }
        }
    }
    return norm;
}


void multi_sv(float** data, int rows, int cols, float* weights) {
    float maxVal = 0;
    int maxIndex = 0;

    std::cout << " Multiplicative Convolution\n     |";
    for (int i = 0; i < cols; i++) {
        std::cout << "  C" << i + 1 << " |";
    }
    std::cout << " Sum\n";

    for (int i = 0; i < rows; i++) {
        std::cout << "Row" << i + 1 << "  |";
        float sum = 1;

        for (int j = 0; j < cols; j++) {
            float value = data[j][i];
            std::cout << " " << value << " |";
            sum *= std::pow(value, weights[j]);
        }

        std::cout << " " << sum << "\n";

        if (sum > maxVal) {
            maxVal = sum;
            maxIndex = i;
        }
    }

    std::cout << "Weights |";
    for (int j = 0; j < cols; j++) {
        std::cout << " " << weights[j] << " |";
    }
    std::cout << "\n\nThe optimal solution is alternative A" << maxIndex + 1 << "\n\n";
}

void addit_sv(float** data, int rows, int cols, float* weights) {
    float maxVal = 0;
    int maxIndex = 0;
    std::cout << " Additive Convolution\n     |";
    for (int i = 0; i < cols; i++) {
        std::cout << "  C" << i + 1 << " |";
    }
    std::cout << " Sum\n";
    for (int i = 0; i < rows; i++) {
        std::cout << "Row" << i + 1 << "  |";
        float sum = 0;
        for (int j = 0; j < cols; j++) {
            float value = data[j][i];
            std::cout << " " << value << " |";
            sum += value * weights[j];
        }
        std::cout << " " << sum << "\n";
        if (sum > maxVal) {
            maxVal = sum;
            maxIndex = i;
        }
    }
    std::cout << "Weights |";
    for (int j = 0; j < cols; j++) {
        std::cout << " " << weights[j] << " |";
    }
    std::cout << "\n\nThe optimal solution is alternative A" << maxIndex + 1 << "\n\n";
}

void maximin(float** data, int rows, int cols) {
    float maxVal = 0;
    int maxIndex = 0;
    std::cout << "  Maximin Method\n     |";
    for (int i = 0; i < cols; i++) {
        std::cout << " Column" << i + 1 << " |";
    }
    std::cout << " Min\n";
    for (int i = 0; i < rows; i++) {
        std::cout << "Row" << i + 1 << "  |";
        int minIndex = 0;
        float minVal = data[minIndex][i];

        for (int j = 0; j < cols; j++) {
            float value = data[j][i];
            std::cout << " " << value << " |";
            if (value < minVal) {
                minVal = value;
                minIndex = j;
            }
        }
        std::cout << " " << minVal << "\n";
        if (minVal > maxVal) {
            maxVal = minVal;
            maxIndex = i;
        }
    }
    std::cout << "\n\nThe optimal solution is alternative A" << maxIndex + 1 << "\n\n";
}

void celProg(float** norm, int nA, int nC, float p) {
    float minim = 10;
    int minI = 0;
    float* f = new float[nC]();

    // Вычисление максимальных значений для каждого столбца
    for (int i = 0; i < nC; i++) {
        float maxVal = norm[i][0];
        for (int j = 1; j < nA; j++) {
            if (norm[i][j] > maxVal) {
                maxVal = norm[i][j];
            }
        }
        f[i] = maxVal;
    }

    std::cout << " Целевое программирование\n     |";
    for (int i = 0; i < nC; i++) {
        std::cout << "  C" << i + 1 << " |";
    }
    std::cout << " sumA\n";

    // Вычисление суммы и определение оптимального решения
    for (int i = 0; i < nA; i++) {
        std::cout << "A" << i + 1 << "  |";
        float sum = 0;
        for (int j = 0; j < nC; j++) {
            std::cout << " " << norm[j][i] << " |";
            sum += std::pow(std::abs(norm[j][i] - f[j]), p);
        }
        sum = std::pow(sum, 1.0 / p);
        std::cout << " " << sum << "\n";
        if (sum < minim) {
            minim = sum;
            minI = i;
        }
    }

    std::cout << "\n\nОптимальным решением является альтернатива А" << minI + 1 << "\n\n";
    delete[] f;
}


void mainCrit(float** alt, int nA, int nC, bool* up, int mainC) {
    std::cout << " Метод главного критерия\n      |";
    for (int i = 0; i < nC; i++)
        std::cout << (up[i] ? "    max  " : "    min  ");
    std::cout << "\n      |";
    for (int i = 0; i < nC; i++) {
        std::cout << "     C" << i + 1 << "  ";
    }
    for (int i = 0; i < nA; i++) {
        std::cout << "\n A" << i + 1 << "  |";
        for (int j = 0; j < nC; j++)
            std::cout << std::setw(7) << alt[j][i] << " |";
    }
    std::cout << "\n Опт  |";
    float* C = new float[nC - 1]();
    for (int i = 0, j = 0; i < nC; i++) {
        if (i != mainC - 1) {
            float ma = *std::max_element(alt[i], alt[i] + nA);
            float mi = *std::min_element(alt[i], alt[i] + nA);
            if (up[i]) {
                std::cout << std::setw(7) << ma << " |";
                C[j] = ma - (ma - mi) / 3;
            }
            else {
                std::cout << std::setw(7) << mi << " |";
                C[j] = mi + (ma - mi) / 3;
            }
            j++;
        }
        else {
            std::cout << "        |";
        }
    }
    int res = -1;
    while (res == -1) {
        float glCr = up[mainC - 1] ? -999999999 : 999999999;
        for (int i = 0; i < nA; i++) {
            bool correct = true;
            for (int j = 0, h = 0; j < nC && correct; j++) {
                if (j != mainC - 1) {
                    if ((up[j] && C[h] > alt[j][i]) || (!up[j] && C[h] < alt[j][i])) {
                        correct = false;
                    }
                    h++;
                }
            }
            if (correct && ((up[mainC - 1] && glCr < alt[mainC - 1][i]) || (!up[mainC - 1] && glCr > alt[mainC - 1][i]))) {
                glCr = alt[mainC - 1][i];
                res = i;
            }
        }

        for (int i = 0, j = 0; i < nC; i++) {
            if (i != mainC - 1) {
                if (up[i]) {
                    std::cout << "\nC" << i + 1 << " >= " << C[j];
                    float mi = *std::min_element(alt[i], alt[i] + nA);
                    C[j] -= (C[j] - mi) / 3;
                }
                else {
                    std::cout << "\nC" << i + 1 << " <= " << C[j];
                    float ma = *std::max_element(alt[i], alt[i] + nA);
                    C[j] += (ma - C[j]) / 3;
                }
                j++;
            }
        }
        std::cout << "\n";
    }
    std::cout << "\nОптимальным решением является альтернатива А" << res + 1 << "\n\n";
    delete[] C;
}

void metodustupok(float** alt, int numAlt, int numCrit, float* priority, bool* isUpward) {
    float* tmp = new float[numCrit];
    int* critInd = new int[numCrit];

    std::cout << "Метод уступок\n      |";
    for (int i = 0; i < numCrit; i++) {
        std::cout << "   C" << i + 1 << "   ";
    }
    for (int i = 0; i < numAlt; i++) {
        std::cout << "\n A" << i + 1 << "  |";
        for (int j = 0; j < numCrit; j++)
            std::cout << std::setw(7) << alt[j][i] << " |";
    }
    std::cout << "\n\n";

    for (int i = 0; i < numCrit; i++)
        tmp[i] = priority[i];

    for (int i = 0; i < numCrit; i++) {
        auto maxIter = std::max_element(tmp, tmp + numCrit);
        int maxIdx = std::distance(tmp, maxIter);
        critInd[i] = maxIdx;
        tmp[maxIdx] = -1;
    }
    delete[] tmp;

    bool* result = new bool[numAlt];
    for (int i = 0; i < numAlt; i++)
        result[i] = true;

    int chosenAlt, unmatchedCount;
    for (int i = 0; i < numCrit; i++) {
        auto minIter = std::min_element(alt[critInd[i]], alt[critInd[i]] + numAlt);
        auto maxIter = std::max_element(alt[critInd[i]], alt[critInd[i]] + numAlt);
        int minIdx = std::distance(alt[critInd[i]], minIter);
        int maxIdx = std::distance(alt[critInd[i]], maxIter);


        if (isUpward[critInd[i]])
            std::cout << "C" << critInd[i] + 1 << " -> max = " << alt[critInd[i]][maxIdx] << "\n";
        else
            std::cout << "C" << critInd[i] + 1 << " -> min = " << alt[critInd[i]][minIdx] << "\n";

        float metust = (alt[critInd[i]][maxIdx] - alt[critInd[i]][minIdx]) / 2.5;
        do {
            metust *= 1.2;
            std::cout << "Назначим уступку Z" << i + 1 << " = " << metust << "\n";
            if (isUpward[critInd[i]])
                std::cout << "C" << critInd[i] + 1 << " >= " << alt[critInd[i]][maxIdx] - metust << " ( ";
            else
                std::cout << "C" << critInd[i] + 1 << " <= " << alt[critInd[i]][minIdx] + metust << " ( ";

            unmatchedCount = 0;
            for (int j = 0; j < numAlt; j++) {
                if (isUpward[critInd[i]]) {
                    if (alt[critInd[i]][j] >= alt[critInd[i]][maxIdx] - metust && result[j]) {
                        std::cout << "A" << j + 1 << " ";
                        chosenAlt = j;
                    }
                    else
                        unmatchedCount++;
                }
                else {
                    if (alt[critInd[i]][j] <= alt[critInd[i]][minIdx] + metust && result[j]) {
                        std::cout << "A" << j + 1 << " ";
                        chosenAlt = j;
                    }
                    else
                        unmatchedCount++;
                }
            }
            std::cout << ")\n";
        } while (unmatchedCount == numAlt);

        if (unmatchedCount == numAlt - 1) {
            std::cout << "\nОптимальным решением является альтернатива А" << chosenAlt + 1 << "\n\n";
            break;
        }

        for (int j = 0; j < numAlt; j++) {
            if (isUpward[critInd[i]]) {
                if (alt[critInd[i]][j] < alt[critInd[i]][maxIdx] - metust)
                    result[j] = false;
            }
            else {
                if (alt[critInd[i]][j] > alt[critInd[i]][minIdx] + metust)
                    result[j] = false;
            }
        }
    }

    if (numAlt - unmatchedCount >= 2) {
        int chosenAltIdx = chosenAlt;
        for (int i = 0; i < numAlt; i++) {
            if (isUpward[critInd[0]]) {
                if (result[i] && alt[critInd[0]][chosenAltIdx] < alt[critInd[0]][i])
                    chosenAltIdx = i;
            }
            else {
                if (result[i] && alt[critInd[0]][chosenAltIdx] > alt[critInd[0]][i])
                    chosenAltIdx = i;
            }
        }
        std::cout << "\nОптимальным решением является альтернатива А" << chosenAltIdx + 1 << "\n\n";
    }

    delete[] result;
    delete[] critInd;
}

void printResult(float** alt, int numAlt, int numCrit) {
    for (int i = 0; i < numCrit; i++) {
        std::cout << "   C" << i + 1 << "  ";
    }
    for (int i = 0; i < numAlt; i++) {
        std::cout << " \n A" << i + 1 << "  |";
        for (int j = 0; j < numCrit; j++)
        {
            std::cout << std::setw(7) << alt[j][i] << " |";
        }
    }
    std::cout << "\n\n";
}

void analiz(float** alt, int numAlt, int numCrit, bool* isUpward, float* priority) {
    std::cout << "Результат:\n      |";
    printResult(alt, numAlt, numCrit);

    float** normalized12 = normalize(alt, numAlt, numCrit, isUpward, 1);
    float** normalized34 = normalize(alt, numAlt, numCrit, isUpward, 2);

    std::cout << "Нормализация (1,2):\n      |";
    printResult(normalized12, numAlt, numCrit);
    std::cout << "Нормализация (3,4):\n      |";
    printResult(normalized34, numAlt, numCrit);

    addit_sv(normalized12, numAlt, numCrit, priority);
    multi_sv(normalized34, numAlt, numCrit, priority);
    maximin(normalized34, numAlt, numCrit);
    celProg(normalized34, numAlt, numCrit, 2);
    mainCrit(alt, numAlt, numCrit, isUpward, 1);
    metodustupok(alt, numAlt, numCrit, priority, isUpward);

    // Освобождение памяти
    delete[] * normalized12;
    delete[] normalized12;
    delete[] * normalized34;
    delete[] normalized34;
}

float ravnrasp(float a, float b) {
    return a + ((float)rand() / RAND_MAX) * (b - a);
}

float normrasp(float mean, float sr_otkl) {
    float lowerBound = (mean - sqrt(3) * sr_otkl) / 10;
    float upperBound = (mean + sqrt(3) * sr_otkl) / 10;
    float value = 0;

    for (int i = 0; i < 10; i++) {
        value += lowerBound + (upperBound - lowerBound) * rand() / RAND_MAX;
    }
    return value;
}

float exponrasp(float mean) {
    float r = (float)rand() / RAND_MAX;
    r = (r < 0.00000000001) ? 0.00000000001 : r;
    return -mean * log(r);
}

float generate_polomka() {
    float p = (float)rand() / RAND_MAX;
    return (p < 0.05) ? exponrasp(30) : 0;
}

int minP(float* mas, int n) {
    int m = 0;
    for (int i = 1; i < n; i++) {
        if (mas[i] < mas[m]) {
            m = i;
        }
    }
    return m;
}

void MinMag(float(*mas)[2], int n, int m, int* tekMag, int* tekTovar) {
    float min = 10000000;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (mas[i][j] < min) {
                min = mas[i][j];
                *tekMag = i;
                *tekTovar = j;
            }
        }
    }
}

void opred_fazy(float** array, int rows, int columns, int* tekStolbec, int* tekStroka) {
    float minValue = 1000000;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            if (array[i][j] < minValue) {
                minValue = array[i][j];
                *tekStroka = i;
                *tekStolbec = j;
            }
        }
    }
}

int main()
{
    setlocale(LC_ALL, "");
    srand(time(NULL));

    std::ofstream file("delivery.txt");
    if (!file.is_open()) {
        std::cerr << "Не удалось открыть файл!\n";
        return -1;
    }

    std::streambuf* stream_buffer_cout = std::cout.rdbuf();
    std::cout.rdbuf(file.rdbuf());
    std::cout.precision(3);
    std::cout << std::fixed;

    int f = 0, fl = 0, tekAuto = 0, currPogr = 0, tekMag = 0, tekTovar = 0, nCar = 0;
    const int count = 11;
    float minn = 0, prostoiSr = 0, prostoiM = 0, prSt[count] = { 90000000 }, prM[count] = { 9000000 },
        zapastov[6][2] = { 0 }, zapastov_copy[6][2] = { 0 }, pogr_osv[2] = { 0 };
    const float infin = 900000000;

    // цикл по количеству машин
    for (nCar = 2; nCar < count; nCar++) {
        f = 0;
        for (int num_mag = 0; num_mag < 6; num_mag++)
            for (int tip = 0; tip < 2; tip++) {
                zapastov[num_mag][tip] = 0;
                zapastov_copy[num_mag][tip] = 0;
            }
        float** vehicles = new float* [2]; // массив на две строки: 1 - приезд на хлебозавод, 2 - приезд в магазин

        // цикл по количеству поездок
        for (int s = 0; s < 2; s++) {
            vehicles[s] = new float[nCar];
        }
        for (int s = 0; s < nCar; s++) {
            vehicles[1][s] = infin;
            vehicles[0][s] = 0;
        }
        int* type_p = new int[nCar] {}; // Массив для хранения типа поездки для каждого автомобиля
        int** mag = new int* [2]; // Массив для хранения номеров магазинов для каждого автомобиля
        int* tovar = new int[nCar]; // Массив для хранения типа товара для каждого автомобиля

        for (int yt = 0; yt < 2; yt++)
            mag[yt] = new int[nCar];

        prostoiSr = 0; prostoiM = 0; pogr_osv[0] = 0; pogr_osv[1] = 0; //время освобождения мест на погрузку
        fl = 0; tekAuto = 0; currPogr = 0; tekMag = 0, tekTovar = 0; minn = 0;
        // по количеству поездок
        std::cout << "Задействовано " << nCar << " машин\n________________________________________________\n";
        int i = 0;
        while (true)
        {
            // выбираем фазу путем поиска минимального числа в массиве vehicles
            opred_fazy(vehicles, 2, nCar, &tekAuto, &f); // нашли минимальный элемент, номер строки - фаза, номер столбца - номер машины
            if (f == 0) {
                // поездка в магазин
                float r = (float)rand() / RAND_MAX;
                if (r < 0.4) type_p[tekAuto] = 1;// тип 1: поездка в один магазин, которому раньше всех нужен один из товаров
                else if (r < 0.6) type_p[tekAuto] = 2;// тип 2: поездка в два магазина, которым в скорем времени нужен один товар (разница между привозом не больше 8 часов (480 минут))
                else type_p[tekAuto] = 3; // тип 3: поездка с двумя товарами в один магазин

                if (type_p[tekAuto] == 1) {
                    std::cout << std::setw(2) << type_p[tekAuto] << " Едет машина " << tekAuto + 1 << "---> " << vehicles[f][tekAuto] << "\n";
                    // выбираем магазин и товар
                    MinMag(zapastov, 6, 2, &tekMag, &tekTovar);

                    double remainingTime = zapastov[tekMag][tekTovar] - vehicles[f][tekAuto];

                    if ((infin - zapastov[tekMag][tekTovar]) < 0.000001)
                    {
                        vehicles[f][tekAuto] += 100;
                        std::cout << "\n1 Погрузка данного типа невозможна\n";
                        break;
                    }
                    // если до нужного времени прибытия больше времени одной поездки, то машина начинает обслуживаться с нового времени
                    if (remainingTime > 60) vehicles[f][tekAuto] = zapastov[tekMag][tekTovar] - (10 + 30 + 20);
                    std::cout << "\nВыбран " << tekMag + 1 << " магазин, " << tekTovar << " товар";
                    std::cout << "\nМашина готова к погрузке с " << vehicles[f][tekAuto] << "мин.\n";
                    mag[0][tekAuto] = tekMag;
                    tovar[tekAuto] = tekTovar;
                    currPogr = minP(pogr_osv, 2); // выбираем место погрузки

                    // если место занято, то машина обслуживается с его освобождения
                    if (pogr_osv[currPogr] > vehicles[f][tekAuto]) {
                        vehicles[f][tekAuto] = pogr_osv[currPogr];
                        std::cout << "Место занято, машина обслуживается с: " << vehicles[f][tekAuto] << "\n";
                    }
                    // если машина начнет обслуживание после необходимого срока, то прекращение работы
                    if (vehicles[f][tekAuto] > 10080) {
                        fl = 1;
                        break;
                    }
                    // если место было свободно, а машины не было долго, то увеличиваем простой мест
                    if (vehicles[f][tekAuto] > pogr_osv[currPogr] && pogr_osv[currPogr] > 0) prostoiM += vehicles[f][tekAuto] - pogr_osv[currPogr];
                    std::cout << "Машина становится на погрузку (место " << currPogr + 1 << ") в " << vehicles[f][tekAuto] << "мин.\n";
                    vehicles[f][tekAuto] += normrasp(10, 3); // погрузка
                    std::cout << "Погрузка завершена " << vehicles[f][tekAuto] << "мин.\n";
                    // место погрузки снова свободно с нового времени
                    pogr_osv[currPogr] = vehicles[f][tekAuto];
                    vehicles[f][tekAuto] += exponrasp(30); // движение
                    vehicles[f][tekAuto] += generate_polomka(); // если сломалась, то + ремонт
                    std::cout << "Прибытие в магазин " << vehicles[f][tekAuto] << "мин.\n";
                    vehicles[1][tekAuto] = vehicles[f][tekAuto];
                    vehicles[f][tekAuto] = infin;
                    zapastov_copy[tekMag][tekTovar] = zapastov[tekMag][tekTovar];
                    zapastov[tekMag][tekTovar] = infin;
                    std::cout << "------------------------------------- \n";
                }
                else if (type_p[tekAuto] == 2) {
                    std::cout << std::setw(2) << type_p[tekAuto] << " Едет машина " << tekAuto + 1 << "---> " << vehicles[f][tekAuto] << "\n";
                    // сразу выбираем два магазина
                    // выбираем первый магазин и общий товар
                    MinMag(zapastov, 6, 2, &tekMag, &tekTovar);
                    // выберем второй магазин
                    minn = zapastov[0][tekTovar];
                    mag[1][tekAuto] = 0;
                    for (int e = 1; e < 6; e++)
                        if (zapastov[e][tekTovar] <= minn && e != tekMag) {
                            minn = zapastov[e][tekTovar];
                            mag[1][tekAuto] = e;
                        }
                    if ((infin - zapastov[tekMag][tekTovar]) < 0.000001 || (infin - zapastov[mag[1][tekAuto]][tekTovar]) < 0.000001) {
                        vehicles[f][tekAuto] += 100;
                        std::cout << "\n2 Погрузка данного типа невозможна\n"; break;
                    }
                    if (zapastov[tekMag][tekTovar] - vehicles[f][tekAuto] > 65)
                        vehicles[f][tekAuto] = zapastov[tekMag][tekTovar] - (15 + 30 + 20);
                    std::cout << "\nВыбраны " << tekMag + 1 << ", " << mag[1][tekAuto] + 1 << " магазины, " << tekTovar << " товар";
                    std::cout << "\nМашина готова к погрузке с " << vehicles[f][tekAuto] << "мин.\n";
                    mag[0][tekAuto] = tekMag;
                    tovar[tekAuto] = tekTovar;
                    currPogr = minP(pogr_osv, 2);
                    if (pogr_osv[currPogr] > vehicles[f][tekAuto]) {
                        vehicles[f][tekAuto] = pogr_osv[currPogr];
                        std::cout << "Место занято, машина обслуживается с: " << vehicles[f][tekAuto] << "\n";
                    }
                    if (vehicles[f][tekAuto] > 10080) { fl = 1; break; }
                    if (vehicles[f][tekAuto] - pogr_osv[currPogr] > 0 && pogr_osv[currPogr] > 0)  prostoiM += vehicles[f][tekAuto] - pogr_osv[currPogr];
                    std::cout << "Машина становится на погрузку (место " << currPogr + 1 << ") в " << vehicles[f][tekAuto] << "мин.\n";
                    vehicles[f][tekAuto] += normrasp(15, 5); // погрузка
                    std::cout << "Погрузка завершена " << vehicles[f][tekAuto] << "мин.\n";
                    pogr_osv[currPogr] = vehicles[f][tekAuto];
                    vehicles[f][tekAuto] += exponrasp(30);  // движение 
                    vehicles[f][tekAuto] += generate_polomka(); // если сломалась то + ремонт
                    std::cout << "Прибытие в магазин " << vehicles[f][tekAuto] << "мин.\n";
                    vehicles[1][tekAuto] = vehicles[f][tekAuto];
                    vehicles[f][tekAuto] = infin;
                    zapastov_copy[mag[0][tekAuto]][tekTovar] = zapastov[mag[0][tekAuto]][tekTovar];
                    zapastov_copy[mag[1][tekAuto]][tekTovar] = zapastov[mag[1][tekAuto]][tekTovar];
                    zapastov[mag[0][tekAuto]][tekTovar] = infin;
                    zapastov[mag[1][tekAuto]][tekTovar] = infin;
                    std::cout << "------------------------------------- \n";
                }
                else if (type_p[tekAuto] == 3) {
                    // выбираем магазин
                    minn = 100000;
                    for (int r = 0; r < 6; r++) {
                        if ((zapastov[r][0] <= minn || zapastov[r][1] <= minn) && fabs(zapastov[r][1] - zapastov[r][0]) < 480) {
                            minn = zapastov[r][0] < zapastov[r][1] ? zapastov[r][0] : zapastov[r][1];
                            tekMag = r;
                        }
                    }
                    if ((infin - zapastov[tekMag][0]) < 0.000001 || (infin - zapastov[tekMag][1]) < 0.000001) { vehicles[f][tekAuto] += 100; std::cout << "\n3 Погрузка данного типа невозможна\n"; break; }
                    std::cout << std::setw(2) << type_p[tekAuto] << " Едет машина " << tekAuto + 1 << "---> " << vehicles[f][tekAuto];
                    std::cout << "\n";
                    if (zapastov[tekMag][0] - vehicles[f][tekAuto] > 75 && zapastov[tekMag][1] - vehicles[f][tekAuto] > 75)
                        vehicles[f][tekAuto] = (zapastov[tekMag][0] < zapastov[tekMag][1] ? zapastov[tekMag][0] : zapastov[tekMag][1]) - (25 + 30 + 20);
                    std::cout << "\nВыбран " << tekMag + 1 << " магазин, " << "0, 1 товар";
                    std::cout << "\nМашина готова к погрузке с " << vehicles[f][tekAuto] << "мин.\n";
                    mag[0][tekAuto] = tekMag;
                    tovar[tekAuto] = tekTovar;
                    currPogr = minP(pogr_osv, 2);
                    if (pogr_osv[currPogr] > vehicles[f][tekAuto]) {
                        vehicles[f][tekAuto] = pogr_osv[currPogr];
                        std::cout << "Место занято, машина обслуживается с: " << vehicles[f][tekAuto] << "\n";
                    }
                    if (vehicles[f][tekAuto] > 10080) { fl = 1; break; }
                    if (vehicles[f][tekAuto] > pogr_osv[currPogr] && pogr_osv[currPogr] > 0)
                        prostoiM += vehicles[f][tekAuto] - pogr_osv[currPogr];
                    std::cout << "Машина становится на погрузку (место " << currPogr + 1 << ") в " << vehicles[f][tekAuto] << "мин.\n";
                    vehicles[f][tekAuto] += normrasp(25, 7); // погрузка
                    std::cout << "Погрузка завершена " << vehicles[f][tekAuto] << "мин.\n";
                    pogr_osv[currPogr] = vehicles[f][tekAuto];
                    vehicles[f][tekAuto] += exponrasp(30);  // движение 
                    vehicles[f][tekAuto] += generate_polomka(); // если сломалась то + ремонт
                    std::cout << "Прибытие в магазин " << vehicles[f][tekAuto] << "мин.\n";
                    vehicles[1][tekAuto] = vehicles[f][tekAuto];
                    vehicles[f][tekAuto] = infin;
                    zapastov_copy[tekMag][0] = zapastov[tekMag][0];
                    zapastov_copy[tekMag][1] = zapastov[tekMag][1];
                    zapastov[tekMag][0] = infin;
                    zapastov[tekMag][1] = infin;
                    std::cout << "------------------------------------- \n";
                }
            }
            else if (f == 1) {
                // 2 фаза: разгрузка и поездка назад
                if (type_p[tekAuto] == 1) {
                    std::cout << std::setw(2) << type_p[tekAuto] << " Разгружается машина " << tekAuto + 1 << "---> " << vehicles[f][tekAuto] << "\n";
                    vehicles[f][tekAuto] += exponrasp(20); // разгрузка
                    // если машина приехала позже нужного времени, то увеличиваем простой магазинов
                    std::cout << "Магазин " << mag[0][tekAuto] + 1 << ", товар " << tovar[tekAuto] << "\n";
                    std::cout << "Доставка товара в " << zapastov_copy[mag[0][tekAuto]][tovar[tekAuto]] << " мин.\n";
                    if (vehicles[f][tekAuto] - zapastov_copy[mag[0][tekAuto]][tovar[tekAuto]] > 0 && zapastov_copy[mag[0][tekAuto]][tovar[tekAuto]] > 0) {
                        prostoiSr += vehicles[f][tekAuto] - zapastov_copy[mag[0][tekAuto]][tovar[tekAuto]];
                        std::cout << "Машина опоздала на " << vehicles[f][tekAuto] - zapastov_copy[mag[0][tekAuto]][tovar[tekAuto]] << "мин.\n";
                    }
                    std::cout << "Машина приехала в " << vehicles[f][tekAuto] << "мин. \n" << "Текущий простой " << prostoiSr << " мин.\n";
                    // новая поставка через 8-12 часов
                    zapastov[mag[0][tekAuto]][tovar[tekAuto]] = vehicles[f][tekAuto] + ravnrasp(8 * 60, 12 * 60);
                    std::cout << "Следующая поставка: " << zapastov[mag[0][tekAuto]][tovar[tekAuto]] << "\n";
                    vehicles[f][tekAuto] += exponrasp(20);  // движение
                    vehicles[f][tekAuto] += generate_polomka(); // если сломалась то + ремонт
                    std::cout << "Машина вернулась на хлебозавод в " << vehicles[f][tekAuto] << " мин.\n";
                    vehicles[0][tekAuto] = vehicles[f][tekAuto];
                    vehicles[f][tekAuto] = infin;
                    std::cout << "------------------------------------- \n";
                }
                else if (type_p[tekAuto] == 2) {
                    std::cout << std::setw(2) << type_p[tekAuto] << " Разгружается машина " << tekAuto + 1 << "---> " << vehicles[f][tekAuto] << "\n";
                    vehicles[f][tekAuto] += exponrasp(20); // разгрузка
                    std::cout << "Магазины " << mag[0][tekAuto] + 1 << ", " << mag[1][tekAuto] + 1 << ", товар " << tovar[tekAuto] << "\n";
                    std::cout << "Доставка товара в " << zapastov_copy[mag[0][tekAuto]][tovar[tekAuto]] << " мин.\n";
                    if (vehicles[f][tekAuto] - zapastov_copy[mag[0][tekAuto]][tovar[tekAuto]] > 0 && zapastov_copy[mag[0][tekAuto]][tovar[tekAuto]] > 0) {
                        prostoiSr += vehicles[f][tekAuto] - zapastov_copy[mag[0][tekAuto]][tovar[tekAuto]];
                        std::cout << "Машина опоздала на " << vehicles[f][tekAuto] - zapastov_copy[mag[0][tekAuto]][tovar[tekAuto]] << "мин.\n";
                    }
                    std::cout << "Машина приехала в магазин в " << vehicles[f][tekAuto] << "мин. \n" << "Текущий простой " << prostoiSr << " мин.\n";
                    zapastov[mag[0][tekAuto]][tovar[tekAuto]] = vehicles[f][tekAuto] + ravnrasp(8 * 60, 12 * 60);
                    std::cout << "Следующая поставка: " << zapastov[mag[0][tekAuto]][tovar[tekAuto]] << "\n";
                    vehicles[f][tekAuto] += exponrasp(10);  // движение в другой магазин
                    vehicles[f][tekAuto] += generate_polomka(); // если сломалась то + ремонт
                    vehicles[f][tekAuto] += exponrasp(20); // разгрузка
                    std::cout << "Доставка товара в " << zapastov_copy[mag[1][tekAuto]][tovar[tekAuto]] << " мин.\n";
                    if (vehicles[f][tekAuto] - zapastov_copy[mag[1][tekAuto]][tovar[tekAuto]] > 0 && zapastov_copy[mag[1][tekAuto]][tovar[tekAuto]] > 0) {
                        prostoiSr += vehicles[f][tekAuto] - zapastov_copy[mag[1][tekAuto]][tovar[tekAuto]];
                        std::cout << "Машина опоздала на " << vehicles[f][tekAuto] - zapastov_copy[mag[1][tekAuto]][tovar[tekAuto]] << "мин.\n";
                    }
                    std::cout << "Машина приехала в магазин в " << vehicles[f][tekAuto] << "мин. \n" << "Текущий простой " << prostoiSr << " мин.\n";
                    zapastov[mag[1][tekAuto]][tovar[tekAuto]] = vehicles[f][tekAuto] + ravnrasp(8 * 60, 12 * 60);
                    std::cout << "Следующая поставка: " << zapastov[mag[1][tekAuto]][tovar[tekAuto]] << "\n";
                    vehicles[f][tekAuto] += exponrasp(20);  // движение
                    vehicles[f][tekAuto] += generate_polomka(); // если сломалась то + ремонт
                    std::cout << "Машина вернулась на хлебозавод в " << vehicles[f][tekAuto] << " мин.\n";
                    vehicles[0][tekAuto] = vehicles[f][tekAuto];
                    vehicles[f][tekAuto] = infin;
                    std::cout << "------------------------------------- \n";
                }
                else if (type_p[tekAuto] == 3) {
                    std::cout << std::setw(2) << type_p[tekAuto] << " Разгружается машина " << tekAuto + 1 << "---> " << vehicles[f][tekAuto] << "\n";
                    vehicles[f][tekAuto] += exponrasp(20); // разгрузка
                    std::cout << "Магазин " << mag[0][tekAuto] + 1 << ", товары 0, 1" << "\n";
                    std::cout << "Доставка товара в " << zapastov_copy[mag[0][tekAuto]][0] << "мин. ," << zapastov_copy[mag[0][tekAuto]][1] << " мин.\n";
                    if (vehicles[f][tekAuto] - zapastov_copy[mag[0][tekAuto]][0] > 0 && zapastov_copy[mag[0][tekAuto]][0] > 0) {
                        prostoiSr += vehicles[f][tekAuto] - zapastov_copy[mag[0][tekAuto]][0];
                        std::cout << "Машина опоздала с доставкой первого товара на " << vehicles[f][tekAuto] - zapastov_copy[mag[0][tekAuto]][0] << "мин.\n";
                    }
                    if (vehicles[f][tekAuto] - zapastov_copy[mag[0][tekAuto]][1] > 0 && zapastov_copy[mag[0][tekAuto]][1] > 0) {
                        prostoiSr += vehicles[f][tekAuto] - zapastov_copy[mag[0][tekAuto]][1];
                        std::cout << "Машина опоздала с доставкой второго товара на " << vehicles[f][tekAuto] - zapastov_copy[mag[0][tekAuto]][1] << "мин.\n";
                    }
                    std::cout << "Машина приехала в магазин в " << vehicles[f][tekAuto] << "мин. \n" << "Текущий простой " << prostoiSr << " мин.\n";
                    // новая поставка через 8-12 часов
                    float postavka = vehicles[f][tekAuto] + ravnrasp(8 * 60, 12 * 60);
                    zapastov[mag[0][tekAuto]][0] = postavka;
                    zapastov[mag[0][tekAuto]][1] = postavka;
                    std::cout << "Следующая поставка товара 0: " << zapastov[mag[0][tekAuto]][0] << "\n";
                    std::cout << "Следующая поставка товара 1: " << zapastov[mag[0][tekAuto]][1] << "\n";
                    vehicles[f][tekAuto] += exponrasp(20);  // движение
                    vehicles[f][tekAuto] += generate_polomka(); // если сломалась то + ремонт
                    std::cout << "Машина вернулась на хлебозавод в " << vehicles[f][tekAuto] << " мин.\n";
                    vehicles[0][tekAuto] = vehicles[f][tekAuto];
                    vehicles[f][tekAuto] = infin;
                    std::cout << "------------------------------------- \n";

                }
            }
            if (fl == 1) { std::cout << "Погрузка переносится на следующую неделю \n------------------------------------------------------------ \n"; break; }
            i++;
        }
        delete[] tovar; delete[] type_p;
        for (int sd = 0; sd < 2; sd++) {
            delete[] vehicles[sd];
            delete[] mag[sd];
        }
        prSt[nCar - 1] = prostoiSr;
        prM[nCar - 1] = prostoiM;
    }
    std::cout.rdbuf(stream_buffer_cout);
    for (int d = 1; d < count - 1; d++) {
        std::cout << "Кол-во машин: " << d + 1 << ", простой магазинов: " << prSt[d] << ", простой мест на складе: " << prM[d] << "\n";
    }
    float** rez = new float* [3];
    for (int x = 0; x < 3; x++)
        rez[x] = new float[count - 1];
    for (int b = 0; b < count - 2; b++) {
        rez[0][b] = b + 2;
        rez[1][b] = prSt[b + 1];
        rez[2][b] = prM[b + 1];
    }
    bool v[3] = { false, false, false };
    float priority[3] = { 0.6, 0.3, 0.1 };
    analiz(rez, count - 2, 3, v, priority);
    file.close();
    for (int sd = 0; sd < 3; sd++) {
        delete[] rez[sd];
    }
    return 0;
}
