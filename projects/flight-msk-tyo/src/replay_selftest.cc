#include <random>
#include <vector>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>
#include <numeric>
#include <set>

// Предполагаем, что эти заголовки доступны
#include "replay.hh"

using namespace std;

double unionMeasure(const vector<double>& ys, double ell, double L) {
    if (ys.empty()) return 0.0;
    vector<pair<double, double>> intervals;
    double halfL = (L - ell) / 2.0;   // допустимые положения прибора: окно целиком в трубке
    double halfEll = ell / 2.0;
    for (double y : ys) {
        double left = y - halfEll;
        double right = y + halfEll;
        if (right < -halfL || left > halfL) continue;
        intervals.emplace_back(max(left, -halfL), min(right, halfL));
    }
    if (intervals.empty()) return 0.0;
    sort(intervals.begin(), intervals.end());
    double total = 0.0;
    double curLeft = intervals[0].first;
    double curRight = intervals[0].second;
    for (size_t i = 1; i < intervals.size(); ++i) {
        if (intervals[i].first <= curRight) {
            curRight = max(curRight, intervals[i].second);
        } else {
            total += curRight - curLeft;
            curLeft = intervals[i].first;
            curRight = intervals[i].second;
        }
    }
    total += curRight - curLeft;
    return total / (L - ell);
}

bool test1() {
    const int N = 4000000;
    vector<vector<double>> events = {
        {0.0},
        {0.0, 10.0, 100.0},
        {140.0, 145.0},
        {-149.0, -100.0, 0.0, 149.0}
    };
    vector<double> expected = {
        unionMeasure({0.0}, 30.0, 300.0),   // ручное значение 0.1
        unionMeasure({0.0, 10.0, 100.0}, 30.0, 300.0),
        unionMeasure({140.0, 145.0}, 30.0, 300.0),
        unionMeasure({-149.0, -100.0, 0.0, 149.0}, 30.0, 300.0)
    };
    vector<string> labels = {"a", "b", "c", "d"};
    bool allPass = true;
    for (size_t i = 0; i < events.size(); ++i) {
        mt19937_64 gen(12345 + static_cast<int>(i));
        uniform_real_distribution<double> dist(0.0, 1.0);
        auto rand01 = [&]() { return dist(gen); };
        vector<PspRec> ev(events[i].size());
        for (size_t j = 0; j < events[i].size(); ++j) {
            ev[j].y_cm = events[i][j];
        }
        double sumW = 0.0;
        double sumW2 = 0.0;
        vector<int> idx;
        double w, zd;
        for (int k = 0; k < N; ++k) {
            if (ChooseGroup(ev, 30.0, 300.0, rand01, idx, w, zd)) {
                int cnt = 0;   // независимый счёт записей в окне ±15 см вокруг выбранной точки
                for (const auto& r : ev) if (fabs(r.y_cm - zd) <= 15.0) ++cnt;
                if (cnt != static_cast<int>(idx.size())) allPass = false;
                sumW += w;
                sumW2 += w * w;
            }
        }
        double mean = sumW / N;
        double var = max(0.0, sumW2 / N - mean * mean);
        double se = sqrt(var / N);
        double diff = abs(mean - expected[i]);
        bool pass = (diff <= 5.0 * se + 1e-9) && (diff <= 0.01 * expected[i]);
        printf("T1 %s expected=%.6f mean=%.6f se=%.6f %s\n",
               labels[i].c_str(), expected[i], mean, se, pass ? "PASS" : "FAIL");
        if (!pass) allPass = false;
    }
    return allPass;
}

bool test2() {
    // Записываем файлы
    {
        PspWriter w;
        if (!w.Open("rs_A.psp")) { printf("cannot open rs_A.psp\n"); return false; }
        vector<PspRec> recs = {{}, {}, {}, {}, {}, {}};
        recs[0].evt = 5; recs[1].evt = 5;
        recs[2].evt = 6;
        recs[3].evt = 7; recs[4].evt = 7; recs[5].evt = 7;
        for (auto& r : recs) w.Write(r);
        w.Close();
    }
    {
        PspWriter w;
        if (!w.Open("rs_B.psp")) { printf("cannot open rs_B.psp\n"); return false; }
        w.Close();
    }
    {
        PspWriter w;
        if (!w.Open("rs_C.psp")) { printf("cannot open rs_C.psp\n"); return false; }
        vector<PspRec> recs = {{}, {}, {}};
        recs[0].evt = 7; recs[1].evt = 7;
        recs[2].evt = 8;
        for (auto& r : recs) w.Write(r);
        w.Close();
    }

    PspEvents ev;
    bool a = ev.AddFile("rs_A.psp");
    bool b = ev.AddFile("rs_missing.psp");
    bool c = ev.AddFile("rs_B.psp");
    bool d = ev.AddFile("rs_C.psp");

    vector<int> sizes;
    vector<PspRec> buf;
    vector<int> firstEvt;   // evt записей ПЕРВОГО события
    while (ev.NextEvent(buf)) {
        if (sizes.empty()) for (auto& r : buf) firstEvt.push_back(r.evt);
        sizes.push_back(static_cast<int>(buf.size()));
    }
    long long count = ev.EventsRead();

    // Удаляем временные файлы
    remove("rs_A.psp");
    remove("rs_B.psp");
    remove("rs_C.psp");

    vector<int> expected = {2, 1, 3, 2, 1};
    bool pass = (sizes == expected) && (count == 5);
    if (firstEvt != vector<int>{5, 5}) pass = false;
    printf("T2 sizes=");
    for (size_t i = 0; i < sizes.size(); ++i) {
        printf("%d", sizes[i]);
        if (i != sizes.size() - 1) printf(",");
    }
    printf(" expected=2,1,3,2,1 %s\n", pass ? "PASS" : "FAIL");
    return pass;
}

int main() {
    try {
        bool t1 = test1();
        bool t2 = test2();
        if (t1 && t2) {
            cout << "SELFTEST PASS" << endl;
            return 0;
        } else {
            cout << "SELFTEST FAIL" << endl;
            return 1;
        }
    } catch (const exception& e) {
        cerr << "Exception: " << e.what() << endl;
        return 2;
    }
}
