#include <iostream>
#include <fstream>

#include "Matrix.h"
#include "TfRunner.h"
#include "FuzzyController.h"

// Simulate open-loop system
void runOpenLoop(TfRunner tf, char* argv0)
{
    const float Te = 0.1F; // [s]
    const float runTime = 1000.F; // [s]
    const int cycles = (int)ceil(runTime / Te);

    // open-loop response
    std::ofstream f("y.txt");
    for (int i = 0; i < cycles; ++i)
    {
        f << tf.step(1.F) << std::endl;
    }

    // script for open-loop plot
    std::string path(argv0);
    path = path.substr(0U, path.find_last_of('\\') + 1U);
    std::ofstream script("plot_open_loop.m");

    script << "load('" << path.c_str() << "y.txt');" << std::endl;
    script << "t=0:0.1:999.9;" << std::endl;
    script << "figure" << std::endl;
    script << "hold on" << std::endl;
    script << "grid on" << std::endl;
    script << "plot(t, y)" << std::endl;
}

// Simulate closed-loop system
// with Fuzzy controller
void runClosedLoop(TfRunner tf, char* argv0)
{
    // vars id
    const int e  = 0;
    const int de = 1;
    const int u  = 2;

    // MFs id
    const int NB = 0;
    const int NS = 1;
    const int ZE = 2;
    const int PS = 3;
    const int PB = 4;

    // var ranges
    //                   min, max, step
    const Range e_range(-1.1F, 1.1F, 0.01F);
    const Range de_range(-1.1F, 1.1F, 0.01F);
    const Range u_range(-1.1F, 1.1F, 0.01F);

    // MFs
    MFVector mfs;  // type,    a,     b,    c,    d
    mfs.emplace_back(MF_TRI, -1.5F, -1.F, -0.5F, 0.F); // NB
    mfs.emplace_back(MF_TRI, -1.F,  -0.5F, 0.F,  0.F); // NS
    mfs.emplace_back(MF_TRI, -0.5F,  0.F,  0.5F, 0.F); // ZE
    mfs.emplace_back(MF_TRI,  0.F,   0.5F, 1.F,  0.F); // PS
    mfs.emplace_back(MF_TRI,  0.5F,  1.F,  1.5F, 0.F); // PB

    RuleVector rules = {
        //de is { NB, NS, ZE, PS, PB }
                { NB, NB, NS, NS, NS }, // e is NB
                { NS, NS, NS, NS, ZE }, // e is NS
                { ZE, ZE, ZE, ZE, ZE }, // e is ZE
                { ZE, PS, PS, PS, PS }, // e is PS
                { PS, PS, PS, PB, PB }  // e is PB
    };

    FuzzyController controller;
    controller.setRange(e, e_range);
    controller.setRange(de, de_range);
    controller.setRange(u, u_range);

    controller.setMFs(e, mfs);
    controller.setMFs(de, mfs);
    controller.setMFs(u, mfs);

    controller.setRules(rules);

    controller.build();

    std::ofstream f("e_de_u.txt");
    for (int i = 10; i < (e_range.size() - 10); ++i)
    {
        for (int j = 10; j < (de_range.size() - 10); ++j)
        {
            f << e_range.at(i) << '\t';
            f << de_range.at(j) << '\t';
            f << controller.calculate(e_range.at(i), de_range.at(j));
            f << std::endl;
        }
    }

    const float set_point = 1.F;
    float err = 1.F;
    float derr = 0.F;

    const float Te = 0.1F; // [s]
    const float runTime = 1000.F; // [s]
    const int cycles = (int)ceil(runTime / Te);

    std::ofstream f2("y2.txt");
    for (int i = 0; i < cycles; ++i)
    {
        const float u = controller.calculate(err, derr);
        const float y = tf.step(u);

        const float e = set_point - y;
        derr = (e - err) / Te;
        err = e;

        f2 << y << std::endl;
    }
}

// Manually precalculated result!
// The result should be: ~8.39
void test(TfRunner tf, char* argv0)
{
    // vars id
    const int e = 0;
    const int de = 1;
    const int u = 2;

    // MFs id
    const int NB = 0;
    const int NM = 1;
    const int NS = 2;
    const int ZE = 3;
    const int PS = 4;
    const int PM = 5;
    const int PB = 6;

    // var ranges
    //                   min, max, step
    const Range e_range(-18.F, 18.F, 0.1F);
    const Range de_range(-13.F, 13.F, 0.1F);
    const Range u_range(-19.F, 19.F, 0.1F);

    FuzzyController controller;

    // e MFs
    MFVector mfs;  // type,    a,     b,    c,    d
    mfs.emplace_back(MF_TRAP, -20.F, -20.F, -18.F, -16.F); // NB
    mfs.emplace_back(MF_TRAP, -18.F, -16.F, -8.5F, -5.5F); // NM
    mfs.emplace_back(MF_TRAP, -8.5F, -5.5F, -1.F, 0.F); // NS
    mfs.emplace_back(MF_TRI, -1.F, 0.F, 1.F, 0.F); // ZE
    mfs.emplace_back(MF_TRAP, 0.F, 1.F, 5.5F, 8.5F); // PS
    mfs.emplace_back(MF_TRAP, 5.5F, 8.5F, 16.F, 18.F); // PM
    mfs.emplace_back(MF_TRAP, 16.F, 18.F, 20.F, 20.F); // PB
    controller.setMFs(e, mfs);

    mfs.clear();
    mfs.emplace_back(MF_TRAP, -15.F, -15.F, -13.F, -11.F); // NB
    mfs.emplace_back(MF_TRAP, -13.F, -11.F, -6.F, -4.F); // NM
    mfs.emplace_back(MF_TRAP, -6.F, -4.F, -0.5F, 0.F); // NS
    mfs.emplace_back(MF_TRI, -0.5F, 0.F, 0.5F, 0.F); // ZE
    mfs.emplace_back(MF_TRAP, 0.F, 0.5F, 4.F, 6.F); // PS
    mfs.emplace_back(MF_TRAP, 4.F, 6.F, 11.F, 13.F); // PM
    mfs.emplace_back(MF_TRAP, 11.F, 13.F, 15.F, 15.F); // PB
    controller.setMFs(de, mfs);

    mfs.clear();
    mfs.emplace_back(MF_TRAP, -20.F, -20.F, -19.F, -17.F); // NB
    mfs.emplace_back(MF_TRAP, -19.F, -17.F, -8.F, -7.F); // NM
    mfs.emplace_back(MF_TRAP, -8.F, -7.F, -1.F, 0.F); // NS
    mfs.emplace_back(MF_TRI, -1.F, 0.F, 1.F, 0.F); // ZE
    mfs.emplace_back(MF_TRAP, 0.F, 1.F, 7.F, 8.F); // PS
    mfs.emplace_back(MF_TRAP, 7.F, 8.F, 17.F, 19.F); // PM
    mfs.emplace_back(MF_TRAP, 17.F, 19.F, 20.F, 20.F); // PB
    controller.setMFs(u, mfs);

    RuleVector rules = {
        //de is { NB, NM, NS, ZE, PS, PM, PB }
                { PB, PM, PS, ZE, NS, NS, NS }, // e is NB
                { PB, PM, PS, ZE, NS, NS, NS }, // e is NM
                { PB, PM, PS, ZE, NS, NS, NM }, // e is NS
                { PB, PM, PS, ZE, NS, NM, NM }, // e is ZE
                { PB, PS, PS, ZE, NS, NM, NM }, // e is PS
                { PM, PS, PS, ZE, NS, NM, NM }, // e is PM
                { PM, PS, PS, ZE, NS, NM, NM }  // e is PB
    };

    controller.setRange(e, e_range);
    controller.setRange(de, de_range);
    controller.setRange(u, u_range);

    controller.setRules(rules);

    controller.build();

    const float y = controller.calculate(-6.F, -5.F);

    // test returned value
    assert(fabsf(y - 8.39583302F) < std::numeric_limits<float>::epsilon());
}

int main(int argc, char* argv[])
{
    // State-space matrices
    Matrix A(2, 2);
    Matrix B(2, 1);
    Matrix C(1, 2);
    Matrix D(1, 1);

    A.at(0, 0) = 0.9931F;
    A.at(0, 1) = -4.152e-05F;
    A.at(1, 0) = 0.09965F;
    A.at(1, 1) = 1.F;

    B.at(0, 0) = 0.09965F;
    B.at(1, 0) = 0.004988F;

    C.at(0, 0) = 0.F;
    C.at(0, 1) = 0.000625F;

    D.at(0, 0) = 0.F;

    TfRunner tf(A, B, C, D);

    runOpenLoop(tf, argv[0]);
    runClosedLoop(tf, argv[0]);
    //test(tf, argv[0]);

    return 0;
}
