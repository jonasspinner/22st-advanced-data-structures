#include <iostream>
#include <chrono>

#include "io.h"
#include "commandline.h"

#include "DynamicBitVector.h"

namespace ads {
void bv_io() {
    auto [bits, operations] = ads::read_bv_input("../data/bv_test_01.txt");
    for (auto b: bits) {
        std::cout << b;
    }
    std::cout << '\n';
    for (auto operation: operations) {
        std::cout << operation << '\n';
    }
}

void bv(const std::string &input_file, const std::string &output_file) {
    auto [bits, operations] = ads::read_bv_input(input_file);
    std::ofstream output(output_file);
    if (!output) {
        throw std::runtime_error("");
    }

    auto t0 = std::chrono::steady_clock::now();

    ads::DynamicBitVector<> bv(bits);

    for (auto operation : operations) {
        switch (operation.kind) {
            case bv_operation_kind::insert:
                bv.insert(operation.i, operation.b);
                break;
            case bv_operation_kind::remove:
                bv.remove(operation.i);
                break;
            case bv_operation_kind::flip:
                bv.flip(operation.i);
                break;
            case bv_operation_kind::rank: {
                auto rank = bv.rank(operation.i, operation.b);
                output << rank << "\n";
                break;
            }
            case bv_operation_kind::select: {
                auto idx = bv.select(operation.i, operation.b);
                output << idx << "\n";
                break;
            }
        }
    }

    auto t1 = std::chrono::steady_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
    std::cout << "RESULT "
    << "algo=bv "
    << "name=<jonas_spinner> "
    << "time=<" << time.count() << "> "
    << "space=<" << bv.required_bits_upperbound() << ">\n";
}

}

int main(int argc, char *argv[]) {
    auto args = ads::parse(argc, argv);
    if (args.algo == ads::ParseResult::AlgoKind::bv) {
        ads::bv(args.input, args.output);
    }
    return 0;
}
