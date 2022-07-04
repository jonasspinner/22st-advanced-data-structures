#include <iostream>
#include <chrono>

#include "io.h"
#include "commandline.h"

#include "DynamicBitVector.h"
#include "BP.h"

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

        ads::DynamicBitVector<uint64_t, 128> bv(bits);

        auto t1 = std::chrono::steady_clock::now();

        for (auto operation: operations) {
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

        auto t2 = std::chrono::steady_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        std::cout << "RESULT"
                  << " algo=bv"
                  << " name=jonas_spinner"
                  << " constructor=" << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()
                  << " time=" << time.count()
                  << " space=" << bv.required_bits_upperbound() << "\n";
    }

    void bp(const std::string &input_file, const std::string &output_file) {
        auto operations = ads::read_bp_input(input_file);
        std::ofstream output(output_file);
        if (!output) {
            throw std::runtime_error("");
        }

        auto t0 = std::chrono::steady_clock::now();

        ads::BP<uint64_t, 64> bp;

        int iterations = 0;
        for (auto operation: operations) {
            if (iterations % 1000 == 0) std::cout << iterations << "\n";
            switch (operation.kind) {
                case bp_operation_kind::insertchild:
                    bp.insertchild(operation.v, operation.i, operation.k);
                    break;
                case bp_operation_kind::deletenode:
                    bp.deletenode(operation.v);
                    break;
            }
            ++iterations;
        }

        auto t1 = std::chrono::steady_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
        std::cout << "RESULT"
                  << " algo=bp"
                  << " name=jonas_spinner"
                  << " time=" << time.count()
                  << " space=" << bp.required_bits_upperbound() << "\n";
    }

}

int main(int argc, char *argv[]) {
    auto args = ads::parse(argc, argv);
    if (args.algo == ads::ParseResult::AlgoKind::bv) {
        ads::bv(args.input, args.output);
    } else if (args.algo == ads::ParseResult::AlgoKind::bp) {
        ads::bp(args.input, args.output);
    }
    return 0;
}
