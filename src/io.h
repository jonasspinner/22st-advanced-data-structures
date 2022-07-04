#ifndef ADVANCED_DATA_STRUCTURES_IO_H
#define ADVANCED_DATA_STRUCTURES_IO_H

#include <iostream>
#include <fstream>
#include <vector>

namespace ads {
    enum class bv_operation_kind {
        insert,
        remove,
        flip,
        rank,
        select,
    };
    struct bv_operation {
        bv_operation_kind kind{};
        int i{};
        bool b{};
    };
}

std::ostream &operator<<(std::ostream &os, ads::bv_operation_kind kind) {
    using ads::bv_operation_kind;
    if (kind == bv_operation_kind::insert) {
        os << "insert";
    } else if (kind == bv_operation_kind::remove) {
        os << "delete";
    } else if (kind == bv_operation_kind::flip) {
        os << "flip";
    } else if (kind == bv_operation_kind::rank) {
        os << "rank";
    } else if (kind == bv_operation_kind::select) {
        os << "select";
    } else {
        os.setstate(std::ios_base::failbit);
    }
    return os;
}

std::istream &operator>>(std::istream &is, ads::bv_operation_kind &kind) {
    using ads::bv_operation_kind;
    std::string s;
    is >> s;
    if (s == "insert") {
        kind = ads::bv_operation_kind::insert;
    } else if (s == "delete") {
        kind = ads::bv_operation_kind::remove;
    } else if (s == "flip") {
        kind = ads::bv_operation_kind::flip;
    } else if (s == "rank") {
        kind = ads::bv_operation_kind::rank;
    } else if (s == "select") {
        kind = ads::bv_operation_kind::select;
    } else {
        is.setstate(std::ios_base::failbit);
    }
    return is;
}

std::ostream &operator<<(std::ostream &os, ads::bv_operation operation) {
    using ads::bv_operation_kind;
    os << operation.kind << " ";
    switch (operation.kind) {
        case bv_operation_kind::insert:
            os << operation.i << " " << operation.b;
            break;
        case bv_operation_kind::remove:
        case bv_operation_kind::flip:
            os << operation.i;
            break;
        case bv_operation_kind::rank:
        case bv_operation_kind::select:
            os << operation.b << " " << operation.i;
            break;
        default:
            os.setstate(std::ios_base::failbit);
    }
    return os;
}

std::istream &operator>>(std::istream &is, ads::bv_operation &operation) {
    using ads::bv_operation_kind;
    is >> operation.kind;
    switch (operation.kind) {
        case bv_operation_kind::insert:
            is >> operation.i >> operation.b;
            break;
        case bv_operation_kind::remove:
        case bv_operation_kind::flip:
            is >> operation.i;
            break;
        case bv_operation_kind::rank:
        case bv_operation_kind::select:
            is >> operation.b >> operation.i;
            break;
        default:
            is.setstate(std::ios_base::failbit);
            break;
    }
    return is;
}


namespace ads {
    std::pair<std::vector<bool>, std::vector<bv_operation>> read_bv_input(const std::string &input_file_name) {
        std::ifstream input_file(input_file_name);

        if (!input_file.is_open()) {
            throw std::runtime_error("file not found");
        }

        size_t n{0};
        input_file >> n;

        std::vector<bool> bits;
        bits.reserve(n);
        for (size_t i = 0; i < n; ++i) {
            bool b;
            input_file >> b;
            bits.push_back(b);
        }

        std::vector<bv_operation> operations;
        bv_operation operation{};
        while (input_file >> operation) {
            operations.push_back(operation);
        }

        if (!input_file.eof() || input_file.bad()) {
            throw std::runtime_error("an error occurred while reading bv input file");
        }
        return {bits, operations};
    }
}


namespace ads {
    enum class bp_operation_kind {
        insertchild,
        deletenode,
    };
    struct bp_operation {
        bp_operation_kind kind{};
        int v{};
        int i{};
        int k{};
    };
}

std::ostream &operator<<(std::ostream &os, ads::bp_operation_kind kind) {
    using ads::bp_operation_kind;
    if (kind == bp_operation_kind::insertchild) {
        os << "insertchild";
    } else if (kind == bp_operation_kind::deletenode) {
        os << "deletenode";
    } else {
        os.setstate(std::ios_base::failbit);
    }
    return os;
}

std::istream &operator>>(std::istream &is, ads::bp_operation_kind &kind) {
    using ads::bp_operation_kind;
    std::string s;
    is >> s;
    if (s == "insertchild") {
        kind = ads::bp_operation_kind::insertchild;
    } else if (s == "deletenode") {
        kind = ads::bp_operation_kind::deletenode;
    } else {
        is.setstate(std::ios_base::failbit);
    }
    return is;
}

std::ostream &operator<<(std::ostream &os, ads::bp_operation operation) {
    using ads::bp_operation_kind;
    os << operation.kind << " " << operation.v;
    switch (operation.kind) {
        case bp_operation_kind::insertchild:
            os << ' ' << operation.i << ' ' << operation.k;
            break;
        case bp_operation_kind::deletenode:
            break;
        default:
            os.setstate(std::ios_base::failbit);
    }
    return os;
}

std::istream &operator>>(std::istream &is, ads::bp_operation &operation) {
    using ads::bp_operation_kind;
    is >> operation.kind >> operation.v;
    switch (operation.kind) {
        case bp_operation_kind::insertchild:
            is >> operation.i >> operation.k;
            break;
        case bp_operation_kind::deletenode:
            break;
        default:
            is.setstate(std::ios_base::failbit);
            break;
    }
    return is;
}


namespace ads {
    std::vector<bp_operation> read_bp_input(const std::string &input_file_name) {
        std::ifstream input_file(input_file_name);

        if (!input_file.is_open()) {
            throw std::runtime_error("file not found");
        }

        std::vector<bp_operation> operations;
        bp_operation operation{};
        while (input_file >> operation) {
            operations.push_back(operation);
        }

        if (!input_file.eof() || input_file.bad()) {
            throw std::runtime_error("an error occurred while reading bp input file");
        }
        return operations;
    }
}

#endif //ADVANCED_DATA_STRUCTURES_IO_H
