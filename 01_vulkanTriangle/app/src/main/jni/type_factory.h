#pragma once

#include "type.h"
#include <map>
#include <unordered_map>
#include <mutex>

class TypeFactory {
public:
    static TypeFactory &get_instance();

    // TODO(type): maybe it makes sense to let each get_X function return X*
    // instead of generic Type*

    Type *get_primitive_type(PrimitiveTypeID id);

    PrimitiveType *get_primitive_int_type(int bits, bool is_signed = true);

    Type *get_vector_type(int num_elements, Type *element);

    Type *get_pointer_type(Type *element, bool is_bit_pointer = false);

    Type *get_custom_int_type(int num_bits, bool is_signed, Type *compute_type);

    Type *get_custom_float_type(Type *digits_type,
                                Type *exponent_type,
                                Type *compute_type,
                                float64 scale);

    Type *get_bit_struct_type(PrimitiveType *physical_type,
                              const std::vector<Type *> &member_types,
                              std::vector<int> member_bit_offsets);

    Type *get_bit_array_type(PrimitiveType *physical_type,
                             Type *element_type,
                             int num_elements);

    static DataType create_vector_or_scalar_type(int width,
                                                 DataType element,
                                                 bool element_is_pointer = false);

private:
    TypeFactory();

    std::unordered_map<PrimitiveTypeID, std::unique_ptr<Type>> primitive_types_;

    // TODO: use unordered map
    std::map<std::pair<int, Type *>, std::unique_ptr<Type>> vector_types_;

    // TODO: is_bit_ptr?
    std::map<std::pair<Type *, bool>, std::unique_ptr<Type>> pointer_types_;

    // TODO: use unordered map
    std::map<std::tuple<int, bool, Type *>, std::unique_ptr<Type>>
            custom_int_types;

    // TODO: use unordered map
    std::map<std::tuple<Type *, Type *, Type *, float64>, std::unique_ptr<Type>>
            custom_float_types;

    // TODO: avoid duplication
    std::vector<std::unique_ptr<Type>> bit_struct_types_;

    // TODO: avoid duplication
    std::vector<std::unique_ptr<Type>> bit_array_types_;

    std::mutex mut_;
};

DataType promoted_type(const DataType &a, const DataType &b);
