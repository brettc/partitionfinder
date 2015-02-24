"""
Definitions for c++ or boost classes
"""
from libcpp.vector cimport vector
# from libcpp.string cimport string
# from libcpp.utility cimport pair
# from libcpp.cast import dynamic_cast, static_cast
# from libcpp.set cimport set as std_set
# from libcpp.map cimport map as std_map

cdef extern from "<boost/dynamic_bitset.hpp>" namespace "boost":
    cdef cppclass dynamic_bitset[T]:
        dynamic_bitset(size_t)
        void resize(size_t)
        void set(size_t)
        void reset(size_t)
        void flip(size_t)

        size_t size()
        bint test(size_t)
        bint empty()
        bint is_subset_of(dynamic_bitset[T]&)
        bint is_proper_subset_of(dynamic_bitset[T]&)
        bint intersects(dynamic_bitset[T]& a)

    # cdef void to_string(dynamic_bitset[size_t], string s)

ctypedef dynamic_bitset[size_t] c_Bitset
ctypedef vector[c_Bitset] c_ColumnBitsets
ctypedef vector[c_ColumnBitsets] c_AlignmentBitsets
    
# cdef class Bitset:
#     cdef c_Bitset _this
#
# cdef class ColumnBitsets:
#     cdef c_ColumnBitsets _this
#
cdef class TigerBase:
    cdef:
        c_AlignmentBitsets _bitsets
        readonly: 
            size_t species_count, column_count
