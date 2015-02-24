# distutils: include_dirs = NUMPY_PATH
# cython: wraparound=False
# cython: cdivision=True
# xxxxcython: boundscheck=False
#
# from _utility cimport dynamic_bitset
cimport numpy as np

import cython
import numpy
import copy
import sys
import numpy
from cython.operator import dereference as deref, preincrement as preinc

cdef class TigerBase:
    def __cinit__(self):
        pass


cdef class TigerDNA(TigerBase):
    def __cinit__(self): 
        self.species_count = 0
        self.column_count = 0

    def from_alignment(self, alignment):
        cdef:
            size_t i, j, k, sp_count, col_count
            np.npy_uint8[:, :] data = alignment.data
            unsigned char c
            c_Bitset *A_bits
            c_Bitset *C_bits
            c_Bitset *G_bits
            c_Bitset *T_bits

        sp_count, col_count = alignment.data.shape
        self.species_count = sp_count
        self.column_count = col_count

        # Create all the bitsets
        # NOTE: would be better if we could rely on emplace_back...
        for i in range(col_count):
            self._bitsets.push_back(c_ColumnBitsets())
            for k in range(4):
                self._bitsets.back().push_back(c_Bitset(sp_count))

        for i in range(col_count):
            A_bits = &self._bitsets[i][0]
            C_bits = &self._bitsets[i][1]
            G_bits = &self._bitsets[i][2]
            T_bits = &self._bitsets[i][3]
            for j in range(sp_count):
                c = data[j, i]
                if c == 'A':
                    A_bits.set(j)
                elif c == 'C':
                    C_bits.set(j)
                elif c == 'G':
                    G_bits.set(j)
                elif c == 'T':
                    T_bits.set(j)
                elif c == 'N' or c == '?' or c == '-':
                    A_bits.set(j)
                    C_bits.set(j)
                    G_bits.set(j)
                    T_bits.set(j)


    def calc_rates(self):
        if self.species_count == 0 or self.column_count == 0:
            return None

        cdef: 
            size_t i, j, i_b, j_b
            double rate, axpi, num
            double denom = <double>self.column_count
            c_Bitset *i_bitset
            c_Bitset *j_bitset
            vector[double] rates

        for i in range(self.column_count):
            rate = 0.0
            for j in range(self.column_count):
                # Don't compare to self
                if i == j:
                    continue

                num = 0.0
                axpi = 0.0
                for j_b in range(4):
                    j_bitset = &self._bitsets[j][j_b]
                    if j_bitset.empty():
                        continue
                    num += 1.0
                    for i_b in range(4):
                        if i_bitset.empty():
                            continue
                        i_bitset = &self._bitsets[i][i_b]
                        if j_bitset.is_subset_of(deref(i_bitset)):
                            axpi += 1.0
                            break

                rate += axpi / num

            rate /= denom
            rates.push_back(rate)

        return rates


    def bitsets_as_array(self):
        if self.species_count == 0 or self.column_count == 0:
            return None

        ret = numpy.zeros((self.column_count, 4, self.species_count), dtype='u1')
        cdef: 
            np.npy_uint8[:, :, :] c_ret = ret
            size_t i, j, k

        for i in range(self.column_count):
            for j in range(4):
                for k in range(self.species_count):
                    if self._bitsets[i][j].test(k):
                        c_ret[i, j, k] = 1
        return ret



# cdef class Bitset:
#     def __cinit__(self, size_t size):
#         self._this.resize(size)
#
#     def test(self, size_t i):
#         return self._this.test(i)
#
#     def __getitem__(self, size_t i):
#         return self._this.test(i)
#
#     property size:
#         def __get__(self):
#             return self._this.size()
#
#     def as_array(self):
#         vals = numpy.zeros(self.size, dtype='u1')
#         cdef: 
#             np.npy_uint8[:] v = vals
#             size_t i
#
#         for i in range(self._this.size()):
#             v[i] = self._this.test(i)
#         return vals
#
#     def set(self, size_t i):
#         self._this.set(i)
#
#     def reset(self, size_t i):
#         self._this.reset(i)
#
#     def flip(self, size_t i):
#         self._this.flip(i)
#
#     def __setitem__(self, size_t i, bint b):
#         if b:
#             self._this.set(i)
#         else:
#             self._this.reset(i)
#
#     def empty(self):
#         return self._this.empty()
#
#     def is_subset_of(self, Bitset other):
#         return self._this.is_subset_of(other._this)

    # def __cmp__(self, ChannelStateFrozen other):
    #     return bitset_cmp(self._this, other._this)

    # def __repr__(self):
    #     return "<ChannelsRO: {}>".format(self.__str__())

    # def __copy__(self):
    #     other = ChannelState()
    #     other.init(self.world, self._this)
    #     return other
    #
    # def copy(self):
    #     return self.__copy__()

    # def merge(self, ChannelStateFrozen other):
    #     self._this |= other._this

