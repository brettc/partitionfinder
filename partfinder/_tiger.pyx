# distutils: include_dirs = NUMPY_PATH
# cython: wraparound=False
# cython: cdivision=True
# cython: boundscheck=False
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

    def calc_rates(self):
        if self.species_count == 0 or self.column_count == 0:
            return None

        rates = numpy.zeros(self.column_count, dtype='f8')
        cdef: 
            size_t i, j, i_b, j_b
            double rate, axpi, num, denom
            c_Bitset *i_bitset
            c_Bitset *j_bitset
            np.npy_double[:] c_rates = rates

        denom = <double>self.column_count - 1.0
        dist_matrix = numpy.empty([self.column_count, self.column_count], dtype=float)
        for i in range(self.column_count):
            rate = 0.0
            for j in range(self.column_count):
                # Don't compare to self
                if i == j:
                    dist_matrix[i][j] = 1

                else:
                    num = 0.0
                    axpi = 0.0
                    for j_b in range(4):
                        j_bitset = &self._bitsets[j][j_b]
                        if j_bitset.none():
                            continue
                        num += 1.0
                        for i_b in range(4):
                            i_bitset = &self._bitsets[i][i_b]
                            if i_bitset.none():
                                continue
                            if j_bitset.is_subset_of(deref(i_bitset)):
                                axpi += 1.0
                                break

                    rate += axpi / num
                    dist_matrix[i][j] = (axpi / num)

            rate /= denom
            c_rates[i] = rate

        return rates, dist_matrix

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

cdef class TigerDNA(TigerBase):
    def __cinit__(self): 
        self.species_count = 0
        self.column_count = 0

    def build_bitsets(self, alignment):
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


