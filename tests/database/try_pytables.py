import tables
import numpy

class Result(tables.IsDescription):
    subset_id = tables.StringCol(32)
    model = tables.StringCol(32)
    # model = tables.FloaInt32Col()
    # generation = tables.Int32Col()

    # self.individuals.cols.id.createIndex(kind='light')
    # self.individuals.cols.parent.createIndex(kind='light')
    # self.individuals.cols.generation.createIndex(kind='light')


