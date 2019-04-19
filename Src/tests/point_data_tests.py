from point_data import PointData
from point_data import PointDataList


def point_data_tests():
    pd = PointData((100.321, 100.125315), 4)
    pd.add((100, 100), 1.111111)
    pd.add((100, 100), 1.32312)
    pd.add((100, 100), 1231.3122)
    pd.add((100, 100), 1123.312)
    pd.add((100, 100), 1321.1233)
    print(pd.to_string())


def point_data_list_tests():
    pdl = PointDataList(4)
    pdl.set((100, 100), 10, 3)
    pdl.set((100, 100), 10, 2)
    pdl.set((82, 82), 123, 0)
    pdl.set((101, 99), 1234, 10)
    print(pdl.to_string())
    print(pdl.to_list())


point_data_list_tests()
